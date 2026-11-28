// Multithreading performance benchmark for GenFit.
//
// Times <totalFits> Kalman fits over deterministic toy events (muons crossing the
// test geometry from test/makeGeom.C, 5-14 measurements each) and prints the wall
// time, fit throughput and a result checksum. Measurement generation is part of
// the timed region. Run it from the directory containing genfitGeom.root.
//
//   mtBenchmark <threads> <totalFits> [tgeo|analytic|noEffects]
//
// Material modes:
//   tgeo     - genfit::TGeoMaterialInterface (default; ROOT's TGeo navigation is
//              not concurrency safe, so its usage is serialized)
//   analytic - ShellMaterialInterface below: a thread-safe, non-TGeo material
//              model of the same geometry, so nothing is serialized
//   noEffects - material effects disabled per worker thread (upper scaling bound)
//
// The events and their measurements are generated deterministically per event
// index, so the checksum is identical for every thread count within a mode.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

#include <TDatabasePDG.h>
#include <TGeoManager.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TRandom3.h>

#include <FieldManager.h>
#include <ConstField.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <MeasurementCreator.h>
#include <MathHelpers.h>
#include <HelixTrackModel.h>
#include <RKTrackRep.h>
#include <KalmanFitterRefTrack.h>
#include <Track.h>
#include <TrackPoint.h>
#include <Exception.h>

namespace {

// Boundaries of the test geometry (test/makeGeom.C): four silicon tubes
// [R[2i], R[2i+1]] of half length zMax, in air.
const double kR[8] = {1.0, 1.05, 2.0, 2.05, 3.0, 3.05, 4.0, 4.05};
const double kZMax = 20.;

bool inSilicon(double x, double y, double z) {
  if (fabs(z) > kZMax) return false;
  const double r2 = x * x + y * y;
  for (int i = 0; i < 4; ++i)
    if (r2 >= kR[2 * i] * kR[2 * i] && r2 <= kR[2 * i + 1] * kR[2 * i + 1]) return true;
  return false;
}

/**
 * @brief Thread-safe analytic AbsMaterialInterface for the test geometry.
 *
 * Serves the same materials as TGeoMaterialInterface for the four silicon tubes
 * in air of test/makeGeom.C, but computes boundaries analytically instead of
 * navigating ROOT's TGeoManager, whose navigation is not concurrency safe.
 * All state is thread_local, so several threads may share one instance and
 * nothing is serialized. The step search mirrors
 * TGeoMaterialInterface::findNextBoundary() (same delta/epsilon/RKPropagate
 * arch control); only the geometry queries are replaced by closed-form
 * distances to cylinders and endcap planes.
 */
class ShellMaterialInterface : public genfit::AbsMaterialInterface {
 public:
  ShellMaterialInterface(const genfit::Material& silicon, const genfit::Material& air)
      : silicon_(silicon), air_(air) {}

  bool initTrack(double posX, double posY, double posZ,
                 double dirX, double dirY, double dirZ) override {
    State& s = state();
    const bool silicon = inSilicon(posX, posY, posZ);
    const bool changed = s.valid && (silicon != s.silicon);
    s.pos[0] = posX; s.pos[1] = posY; s.pos[2] = posZ;
    s.dir[0] = dirX; s.dir[1] = dirY; s.dir[2] = dirZ;
    s.silicon = silicon;
    s.valid = true;
    return changed;
  }

  genfit::Material getMaterialParameters() override {
    return state().silicon ? silicon_ : air_;
  }

  double findNextBoundary(const genfit::RKTrackRep* rep,
                          const genfit::M1x7& stateOrig,
                          double sMax, // signed
                          bool varField) override {
    const double delta(1.E-2); // cm, distance limit beneath which straight-line steps are taken.
    const double epsilon(1.E-1); // cm, allowed upper bound on arch deviation

    genfit::M1x3 SA;
    genfit::M1x7 state7, oldState7;
    oldState7 = stateOrig;

    int stepSign(sMax < 0 ? -1 : 1);

    double s = 0;  // trajectory length to boundary

    const unsigned maxIt = 300;
    unsigned it = 0;

    // Initialize to the current location (set by the caller via initTrack).
    double safety, slDist;
    safetyStep(fabs(sMax) - s, safety, slDist);
    if (slDist < safety)
      slDist = safety;
    double step = slDist;

    while (1) {
      if (++it > maxIt) {
        genfit::Exception exc("ShellMaterialInterface::findNextBoundary ==> maximum number of iterations exceeded", __LINE__, __FILE__);
        exc.setFatal();
        throw exc;
      }

      // No boundary in sight?
      if (s + safety > fabs(sMax)) {
        return stepSign * (s + safety);
      }

      // We are at the boundary?
      if (slDist < delta) {
        return stepSign * (s + slDist);
      }

      // Follow curved arch, then see if we may have missed a boundary.
      // Always propagate complete way from original start to avoid
      // inconsistent extrapolations.
      state7 = stateOrig;
      rep->RKPropagate(state7, nullptr, SA, stepSign * (s + step), varField);

      // Straight line distance^2 between extrapolation finish and
      // the end of the previously determined safe segment.
      double dist2 = (genfit::square(state7[0] - oldState7[0])
          + genfit::square(state7[1] - oldState7[1])
          + genfit::square(state7[2] - oldState7[2]));
      // Maximal lateral deviation^2.
      double maxDeviation2 = 0.25 * (step * step - dist2);

      if (step > safety && maxDeviation2 > epsilon * epsilon) {
        // Take a shorter step to reliably estimate material, never shorter
        // than safety.
        step = std::max(step / 2, safety);
      } else {
        const State saved = state();
        bool volChanged = initTrack(state7[0], state7[1], state7[2],
            stepSign * state7[3], stepSign * state7[4], stepSign * state7[5]);

        if (volChanged) {
          state() = saved; // move back to start
          if (step <= safety) {
            return stepSign * (s + step);
          }
          step = std::max(step / 2, safety);
        } else {
          // we are in the new place, the step was safe, advance
          s += step;
          oldState7 = state7;

          safetyStep(fabs(sMax) - s, safety, slDist);
          if (slDist < safety)
            slDist = safety;
          step = slDist;
        }
      }
    }
  }

 private:
  struct State {
    double pos[3] = {0., 0., 0.};
    double dir[3] = {0., 0., 1.};
    bool silicon = false;
    bool valid = false;
  };

  // Per-thread navigator state (initTrack/findNextBoundary/getMaterialParameters
  // form one usage sequence per thread, as documented in AbsMaterialInterface).
  static State& state() {
    static thread_local State s;
    return s;
  }

  //! Lower bound of the distance to any boundary (any surface) and the distance
  //! to the nearest boundary surface along the current direction, capped at maxDist.
  void safetyStep(double maxDist, double& safety, double& step) const {
    const State& s = state();
    const double x = s.pos[0], y = s.pos[1], z = s.pos[2];
    const double dx = s.dir[0], dy = s.dir[1], dz = s.dir[2];
    const double r = sqrt(x * x + y * y);

    // safety: distance to the nearest surface (never overestimates the true
    // distance to the boundary of the current region)
    safety = 1.E30;
    for (double R : kR) safety = std::min(safety, fabs(r - R));
    safety = std::min(safety, fabs(kZMax - fabs(z)));

    // step: nearest positive surface crossing along the direction
    double nearest = 1.E30;
    const double a = dx * dx + dy * dy;
    const double b = 2. * (x * dx + y * dy);
    for (double R : kR) {
      if (a < 1.E-30) continue;
      const double c = x * x + y * y - R * R;
      const double disc = b * b - 4. * a * c;
      if (disc <= 0.) continue;
      const double sq = sqrt(disc);
      for (double t : {(-b - sq) / (2. * a), (-b + sq) / (2. * a)}) {
        if (t <= 1.E-10) continue;
        if (fabs(z + t * dz) <= kZMax) nearest = std::min(nearest, t);
      }
    }
    if (fabs(dz) > 1.E-30) {
      for (double zPlane : {-kZMax, kZMax}) {
        const double t = (zPlane - z) / dz;
        if (t <= 1.E-10) continue;
        const double rCross = sqrt(genfit::square(x + t * dx) + genfit::square(y + t * dy));
        for (int i = 0; i < 4; ++i) {
          if (rCross >= kR[2 * i] && rCross <= kR[2 * i + 1]) {
            nearest = std::min(nearest, t);
            break;
          }
        }
      }
    }
    step = std::min(nearest, maxDist);
  }

  genfit::Material silicon_;
  genfit::Material air_;
};

} // namespace

static void eventParams(int idx, TVector3& mom, unsigned& nMeas) {
  mom.SetXYZ(1., 0., 0.);
  mom.SetPhi(0.062831853 * (idx % 100));
  mom.SetTheta(0.45 * TMath::Pi() + 0.001 * (idx % 17));
  mom.SetMag(0.2 + 0.008 * (idx % 100));
  nMeas = 5 + (idx % 10);
}

static double fitEvent(int idx, double charge, genfit::MeasurementCreator& creator,
                       genfit::AbsKalmanFitter* fitter) {
  TVector3 pos(0, 0, 0), mom;
  unsigned nMeasurements;
  eventParams(idx, mom, nMeasurements);

  const int pdg = 13;
  genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
  creator.setTrackModel(helix);

  TVector3 posM(pos), momM(mom);
  posM.SetX(posM.X() + 0.01 * ((idx % 7) - 3));
  posM.SetY(posM.Y() + 0.01 * ((idx % 5) - 2));
  momM.SetPhi(momM.Phi() + 0.002 * ((idx % 3) - 1));
  momM.SetMag(momM.Mag() * (1. + 0.01 * ((idx % 4) - 2)));

  TMatrixDSym covM(6);
  const double resolution = 0.01;
  for (int i = 0; i < 3; ++i) covM(i, i) = resolution * resolution;
  for (int i = 3; i < 6; ++i) covM(i, i) = genfit::square(resolution / nMeasurements / sqrt(3.));

  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
  genfit::MeasuredStateOnPlane stateSmeared(rep);
  stateSmeared.setPosMomCov(posM, momM, covM);
  TVectorD seedState(6);
  TMatrixDSym seedCov(6);
  stateSmeared.get6DStateCov(seedState, seedCov);
  genfit::Track fitTrack(rep, seedState, seedCov);

  for (unsigned i = 0; i < nMeasurements; ++i) {
    std::vector<genfit::AbsMeasurement*> measurements =
        creator.create(genfit::EMeasurementType((idx + i) % 8), i * 5.);
    fitTrack.insertPoint(new genfit::TrackPoint(measurements, &fitTrack));
  }
  fitter->processTrack(&fitTrack);

  TVector3 fitPos, fitMom;
  fitTrack.getFittedState(0, rep).getPosMom(fitPos, fitMom);
  return fitMom.Mag();
}

int main(int argc, char** argv) {
  const int nThreads = (argc > 1) ? atoi(argv[1]) : 1;
  const int totalFits = (argc > 2) ? atoi(argv[2]) : 128;
  const std::string mode = (argc > 3) ? argv[3] : "tgeo";

  gRandom->SetSeed(42);
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., 15.));

  genfit::Material siliconMat, airMat;
  {
    // sample the materials of the test geometry once, single-threaded
    genfit::TGeoMaterialInterface probe;
    probe.initTrack(1.025, 0., 0., 1., 0., 0.);
    siliconMat = probe.getMaterialParameters();
    probe.initTrack(0., 0., 0., 1., 0., 0.);
    airMat = probe.getMaterialParameters();
  }
  if (mode == "analytic")
    genfit::MaterialEffects::getInstance()->init(new ShellMaterialInterface(siliconMat, airMat));
  else
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());

  // one creator per event, seeded per event outside the timed region: results
  // are identical for every thread count
  std::vector<genfit::MeasurementCreator*> creators(totalFits);
  for (int i = 0; i < totalFits; ++i) {
    gRandom->SetSeed(1000 + i);
    creators[i] = new genfit::MeasurementCreator();
  }

  std::vector<double> partial(nThreads, 0.); // one slot per worker: no shared writes
  const int perThread = (totalFits + nThreads - 1) / nThreads;
  // TDatabasePDG is a ROOT singleton: look it up once, on the main thread.
  const double charge = TDatabasePDG::Instance()->GetParticle(13)->Charge() / (3.);

  auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> threads;
  for (int k = 0; k < nThreads; ++k) {
    threads.emplace_back([k, perThread, totalFits, charge, mode, &partial, &creators]() {
      // setNoEffects is per-thread (per-thread MaterialEffects instances)
      if (mode == "noEffects") genfit::MaterialEffects::getInstance()->setNoEffects(true);
      genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
      double sum = 0.;
      for (int j = 0; j < perThread; ++j) {
        const int idx = k * perThread + j;
        if (idx >= totalFits) break;
        sum += fitEvent(idx, charge, *creators[idx], fitter);
      }
      partial[k] = sum;
      delete fitter;
    });
  }
  for (auto& t : threads) t.join();
  auto t1 = std::chrono::steady_clock::now();

  double checksum = 0.;
  for (double p : partial) checksum += p;

  const double wall = std::chrono::duration<double>(t1 - t0).count();
  printf("threads=%2d fits=%4d wall=%7.3f s  rate=%8.1f fits/s  checksum=%.6f\n",
         nThreads, totalFits, wall, totalFits / wall, checksum);
  return 0;
}
