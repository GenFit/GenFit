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
//   analytic - genfit::ShellMaterialInterface: a thread-safe, non-TGeo material
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
#include <ShellMaterialInterface.h>
#include <MeasurementCreator.h>
#include <MathHelpers.h>
#include <HelixTrackModel.h>
#include <RKTrackRep.h>
#include <KalmanFitterRefTrack.h>
#include <Track.h>
#include <TrackPoint.h>
#include <Exception.h>

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
    genfit::MaterialEffects::getInstance()->init(new genfit::ShellMaterialInterface(airMat, siliconMat, {1.0, 1.05, 2.0, 2.05, 3.0, 3.05, 4.0, 4.05}, 20.));
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
