#pragma once

#include <Math/Vector3D.h>
#include <TVector3.h>


namespace genfit {
  namespace VectorUtils {

    /// @brief Get the individual coordinates of a 3-vector by index
    /// @param[in] a input vector
    /// @param[in] i index
    /// @return a[i]
    static inline double getCoord(const ROOT::Math::XYZVector& a, const size_t i) {
      switch (i) {
        case 0:
          return a.X();
        case 1:
          return a.Y();
        case 2:
          return a.Z();
        default:
          return 0.; // To make compiler happy, can never happen
      }
    }

    /**
     * @brief Helper function to calculate an orthogonal vector, same logic as in TVector3
     * @param[in] a vector for which to calculate the orthogonal
     * @return vector orthognoal to a
     */
    inline ROOT::Math::XYZVector Orthogonal(const ROOT::Math::XYZVector& a) {
      Double_t xx = a.X() < 0.0 ? -a.X() : a.X();
      Double_t yy = a.Y() < 0.0 ? -a.Y() : a.Y();
      Double_t zz = a.Z() < 0.0 ? -a.Z() : a.Z();
      if (xx < yy) {
          return xx < zz ? ROOT::Math::XYZVector(0,a.Z(),-a.Y()) : ROOT::Math::XYZVector(a.Y(),-a.X(),0);
      } else {
          return yy < zz ? ROOT::Math::XYZVector(-a.Z(),0,a.X()) : ROOT::Math::XYZVector(a.Y(),-a.X(),0);
      }
    }

    /**
     * @brief Helper function to convert XYZVector to TVector3
     * @param[in] a XYZVector to convert to TVector3
     */
    static constexpr auto XYZToTVector = [](const ROOT::Math::XYZVector& a)
    {
      return TVector3(a.X(), a.Y(), a.Z());
    };

    /**
     * @brief Helper function to print the vector coordinates
     * @param[in] a Vector to print
     * @param[in] printOut stream to write to
     */
    static constexpr auto PrintVec = [](const ROOT::Math::XYZVector& a, std::ostream& printOut) {
      printOut << "(x,y,z)=(" << a.X() << "," << a.Y() << "," << a.Z() << "), (r,theta,phi)=("
               << a.R() << "," << a.Theta() * 180. / M_PI << "," << a.Phi() * 180. / M_PI << ")";
    };
    
    /**
     * @brief Set vector by polar coordinates.
     * @param[out] vector Vector.
     * @param[in]  mag    Magnitude.
     * @param[in]  theta  Polar angle.
     * @param[in]  phi    Azimuthal angle.
     */
    inline void SetMagThetaPhi(ROOT::Math::XYZVector& vector,
                               double mag, double theta, double phi)
    {
      const double amag = std::abs(mag);
      const double sinTheta = std::sin(theta);
      const double x = amag * sinTheta * std::cos(phi);
      const double y = amag * sinTheta * std::sin(phi);
      const double z = amag * std::cos(theta);
      vector.SetXYZ(x, y, z);
    }

    /**
     * @brief Set vector magnitude mag
     * @param[inout] vector Vector
     * @param[in]    mag    Magnitude
     */
    inline void SetMag(ROOT::Math::XYZVector& vector, double mag)
    {
      SetMagThetaPhi(vector, mag, vector.Theta(), vector.Phi());
    }

    /**
     * Set vector azimuthal angle theta
     * @param[inout] vector Vector
     * @param[in]    theta  Azimuthal angle
     */
    inline void SetTheta(ROOT::Math::XYZVector& vector, double theta)
    {
      SetMagThetaPhi(vector, vector.R(), theta, vector.Phi());
    }

    /**
     * Set vector polar angle phi
     * @param[inout] vector Vector
     * @param[in]    phi    Polar angle
     */
    inline void SetPhi(ROOT::Math::XYZVector& vector, double phi)
    {
      SetMagThetaPhi(vector, vector.R(), vector.Theta(), phi);
    }
  }
}
