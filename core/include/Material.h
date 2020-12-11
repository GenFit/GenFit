#ifndef genfit_Material_h
#define genfit_Material_h

#include <TObject.h>

namespace genfit {

    struct Material {
        double density;  /// Density in g / cm^3
        double Z;  /// Atomic number
        double A;  /// Mass number in g / mol
        double radiationLength;  /// Radiation Length in cm
        double mEE;  /// Mean excitaiton energy in eV

        Material() : density(0), Z(0), A(0), radiationLength(0), mEE(0) {}

        Material(double density_, double Z_, double A_, double radiationLength_, double mEE_) :
                density(density_), Z(Z_), A(A_), radiationLength(radiationLength_), mEE(mEE_) {}

        Material(const Material &material) = default;

        virtual ~Material() {};

        void Print(const Option_t* = "") const;

        ClassDef(Material, 1)
    };

    bool operator==(const Material &lhs, const Material &rhs);

    bool operator!=(const Material &lhs, const Material &rhs);

}

#endif
