#ifndef genfit_Material_h
#define genfit_Material_h

#include <TObject.h>

namespace genfit {
    typedef double Scalar;

    struct Material {
        Scalar density;  /// Density in g / cm^3
        Scalar Z;  /// Atomic number
        Scalar A;  /// Mass number in g / mol
        Scalar radiationLength;  /// Radiation Length in cm
        Scalar mEE;  /// Mean excitaiton energy in eV

        Material() : density(0), Z(0), A(0), radiationLength(0), mEE(0) {}

        Material(Scalar density_, Scalar Z_, Scalar A_, Scalar radiationLength_, Scalar mEE_) :
                density(density_), Z(Z_), A(A_), radiationLength(radiationLength_), mEE(mEE_) {}

        virtual ~Material() {};

        void Print(const Option_t* = "") const;

        ClassDef(Material, 1)
    };

    bool operator==(const Material &lhs, const Material &rhs);

    bool operator!=(const Material &lhs, const Material &rhs);

}

#endif