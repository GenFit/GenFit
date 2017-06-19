#include "Material.h"

#include "IO.h"

namespace genfit {

    bool operator== (const Material& lhs, const Material& rhs){
        if (&lhs == &rhs)
            return true;

        return !(lhs.density != rhs.density or
                 lhs.Z != rhs.Z or
                 lhs.A != rhs.A or
                 lhs.radiationLength != rhs.radiationLength or
                 lhs.mEE != rhs.mEE);

    }

    bool operator!= (const Material& lhs, const Material& rhs) {
        return !(lhs==rhs);
    }

    void Material::Print(const Option_t*) const {
        printOut << "Density = " << density << ", \t"
                 << "Z = " << Z << ", \t"
                 << "A = " << A << ", \t"
                 << "radiationLength = " << radiationLength << ", \t"
                 << "mEE = " << mEE << "\n";
    }
    
}