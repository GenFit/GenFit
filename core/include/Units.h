#pragma once

#include "EigenMatrixTypedefs.h"

namespace genfit {

    namespace Unit {

        constexpr Scalar GeV = 1;
        constexpr Scalar MeV = 1e-3 * Unit::GeV;
        constexpr Scalar keV = 1e-6 * Unit::GeV;
        constexpr Scalar eV = 1e-9 * Unit::GeV;

    }

}