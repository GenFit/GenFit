#pragma once

#include "EigenMatrixTypedefs.h"

namespace genfit {

    namespace Unit {

        constexpr Scalar GeV = 1;
        constexpr Scalar MeV = 10e-3 * Unit::GeV;
        constexpr Scalar keV = 10e-6 * Unit::GeV;
        constexpr Scalar eV = 10e-9 * Unit::GeV;

    }

}