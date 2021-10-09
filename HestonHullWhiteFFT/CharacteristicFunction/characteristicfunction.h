// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef HESTONHULLWHITEFFT_CHARACTERISTICFUNCTION_CHARACTERISTICFUNCTION_H
#define HESTONHULLWHITEFFT_CHARACTERISTICFUNCTION_CHARACTERISTICFUNCTION_H

#include <cmath>
#include <complex>

#include "HestonHullWhiteFFT/HestonHullWhiteFFTInputs/HestonHullWhiteInputs.h"

namespace heston_hull_white {
    namespace fourier {

        namespace internal {

            // Calculates the characteristic function of the price under the heston model at the point evaluation_point. Uses a scaling factor hhw_inputs.itsScale_.

            // theorem 3, page 8
            std::complex<double> HestonCharacteristicFunction(const HestonHullWhiteInputs& hhw_inputs, std::complex<double> evaluation_point);

            // Calculates the characteristic function of the integral of the interest rate under the Hull-White model at the point evaluation_point. 

            // the "of integral" part in function definition is important here to explain the evaluation_point+std::complex<double>{0,1} in the hestonhullwhite characteristic function
            std::complex<double> HullWhiteCharacteristicFunctionOfIntegral(const HestonHullWhiteInputs& hhw_inputs, std::complex<double> evaluation_point);

        }

        // Calculates the characteristic function of the price under the heston hull-white model at the point evaluation_point. Uses a scaling factor hhw_inputs.itsScale_.

        std::complex<double> HestonHullWhiteCharacteristicFunction(const HestonHullWhiteInputs& hhw_inputs, std::complex<double> evaluation_point);


    }
}

#endif