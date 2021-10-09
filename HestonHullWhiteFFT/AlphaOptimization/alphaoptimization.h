// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef HESTONHULLWHITEFFT_ALPHAOPTIMIZATION_ALPHAOPTIMIZATION_H
#define HESTONHULLWHITEFFT_ALPHAOPTIMIZATION_ALPHAOPTIMIZATION_H

#include <cmath>
#include <complex>

#include "HestonHullWhiteFFT/CharacteristicFunction/characteristicfunction.h"
#include "HestonHullWhiteFFT/HestonHullWhiteFFTInputs/HestonHullWhiteInputs.h"
#include "HestonHullWhiteFFT/util/lbfgsb/ap/ap.h"
#include "HestonHullWhiteFFT/util/lbfgsb/lbfgsb.h"

namespace heston_hull_white {
    namespace fourier {


		namespace internal {

			// Function to be maximized by the AlphaOptimizatizer
			double AlphaObjectiveFunction(const HestonHullWhiteInputs& hhw_inputs, const ap::real_1d_array& x);

		}

		// Computes optimal alpha offset for Fourier Transform.
		double AlphaOptimizatizer(const HestonHullWhiteInputs& hhw_inputs);

    }
}

#endif