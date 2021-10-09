// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef HESTONHULLWHITEFFT_EUROPEANCALLPRICER_EUROPEANCALLPRICER_H
#define HESTONHULLWHITEFFT_EUROPEANCALLPRICER_EUROPEANCALLPRICER_H

// TODO: Use itsScale
#include <cassert>
#include <complex>
#include <numbers>
#include <vector>

#include "HestonHullWhiteFFT/HestonHullWhiteFFTInputs/HestonHullWhiteInputs.h"
#include "HestonHullWhiteFFT/HestonHullWhiteFFTInputs/FFTInputs.h"
#include "HestonHullWhiteFFT/AlphaOptimization/alphaoptimization.h"
#include "HestonHullWhiteFFT/util/fft/fft.h"

namespace heston_hull_white {
	namespace fourier {
		
		class EuropeanCallPricer {
		public:
			EuropeanCallPricer(const HestonHullWhiteInputs& hhw_inputs, const FFTInputs& fft_inputs/*, double itsScale = 1.0 */);

			double operator()(double strike_);


		private:
			double fft_b_;
			double fft_lambda_;
			double itsScale_;
			double alpha_;
			// need some kind of vector to store prices
			// need to store constants that speed up price accesses in operator()
			std::vector<double> prices_;
			//std::vector<double> prices_;

		};

	}
}


#endif