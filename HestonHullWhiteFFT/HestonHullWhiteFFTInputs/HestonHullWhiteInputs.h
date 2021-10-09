// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef HESTONHULLWHITEFFT_HESTONHULLWHITEFFTINPUTS_HESTONHULLWHITEINPUTS_H
#define HESTONHULLWHITEFFT_HESTONHULLWHITEFFTINPUTS_HESTONHULLWHITEINPUTS_H

#include <memory>

#include "HestonHullWhiteFFT/util/interestcurve.h"

namespace heston_hull_white {
	namespace fourier {

		struct HestonHullWhiteInputs {
			// call option inputs
			const double maturity_;
			const double spot_price_;
			const double starting_volatility_;
			const double strike_; // pricer gives prices for multiple strikes

			double itsScale_;  // usually strike

			// constant input variables
			const double kappa_; // mean reversion speed of volatility
			const double mean_volatility_; // mean vol?
			const double sigma_; // vol of vol?
			const double lambda_; // mean reversion speed of interest rate
			const double eta_; // vol of interest rate?

			const double rho_s_v_; // correlation between W_s(t) and W_v(t)

			const std::shared_ptr<heston_hull_white::fourier::util::InterestCurve> interest_curve_;

		};
	}
}

#endif