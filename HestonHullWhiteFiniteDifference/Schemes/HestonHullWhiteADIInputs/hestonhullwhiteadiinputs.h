#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_HESTONHULLWHITEADIINPUTS_HESTONHULLWHITEADIINPUTS_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_HESTONHULLWHITEADIINPUTS_HESTONHULLWHITEADIINPUTS_H

#include <memory>

#include "HestonHullWhiteFiniteDifference/util/InterestCurve/interestcurve.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			struct HestonHullWhiteADIInputs {

				// constant input variables
				const double kappa_; // mean reversion speed of volatility
				const double mean_volatility_; // mean vol?
				const double sigma_; // vol of vol?
				const double lambda_; // mean reversion speed of interest rate
				const double eta_; // vol of interest rate?

				const double rho_s_v_; // correlation between W_s(t) and W_v(t)
				const double rho_s_r_; // correlation between W_s(t) and W_r(t)
				const double rho_v_r_; // correlation between W_v(t) and W_r(t)

				const std::shared_ptr<util::InterestCurve> interest_curve_;

			};

		}
	}
}

#endif