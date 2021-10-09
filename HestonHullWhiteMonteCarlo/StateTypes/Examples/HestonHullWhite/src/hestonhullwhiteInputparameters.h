#ifndef HESTONHULLWHITEMONTECARLO_STATETYPES_EXAMPLES_HESTONHULLWHITE_SRC_HESTONHULLWHITEINPUTPARAMETERS_H
#define HESTONHULLWHITEMONTECARLO_STATETYPES_EXAMPLES_HESTONHULLWHITE_SRC_HESTONHULLWHITEINPUTPARAMETERS_H

#include <memory>

#include "HestonHullWhiteMonteCarlo/util/InterestCurve/interestcurve.h"

namespace heston_hull_white {
	namespace monte_carlo {

		struct HestonHullWhiteInputParameters {
		public:
			HestonHullWhiteInputParameters(double starting_price_in, double starting_volatility_in, double starting_short_rate_in, \
				double kappa_in, double mean_volatility_in, double sigma_in, double lambda_in, double eta_in, \
				double rho_s_v_in, double rho_s_r_in, double rho_v_r_in, std::shared_ptr<util::InterestCurve> interest_curve_in/*pass by value?*/, double gamma_1_in = .5, double gamma_2_in = .5, double gamma_3_in = .5) :
			starting_price(starting_price_in), starting_volatility(starting_volatility_in), starting_short_rate(starting_short_rate_in), \
				kappa(kappa_in), mean_volatility(mean_volatility_in), sigma(sigma_in), lambda(lambda_in), eta(eta_in), \
				rho_s_v(rho_s_v_in), rho_s_r(rho_s_r_in), rho_v_r(rho_v_r_in), \
				interest_curve(std::move(interest_curve_in)), gamma_1(gamma_1_in), gamma_2(gamma_2_in), gamma_3(gamma_3_in) {}
			double starting_price;
			double starting_volatility;
			double starting_short_rate;
			double kappa;
			double mean_volatility;
			double sigma;
			double lambda;
			double eta;
			double rho_s_v;
			double rho_s_r;
			double rho_v_r;
			std::shared_ptr<util::InterestCurve> interest_curve;/*pass by value?*/
			double gamma_1 = .5;
			double gamma_2 = .5;
			double gamma_3 = .5;

		};
	}
}

#endif