// Implements Theorem 5.4, page 39 of https://media-exp1.licdn.com/dms/document/C4D2DAQH3DQXdHV78BA/profile-treasury-document-pdf-analyzed/0/1592809800927?e=1624611600&v=beta&t=1beaP4zemdBLAEvqJsVoOmtQr8KQLX4dmOKBE4-jmMg
// "Approximating the Heston-Hull-White Model"
// by Riaz Patel

// TO DO: I could slightly speed up each construction of a HestonHullWhiteInstance by creating temporary values that allow me to more efficiently construct the 
// "advance_log_asset_price_" and "advance_short_rate_" constants. That way I am not having to recompute divisions. This should be done using a delegating constructor.
// This would provide a tiny, tiny speed boost.

// TO DO: I currently approximate the integrals int_t^{t+dt} volatility(u) du by gamma*v(t)+(1-gamma)*v(t) and int_t^{t+dt} short_rate(u) du by gamma'*short_rate(t)+(1-gamma')*short_rate(t).
// There are ways to sample these integrals better. I can sample int_t^{t+dt} volatility(u) du GIVEN v(t) and v(t+dt). Exact method found here: 
// Broadie, M. and O. Kaya (2006), “Exact simulation of stochastic volatility and other affine jump diffusion processes, ” Operations Research, vol. 54, no. 2.
// That method is slow, but apparently there are "moment matching" techniques found here:
// Dufresne, D. (2001), “The integrated square-root process,” Working paper, University of Montreal.
// Andersen discusses these issues in small detail in Efficient Simulation of the Heston Stochastic Volatility Model.

// TO DO: I am currently approximating int_t ^{t+dt} \eta d W_v by inverting the Euler-Maruyama discretization for the V(t) process. This may or may not be the best option.
// I spent a long time thinking about it though, so there may not be other options. The essential problem is to sample int_t ^{t+dt} 1 d W_v given v(t) and v(t+dt), or alternatively,
// simulate the CIR process for v with a brownian bridge given int_t ^{t+dt} 1 d W_v. I do not believe this is possible in the QE framework. It may not even be a reasonable problem in CIR,
// and may require a differend model altogether.

// TO DO: Check what member variables aren't used after the constructor. Some are just used to compute constants and are not reused.

// TO DO: Replace QE constants with 4 separate constants to reduce multiplications.

#ifndef HESTONHULLWHITEMONTECARLO_STATETYPES_EXAMPLES_HESTONHULLWHITE_HESTONHULLWHITE_H
#define HESTONHULLWHITEMONTECARLO_STATETYPES_EXAMPLES_HESTONHULLWHITE_HESTONHULLWHITE_H

#include <cmath>
#include <concepts>
#include <memory>
#include <random>

#include "HestonHullWhiteMonteCarlo/StateTypes/Examples/HestonHullWhite/src/hestonhullwhiteInputparameters.h"
#include "HestonHullWhiteMonteCarlo/util/inversecumulativenormal.h"
#include "HestonHullWhiteMonteCarlo/util/InterestCurve/interestcurve.h"

namespace heston_hull_white {
	namespace monte_carlo {

		class HestonHullWhite {
		private:
			// starting random variables
			double log_asset_price_;
			double volatility_;
			double short_rate_;
			double integral_of_short_rate_; // e^{int_0 ^ T r(u) du}
			double custom_advance_total_time_; // don't overuse. Default to using standard advance.
			size_t standard_advance_steps_taken_; // used to calculate the current time s = steps_taken*epsilon for the short rate calculation. Used to avoid numerical error from repeated additions.


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

			//
			const double epsilon_; // derived from number of trials;

			// Saved for CustomAdvance
			const double gamma_1_;
			const double gamma_2_;
			const double gamma_3_;

			// cached constants. If epsilon can be variable, consider caching constants for each possible epsilon.
			//const double qe_constant_1_; // this is a constant that appears in the QE step. It is a function of epsilon. This is slightly different from what is found in H. Kammeyer's code.
			//const double qe_constant_2_; // this is a constant that appears in the QE step. It is a function of epsilon.
			//const double qe_constant_3_; // this is a constant that appears in the QE step. It is a function of epsilon.

			// cached constants. If epsilon can be variable, consider caching constants for each possible epsilon.
			const double qe_constant_0_;// this is a constant that appears in the QE step. It is a function of epsilon. This is slightly different from what is found in H. Kammeyer's code.
			const double qe_constant_1_; // this is a constant that appears in the QE step. It is a function of epsilon. 
			const double qe_constant_2_; // this is a constant that appears in the QE step. It is a function of epsilon.
			const double qe_constant_3_; // this is a constant that appears in the QE step. It is a function of epsilon.

			// constant needed for InvertEulerMaruyama_Volatility member function
			const double invert_euler_maruyama_constant_;

			// constants needed to advancing log asset price
			const double advance_log_asset_price_constant_0_;
			const double advance_log_asset_price_constant_1_;
			const double advance_log_asset_price_constant_2_;
			const double advance_log_asset_price_constant_3_;
			const double advance_log_asset_price_constant_4_;
			const double advance_log_asset_price_constant_5_;
			const double advance_log_asset_price_constant_6_;

			// constants needed to advancing short rate
			const double advance_short_rate_constant_0_;
			const double advance_short_rate_constant_1_;
			const double advance_short_rate_constant_2_;
			const double advance_short_rate_constant_3_;

			const double psi_constant_;

		public:
			//HestonHullWhite(double epsilon, double starting_price, double starting_volatility, double starting_short_rate, \
			//	double kappa, double mean_volatility, double sigma, double lambda, double eta, \
			//	double rho_s_v, double rho_s_r, double rho_v_r, std::shared_ptr<util::InterestCurve> interest_curve/*pass by value?*/, double gamma_1 = .5, double gamma_2 = .5, double gamma_3 = .5) : \
			//	log_asset_price_(std::log(starting_price)), volatility_(starting_volatility), short_rate_(starting_short_rate), integral_of_short_rate_(0), custom_advance_total_time_(0), standard_advance_steps_taken_(0), \
			//	kappa_(kappa), mean_volatility_(mean_volatility), sigma_(sigma), lambda_(lambda), eta_(eta), \
			//	rho_s_v_(rho_s_v), rho_s_r_(rho_s_r), rho_v_r_(rho_v_r), \
			//	interest_curve_(std::move(interest_curve)), \
			//	epsilon_(epsilon), \
			//	gamma_1_(gamma_1), gamma_2_(gamma_2), gamma_3_(gamma_3), \
			//	qe_constant_1_(std::exp(-kappa * epsilon)), qe_constant_2_(mean_volatility* (1.0 - qe_constant_1_)), qe_constant_3_((sigma* sigma / kappa)* (1.0 - qe_constant_1_)), \
			//	invert_euler_maruyama_constant_(kappa*epsilon), \
			//	advance_log_asset_price_constant_0_(gamma_1*epsilon), advance_log_asset_price_constant_1_((1.0-gamma_1)*epsilon), advance_log_asset_price_constant_2_(-kappa*rho_s_v*mean_volatility*epsilon/sigma), \
			//	advance_log_asset_price_constant_3_((gamma_2* ((((kappa* rho_s_v) / sigma) - .5)* epsilon)) - (rho_s_v / sigma)), advance_log_asset_price_constant_4_(((1-gamma_2)* ((((kappa* rho_s_v) / sigma) - .5)* epsilon)) + (rho_s_v / sigma)), \
			//	advance_log_asset_price_constant_5_(gamma_3*epsilon*(1.0-rho_s_v* rho_s_v)), advance_log_asset_price_constant_6_((1.0-gamma_3)* epsilon* (1.0 - rho_s_v * rho_s_v)), \
			//	advance_short_rate_constant_0_(std::exp(-lambda*epsilon)), advance_short_rate_constant_1_(eta*rho_v_r*std::sqrt((1.0/(2.0*lambda))*(1.0-std::exp(-2.0*lambda*epsilon)))), /*should use delegating constructor because (1.0-std::exp(-2.0*lambda*epsilon) appears a lot*/ \
			//	advance_short_rate_constant_2_(eta*(rho_s_r-rho_v_r*rho_s_v)*std::sqrt(epsilon*(1.0 / (2.0 * lambda))* (1.0 - std::exp(-2.0 * lambda * epsilon)))/(std::sqrt(1-rho_s_v* rho_s_v))), \
			//	advance_short_rate_constant_3_(eta* (std::sqrt(1.0-rho_v_r*rho_v_r-(std::pow(rho_s_r - rho_v_r * rho_s_v, 2)/(1.0-rho_s_v*rho_s_v))))* std::sqrt(epsilon* (1.0 / (2.0 * lambda))* (1.0 - std::exp(-2.0 * lambda * epsilon)))), \
			//	psi_constant_((eta* eta) / (2.0 * lambda * lambda)) {}
			HestonHullWhite(double epsilon, HestonHullWhiteInputParameters& hhw_inp) : \
				log_asset_price_(std::log(hhw_inp.starting_price)), volatility_(hhw_inp.starting_volatility), short_rate_(hhw_inp.starting_short_rate), integral_of_short_rate_(0), custom_advance_total_time_(0), standard_advance_steps_taken_(0), \
				kappa_(hhw_inp.kappa), mean_volatility_(hhw_inp.mean_volatility), sigma_(hhw_inp.sigma), lambda_(hhw_inp.lambda), eta_(hhw_inp.eta), \
				rho_s_v_(hhw_inp.rho_s_v), rho_s_r_(hhw_inp.rho_s_r), rho_v_r_(hhw_inp.rho_v_r), \
				interest_curve_(hhw_inp.interest_curve), \
				epsilon_(epsilon), \
				gamma_1_(hhw_inp.gamma_1), gamma_2_(hhw_inp.gamma_2), gamma_3_(hhw_inp.gamma_3), \
				/*qe_constant_1_(std::exp(-hhw_inp.kappa * epsilon)), qe_constant_2_(hhw_inp.mean_volatility* (1.0 - qe_constant_1_)), qe_constant_3_((hhw_inp.sigma* hhw_inp.sigma / hhw_inp.kappa)* (1.0 - qe_constant_1_)), */\
				qe_constant_0_(std::exp(-hhw_inp.kappa * epsilon)), qe_constant_1_(hhw_inp.mean_volatility* (1.0 - qe_constant_0_)), \
				qe_constant_2_(std::exp(-hhw_inp.kappa * epsilon)*((hhw_inp.sigma* hhw_inp.sigma / hhw_inp.kappa)* (1.0 - qe_constant_0_))), qe_constant_3_(hhw_inp.mean_volatility* (1.0 - qe_constant_0_)*((hhw_inp.sigma* hhw_inp.sigma / hhw_inp.kappa)* (1.0 - qe_constant_0_))*(0.5)), \
				invert_euler_maruyama_constant_(hhw_inp.kappa* epsilon), \
				advance_log_asset_price_constant_0_(hhw_inp.gamma_1* epsilon), advance_log_asset_price_constant_1_((1.0 - hhw_inp.gamma_1)* epsilon), advance_log_asset_price_constant_2_(-hhw_inp.kappa * hhw_inp.rho_s_v * hhw_inp.mean_volatility * epsilon / hhw_inp.sigma), \
				advance_log_asset_price_constant_3_((hhw_inp.gamma_2* ((((hhw_inp.kappa* hhw_inp.rho_s_v) / hhw_inp.sigma) - .5)* epsilon)) - (hhw_inp.rho_s_v / hhw_inp.sigma)), advance_log_asset_price_constant_4_(((1 - hhw_inp.gamma_2)* ((((hhw_inp.kappa* hhw_inp.rho_s_v) / hhw_inp.sigma) - .5)* epsilon)) + (hhw_inp.rho_s_v / hhw_inp.sigma)), \
				advance_log_asset_price_constant_5_(hhw_inp.gamma_3* epsilon* (1.0 - hhw_inp.rho_s_v * hhw_inp.rho_s_v)), advance_log_asset_price_constant_6_((1.0 - hhw_inp.gamma_3)* epsilon* (1.0 - hhw_inp.rho_s_v * hhw_inp.rho_s_v)), \
				advance_short_rate_constant_0_(std::exp(-hhw_inp.lambda * epsilon)),
				/*advance_short_rate_constant_1_(hhw_inp.eta* hhw_inp.rho_v_r* std::sqrt((1.0 / (2.0 * hhw_inp.lambda))* (1.0 - std::exp(-2.0 * hhw_inp.lambda * epsilon)))), */ /*should use delegating constructor because (1.0-std::exp(-2.0*lambda*epsilon) appears a lot*/ \
				/*advance_short_rate_constant_2_(hhw_inp.eta* (hhw_inp.rho_s_r - hhw_inp.rho_v_r * hhw_inp.rho_s_v)* std::sqrt(epsilon* (1.0 / (2.0 * hhw_inp.lambda))* (1.0 - std::exp(-2.0 * hhw_inp.lambda * epsilon))) / (std::sqrt(1 - hhw_inp.rho_s_v * hhw_inp.rho_s_v))),*/ \
				/*advance_short_rate_constant_3_(0.0141421), */\
				/* advance_short_rate_constant_3_(hhw_inp.eta* (std::sqrt(1.0 - hhw_inp.rho_v_r * hhw_inp.rho_v_r - (std::pow(hhw_inp.rho_s_r - hhw_inp.rho_v_r * hhw_inp.rho_s_v, 2) / (1.0 - hhw_inp.rho_s_v * hhw_inp.rho_s_v))))* std::sqrt(epsilon* (1.0 / (2.0 * hhw_inp.lambda))* (1.0 - std::exp(-2.0 * hhw_inp.lambda * epsilon)))), */ \
				/*advance_short_rate_constant_3_(std::sqrt(hhw_inp.eta* hhw_inp.eta* (std::sqrt(1.0 - hhw_inp.rho_v_r * hhw_inp.rho_v_r - (std::pow(hhw_inp.rho_s_r - hhw_inp.rho_v_r * hhw_inp.rho_s_v, 2) / (1.0 - hhw_inp.rho_s_v * hhw_inp.rho_s_v))))*  (1.0 / (2.0 * hhw_inp.lambda))* (1.0 - std::exp(-2.0 * hhw_inp.lambda * epsilon)))*/ /*std::sqrt(epsilon)),*/ \
				advance_short_rate_constant_1_(hhw_inp.eta* hhw_inp.rho_v_r*std::sqrt((1.0/(epsilon*/*extra epsilon here, see advance short rate in standard advance*/2.0* hhw_inp.lambda))*(1.0-std::exp(-2.0 * hhw_inp.lambda * epsilon)))), \
				advance_short_rate_constant_2_(hhw_inp.eta* (hhw_inp.rho_s_r - hhw_inp.rho_v_r * hhw_inp.rho_s_v) * std::sqrt(((1.0 / (2.0 * hhw_inp.lambda))* (1.0 - std::exp(-2.0 * hhw_inp.lambda * epsilon))/(1.0- hhw_inp.rho_s_v* hhw_inp.rho_s_v)))), \
				advance_short_rate_constant_3_(hhw_inp.eta*std::sqrt((1.0-hhw_inp.rho_v_r* hhw_inp.rho_v_r-(std::pow(hhw_inp.rho_s_r - hhw_inp.rho_v_r * hhw_inp.rho_s_v, 2)/(1.0 - hhw_inp.rho_s_v* hhw_inp.rho_s_v)))*((1.0 / (2.0 * hhw_inp.lambda))* (1.0 - std::exp(-2.0 * hhw_inp.lambda * epsilon))))), \
				psi_constant_((hhw_inp.eta* hhw_inp.eta) / (2.0 * hhw_inp.lambda * hhw_inp.lambda)) {}
				/*advance_log_asset_price_constant_0_(-(kappa * rho_s_v * mean_volatility * epsilon) / sigma), advance_log_asset_price_constant_1_(((((kappa* rho_s_v) / sigma) - .5)* epsilon) - (rho_s_v / sigma)), advance_log_asset_price_constant_2_(rho_s_v / sigma), advance_log_asset_price_constant_3_((1 - (rho_s_v * rho_s_v))* epsilon), \
				advance_short_rate_constant_1_(eta* std::sqrt(epsilon)*(rho_s_r/(std::sqrt(1.0-(rho_s_v*rho_s_v))))), advance_short_rate_constant_2_(eta*std::sqrt(epsilon*(1.0-((rho_s_r*rho_s_r)/(1.0-rho_s_v*rho_s_v))) )) {}*/

			double Value() const {
				return std::exp(log_asset_price_);
			}
			
			double DiscountFactor() const {
				return std::exp(-integral_of_short_rate_);
			}

			template <std::uniform_random_bit_generator _Engine>
			void StandardAdvance(_Engine& rng_engine) {
				// I believe I should do this first and not last, correct?
				integral_of_short_rate_ += short_rate_ * epsilon_;
				
				//std::cout << integral_of_short_rate_ << "\n";

				// QE step
				const double previous_volatility = volatility_;

				QuadraticExponential<_Engine>(rng_engine);
				// Advance stock price and short rate. This can be moved to a separate function, but it requires
				// the previous volatility.

				// First the asset price

				// get my random variables
				double implied_z_v_times_sqrt_epsilon_ = InvertEulerMaruyama_Volatility(previous_volatility);
				std::normal_distribution<double> a_normal_distribution{};
				double z_x = a_normal_distribution(rng_engine);
				double z_r = a_normal_distribution(rng_engine);

				// Save old values
				const double previous_log_asset_price = log_asset_price_;
				const double previous_short_rate = short_rate_;

				//ISSUE: The paper mentions a Z_v, but doesn't use it. I can understand it not being used since the randomness appears
				// in the QE step, but why is it mentioned?

				// compute some simplifying constants. WAIT. SOME OF THESE THINGS ARE BEING RECALCULATED EVERY LOOP. I SHOULD CACHE THEM.
				/*const double j = rho_s_v_ / sigma_;
				const double k = 1.0 - (rho_s_v_ * rho_s_v_);
				const double l = eta_ * std::sqrt(epsilon_);*/

				// old version
				//log_asset_price_ = previous_log_asset_price + previous_short_rate * epsilon_ + advance_log_asset_price_constant_0_ + advance_log_asset_price_constant_1_ * previous_volatility + advance_log_asset_price_constant_2_ * volatility_ + std::sqrt(advance_log_asset_price_constant_3_ * previous_volatility)*z_s;

				// old version
				//short_rate_ = previous_short_rate + lambda_ * (theta - previous_short_rate) * epsilon_ + advance_short_rate_constant_1_ * z_s + advance_short_rate_constant_2_ * z_r;
				
				// fixed version below
				// Note that advanced_short_rate_constant_1 is divided by std::sqrt(epsilon) to remove the additional sqrt(epsilon) factor present in implied_z_v_times_sqrt_epsilon_
				const double current_time = standard_advance_steps_taken_ * epsilon_ + custom_advance_total_time_;
				short_rate_ = advance_short_rate_constant_0_ * previous_short_rate + Psi(current_time+epsilon_) - Psi(current_time) * advance_short_rate_constant_0_ + \
					advance_short_rate_constant_1_ * implied_z_v_times_sqrt_epsilon_ + advance_short_rate_constant_2_ * z_x + advance_short_rate_constant_3_ * z_r;

				// fixed version below
				log_asset_price_ = previous_log_asset_price + advance_log_asset_price_constant_0_ * previous_short_rate + advance_log_asset_price_constant_1_ * short_rate_ + \
					advance_log_asset_price_constant_2_ + advance_log_asset_price_constant_3_ * previous_volatility + advance_log_asset_price_constant_4_ * volatility_ + \
					std::sqrt(advance_log_asset_price_constant_5_ * previous_volatility + advance_log_asset_price_constant_6_ * volatility_) * z_x;

				//log_asset_price_ = previous_log_asset_price + previous_short_rate * epsilon_ - (kappa_ * j * mean_volatility_ * epsilon_) + ((kappa_ * j) - (.5)) * epsilon_ - j)* (previous_volatility)+j * volatility_ + std::sqrt(k * previous_volatility); // from theorem 5.4
				//short_rate_ = previous_short_rate + lambda_ * (theta - previous_short_rate) * epsilon_ + l *( rho_s_r_ / (std::sqrt(k)) * z_s + std::sqrt(1.0 - ((rho_s_r_ * rho_s_r_) / k))*z_r ); // from theorem 5.4

				//WAIT. SOME OF THESE THINGS ARE BEING RECALCULATED EVERY LOOP. I SHOULD CACHE THEM.
				standard_advance_steps_taken_ += 1;
				return;
			}

			// prefer standard advance over this.
			template <std::uniform_random_bit_generator _Engine>
			void CustomAdvance(double epsilon, _Engine& rng_engine) {
				// I use a workaround to avoid having to essentially copy the rest of the code into here.
				// copy this wih different epsilon
				HestonHullWhite copy_of_this{std::exp(log_asset_price_), volatility_, short_rate_, kappa_, mean_volatility_, sigma_, lambda_, eta_, rho_s_v_, rho_s_r_, rho_v_r_, epsilon, interest_curve_, gamma_1_, gamma_2_, gamma_3_};
				copy_of_this.standard_advance_steps_taken_ = this->standard_advance_steps_taken_;
				copy_of_this.custom_advance_total_time_ = this->custom_advance_total_time_;
				copy_of_this.integral_of_short_rate_ = this->integral_of_short_rate_;
				// advance
				copy_of_this.StandardAdvance(rng_engine);
				// copy back
				this->log_asset_price_ = copy_of_this.log_asset_price_;
				this->volatility_ = copy_of_this.volatility_;
				this->short_rate_ = copy_of_this.short_rate_;
				this->integral_of_short_rate_ = copy_of_this.integral_of_short_rate_;
				// step forward
				custom_advance_total_time_ += epsilon;
			}

		private:

			// Computes implied W_v(t+dt)-W_v(dt) by inverting the Euler-Maruyama discretization
			double InvertEulerMaruyama_Volatility(double previous_volatility) {
				if (volatility_ == 0) {
					return 0;
				}
				else {
					return ((volatility_ - previous_volatility) - invert_euler_maruyama_constant_ * (mean_volatility_ - previous_volatility)) / (sigma_ * std::sqrt(previous_volatility));
				}
			}

			// I don't think this is possible because Milstein is quadratic
	/*		double InvertMilsteinScheme_Volatility(double previous_volatility) {
				
			}*/

			double Psi(double t) const //See documentation, equation 2.6, psi
			{
				return interest_curve_->InstForward(t) + psi_constant_ * std::pow(1.0 - std::exp(-lambda_ * t), 2);
				//return interest_curve_->InstForward(t) + ((itsEta * itsEta) / (2.0 * itsLambda * itsLambda)) * (1.0 - exp(-itsLambda * t)) * (1.0 - exp(-itsLambda * t));
			}



			// essentially advances the volatility
			template <std::uniform_random_bit_generator _Engine>
			void QuadraticExponential(_Engine& _Eng) {

				// compute m and S2 from theorem 2.1. (Note: don't get confused here just because I compute Expectation and Variance. volatility isn't normal.)
				//const double conditional_expectation_of_future_volatility = volatility_ * qe_constant_1_ + qe_constant_2_; // maybe I should cache the RHS of this equation
				//const double conditional_variance_of_future_volatility = volatility_ * qe_constant_1_ * qe_constant_3_ + .5 * qe_constant_2_ * qe_constant_3_; // maybe I should cache the RHS of this equation
				const double conditional_expectation_of_future_volatility = volatility_ * qe_constant_0_ + qe_constant_1_; // maybe I should cache the RHS of this equation
				const double conditional_variance_of_future_volatility = volatility_ * qe_constant_2_ + qe_constant_3_; // maybe I should cache the RHS of this equation
				// compute Phi
				const double Psi = conditional_variance_of_future_volatility / conditional_expectation_of_future_volatility * conditional_expectation_of_future_volatility;
				const double Psi_c = 1.5;
			
				double random_value_U_v = std::generate_canonical<double, static_cast<size_t>(std::numeric_limits<double>::digits)>(_Eng);

				if (Psi > Psi_c) {
					const double p_const = (Psi-1.0)/(Psi+1.0); // 

					if (random_value_U_v < p_const) {
						volatility_ = 0;
						return;
					}
					else {
						volatility_ = (conditional_expectation_of_future_volatility / (1.0- p_const)) * std::log((1 - p_const) / (1 - random_value_U_v)); // is this the most numerically stable version?
						return;
					}
				}
				else {
					const double b1 = (2.0 / (Psi));
					const double b = b1 - 1.0 + std::sqrt(b1 * (b1 - 1)); // since Psi < 1.5, we have b1 > 1, so b >0 and computing the max is unnecessary.
					//double b = std::max(b1-1.0+std::sqrt(b1*(b1-1)), 0); //This is the computation from Kienitz code. max is unnecessary.
					const double shifted_normal = b + internal::InverseCumulativeNormal(0, 1)(random_value_U_v); // this seems unnecessary. Why not just sample a normal value?
					//double a = conditional_expectation_of_future_volatility/(1.0+b*b);
					//volatility_ = a * (shifted_normal * shifted_normal);
					volatility_ = (conditional_expectation_of_future_volatility / (1.0 + (b * b))) * (shifted_normal * shifted_normal);
					return;
				}



			}
		};
	}
}

#endif