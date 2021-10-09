#include "HestonHullWhiteTests/catch.hpp"
#include "HestonHullWhiteMonteCarlo/European/Include/europeanmontecarlo.h"

#include <array>
#include <memory>
#include <random>
#include <vector>

// TODO: just simplify all this with "using"
TEST_CASE("Compilation/Accuracy check for MonteCarlo Pricer") {
	// step 1: create a random number generator
	// std::mt19937 mt{ std::random_device()() }; // this is what you should usually do
	std::mt19937 mt{ 0 }; // for deterministic outcomes
	// step 2: created a CombSim
		// I choose MultipleDerivativesSimulator
 
		// step i: specify the interest curve
	std::vector < std::array<double, 2>> discount_curve;
	discount_curve.push_back({ 0, 1 });
	discount_curve.push_back({ 0.003, 0.999884333380315 });
	discount_curve.push_back({ 0.083, 0.996803132736937 });
	discount_curve.push_back({ 0.167, 0.993568709230647 });
	discount_curve.push_back({ 0.25, 0.990285301195274 });
	discount_curve.push_back({ 0.333, 0.986945903402709 });
	discount_curve.push_back({ 0.417, 0.983557350486521 });
	discount_curve.push_back({ 0.5, 0.980185549124449 });
	discount_curve.push_back({ 0.583, 0.976782934344041 });
	discount_curve.push_back({ 0.667, 0.973361992614499 });
	discount_curve.push_back({ 0.75, 0.96997679330522 });
	discount_curve.push_back({ 0.833, 0.966616749933289 });
	discount_curve.push_back({ 0.917, 0.96291431795816 });
	discount_curve.push_back({ 1, 0.959904777446077 });
	discount_curve.push_back({ 2, 0.920091903961326 });
	discount_curve.push_back({ 3, 0.882870065420196 });
	discount_curve.push_back({ 4, 0.847186544281939 });
	discount_curve.push_back({ 5, 0.812742515687365 });
	discount_curve.push_back({ 6, 0.779459552415061 });
	discount_curve.push_back({ 7, 0.747152463119429 });
	discount_curve.push_back({ 8, 0.715745016074346 });
	discount_curve.push_back({ 9, 0.68513872380846 });
	discount_curve.push_back({ 10, 0.655753392359115 });
	discount_curve.push_back({ 11, 0.627333845297308 });
	discount_curve.push_back({ 12, 0.599226698198774 });
	discount_curve.push_back({ 13, 0.572763319281569 });
	discount_curve.push_back({ 14, 0.547259133751455 });
	discount_curve.push_back({ 15, 0.52344199625308 });
	discount_curve.push_back({ 16, 0.499646068368557 });
	discount_curve.push_back({ 17, 0.477507905873099 });
	discount_curve.push_back({ 18, 0.456481811728753 });
	discount_curve.push_back({ 19, 0.436385788738282 });
	discount_curve.push_back({ 20, 0.41735025383105 });
	discount_curve.push_back({ 21, 0.399187111819286 });
	discount_curve.push_back({ 22, 0.381865611666566 });
	discount_curve.push_back({ 23, 0.365435617455498 });
	discount_curve.push_back({ 24, 0.349786183601181 });
	discount_curve.push_back({ 25, 0.334806921914717 });
	discount_curve.push_back({ 26, 0.320548897004994 });
	discount_curve.push_back({ 27, 0.306983265264429 });
	discount_curve.push_back({ 28, 0.29408180091705 });
	discount_curve.push_back({ 29, 0.282443547729164 });
	discount_curve.push_back({ 30, 0.269929224010243 });

	// TODO: This line is pointless, right?
	heston_hull_white::monte_carlo::util::InterestCurve interest_curve{discount_curve};

	std::shared_ptr< heston_hull_white::monte_carlo::util::InterestCurve > interest_curve_ptr = std::make_shared< heston_hull_white::monte_carlo::util::InterestCurve>(discount_curve);

	// step ii: specify Heston Hull White input parameters
	
	heston_hull_white::monte_carlo::HestonHullWhiteInputParameters heston_hull_white_input_parameters(105 /*starting_price*/, .018 /*starting_volatility*/, /*interest_curve.YieldCurve(0.0)*/ - std::log(0.9998843333803) / (0.003)/*starting_short_rate*/, \
		.4 /*kappa*/, .02 /*mean_volatility*/, .2 /*sigma*/, .7 /*lambda*/, .2 /*eta*/, \
		- .3 /*rho_s_v*/, 0 /* rho_s_r*/, 0/*rho_v_r = .05?*/, interest_curve_ptr/*std::shared_ptr<util::InterestCurve> interest_curve*/);

	//heston_hull_white::monte_carlo::HestonHullWhiteInputParameters heston_hull_white_input_parameters(110 /*starting_price*/, .02 /*starting_volatility*/, interest_curve.YieldCurve(0.0) /*-std::log(0.9998843333803) / (0.003) starting_short_rate */ , \
	//	.5 /*kappa*/, .015 /*mean_volatility*/, .3 /*sigma*/, .5 /*lambda*/, .3 /*eta*/, \
	//	- .6 /*rho_s_v*/, 0 /* rho_s_r*/, 0/*rho_v_r = .05?*/, interest_curve_ptr/*std::shared_ptr<util::InterestCurve> interest_curve*/);

	//heston_hull_white::monte_carlo::HestonHullWhiteInputParameters heston_hull_white_input_parameters(105 /*starting_price*/, .018 /*starting_volatility*/, -std::log(0.999884) / (0.003)/*starting_short_rate*/, \
	//	.4 /*kappa*/, .02 /*mean_volatility*/, .2 /*sigma*/, .7 /*lambda*/, .2 /*eta*/, \
	//	- .3 /*rho_s_v*/, -.1 /* rho_s_r*/, .05/*rho_v_r*/, std::make_shared<heston_hull_white::util::InterestCurve>(interest_curve)/*std::shared_ptr<util::InterestCurve> interest_curve*/);

	//HestonHullWhite(double starting_price, double starting_volatility, double starting_short_rate, \
	//	double kappa, double mean_volatility, double sigma, double lambda, double eta, \
	//	double rho_s_v, double rho_s_r, double rho_v_r, double epsilon, std::shared_ptr<util::InterestCurve> interest_curve/*pass by value?*/, double gamma_1 = .5, double gamma_2 = .5, double gamma_3 = .5)

	//heston_hull_white::monte_carlo::HestonHullWhite state{ .01 /*epsilon*/, 105 /*starting_price*/, .018 /*starting_volatility*/, -std::log(0.999884) / (0.003)/*starting_short_rate*/, \
	//	.4 /*kappa*/, .02 /*mean_volatility*/, .2 /*sigma*/, .7 /*lambda*/, .2 /*eta*/, \
	//	-.3 /*rho_s_v*/, -.1 /* rho_s_r*/, .05/*rho_v_r*/, std::make_shared<heston_hull_white::util::InterestCurve>(interest_curve)/*std::shared_ptr<util::InterestCurve> interest_curve*/};


		// step iii: specify path independent
	heston_hull_white::monte_carlo::examples::VanillaCall1<heston_hull_white::monte_carlo::HestonHullWhite> vanilla_european_call_1{ 80 };
	heston_hull_white::monte_carlo::examples::VanillaCall2<heston_hull_white::monte_carlo::HestonHullWhite> vanilla_european_call_2(std::vector<double>{80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135});
	heston_hull_white::monte_carlo::examples::VanillaCall3<heston_hull_white::monte_carlo::HestonHullWhite, 12> vanilla_european_call_3(std::array<double, 12>{ 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135 });
	heston_hull_white::monte_carlo::PathIndependents<heston_hull_white::monte_carlo::HestonHullWhite,
		heston_hull_white::monte_carlo::examples::VanillaCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::VanillaCall2<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::VanillaCall3<heston_hull_white::monte_carlo::HestonHullWhite, 12>> \
		my_path_independents{ vanilla_european_call_1, vanilla_european_call_2, vanilla_european_call_3 };

	//PathIndependents()
		// step iv: specify path dependents
	heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite> asian_european_call_1{ 80 };
	heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite> asian_european_call_1_copy = asian_european_call_1;
	heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite> asian_european_call_2{ 100 };
	heston_hull_white::monte_carlo::PathDependents< \
		heston_hull_white::monte_carlo::HestonHullWhite, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>> \
		my_path_dependents{ asian_european_call_1, asian_european_call_1_copy, asian_european_call_2 };

		// step v: Call the MultipleDerivativesSimulator
	heston_hull_white::monte_carlo::european::MultipleDerivativesSimulator<heston_hull_white::monte_carlo::HestonHullWhite, \
		heston_hull_white::monte_carlo::PathIndependents<heston_hull_white::monte_carlo::HestonHullWhite,
		heston_hull_white::monte_carlo::examples::VanillaCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::VanillaCall2<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::VanillaCall3<heston_hull_white::monte_carlo::HestonHullWhite, 12>>, \
		heston_hull_white::monte_carlo::PathDependents< \
		heston_hull_white::monte_carlo::HestonHullWhite, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>>, \
		heston_hull_white::monte_carlo::HestonHullWhiteInputParameters>
		my_simulator(my_path_independents, my_path_dependents, 5 /*maturity_time*/, 1000 /*num_subintervals*/, heston_hull_white_input_parameters);
		//my_simulator(my_path_independents, my_path_dependents, 5 /*maturity_time*/, 1000 /*num_subintervals*/, heston_hull_white_input_parameters);

			
		
	// Step 3: call MonteCarloSimulation
	heston_hull_white::monte_carlo::european::MonteCarloSimulation simulation_result(mt, std::make_shared<heston_hull_white::monte_carlo::european::MultipleDerivativesSimulator<heston_hull_white::monte_carlo::HestonHullWhite, \
		heston_hull_white::monte_carlo::PathIndependents<heston_hull_white::monte_carlo::HestonHullWhite,
		heston_hull_white::monte_carlo::examples::VanillaCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::VanillaCall2<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::VanillaCall3<heston_hull_white::monte_carlo::HestonHullWhite, 12>>, \
		heston_hull_white::monte_carlo::PathDependents< \
		heston_hull_white::monte_carlo::HestonHullWhite, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>, \
		heston_hull_white::monte_carlo::examples::AsianCall1<heston_hull_white::monte_carlo::HestonHullWhite>>, \
		heston_hull_white::monte_carlo::HestonHullWhiteInputParameters>> \
		(my_simulator), 10000);
	std::vector<double> result = simulation_result.MeanOutcomes();

	// accuracy checks taken from Kienitz paper part ii

	REQUIRE(result[0] == Approx(45.3648403933267).margin(.2));
	REQUIRE(result[1] == result[0]);
	REQUIRE(result[2] == Approx(42.5401009008153).margin(.2));
	REQUIRE(result[3] == Approx(39.8623300348626).margin(.2));
	REQUIRE(result[4] == Approx(37.3299192484043).margin(.2));
	REQUIRE(result[5] == Approx(34.9400403063771).margin(.2));
	REQUIRE(result[6] == Approx(32.6888896067872).margin(.2));
	REQUIRE(result[7] == Approx(30.5719092871391).margin(.2));
	REQUIRE(result[8] == Approx(28.5839821655106).margin(.2));
	REQUIRE(result[9] == Approx(26.7195999278723).margin(.2));
	REQUIRE(result[10] == Approx(24.973005503804).margin(.2));
	REQUIRE(result[11] == Approx(23.3383114921657).margin(.2));
	REQUIRE(result[12] == Approx(21.8095969842613).margin(.2));
	REQUIRE(result[13] == result[1]);
	REQUIRE(result[14] == result[2]);
	REQUIRE(result[15] == result[3]);
	REQUIRE(result[16] == result[4]);
	REQUIRE(result[17] == result[5]);
	REQUIRE(result[18] == result[6]);
	REQUIRE(result[19] == result[7]);
	REQUIRE(result[20] == result[8]);
	REQUIRE(result[21] == result[9]);
	REQUIRE(result[22] == result[10]);
	REQUIRE(result[23] == result[11]);
	REQUIRE(result[24] == result[12]);

	// values computed using this algorithm

	REQUIRE(result[25] == Approx(28.9138).margin(.2));
	REQUIRE(result[26] == result[25]);
	REQUIRE(result[27] == Approx(16.9952).margin(.2));

}