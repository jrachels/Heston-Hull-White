#include <algorithm>

#include "HestonHullWhiteTests/catch.hpp"
#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"

TEST_CASE("Basic Accuracy test for ModifiedCraigSneyd") {

	// This trial uses s_0 = 105, K = 110, t = 5 from the Kammeyer and Kienitz paper

	double K = 110;
	double T = 5;

	// specifying the interest curve

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

	std::shared_ptr< heston_hull_white::finite_difference::util::InterestCurve > interest_curve_ptr = std::make_shared< heston_hull_white::finite_difference::util::InterestCurve>(discount_curve);


	heston_hull_white::finite_difference::examples::HestonHullWhiteADIInputs heston_hull_white_adi_input_parameters{ .4 /*kappa*/, .02 /*mean_volatility*/,  .2 /*sigma*/, .7 /*lambda*/, .2 /*eta*/, -.3 /*-.3*/ /*rho_s_v*/, 0 /* rho_s_r*/, 0 /*rho_v_r = .05?*/, interest_curve_ptr /*std::shared_ptr<util::InterestCurve> interest_curve*/ };

	// values chosen based on recommendations from page 4 and 5 of ADI paper
	heston_hull_white::finite_difference::examples::ADIDiscretizationInputs discretization_inputs{ \
		40 /* m1_ */, K / 20.0 /* d_1 */, std::max(.5, std::exp(-T / (4.0))) * K /* S_left */, K /* S_right */, 14 * K /* S_max*/, \
		40 /* m_2 */, .02 /*d_2 = V_max/500*/, 10 /* V_max */, \
		40 /* m_3 */, 0.0025 /* d_3 = R_max/400 */, .1 /* c */, 1 /* R_max */ };

	// my intrinsic pricer:
	const heston_hull_white::finite_difference::examples::VanillaCall vanilla_call_option_pricer_110{ K };

	heston_hull_white::finite_difference::examples::ADIDiscretization discretization{ discretization_inputs };


	heston_hull_white::finite_difference::FiniteDifference< heston_hull_white::finite_difference::examples::ModifiedCraigSneyd<heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>>, \
		heston_hull_white::finite_difference::examples::ADIDiscretization, heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> \
		fd{ 5, 0, 1000, heston_hull_white_adi_input_parameters, discretization_inputs, vanilla_call_option_pricer_110 };

	auto [s_index, v_index, r_index] = fd.GetIndex<double, double, double>(110.0, 0.018, interest_curve_ptr->ShortRate());
	double actual_price = 30.5719092871391;

	REQUIRE(fd.PriceByIndex(999, s_index - 2, v_index - 2, r_index - 2) < actual_price);
	REQUIRE(fd.PriceByIndex(999, s_index, v_index, r_index) > actual_price);


}