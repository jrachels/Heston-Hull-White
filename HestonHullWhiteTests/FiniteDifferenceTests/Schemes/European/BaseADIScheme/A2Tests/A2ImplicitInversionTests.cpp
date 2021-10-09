#include "HestonHullWhiteTests/catch.hpp"

#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"
#include "HestonHullWhiteTests/FiniteDifferenceTests/Schemes/European/BaseADIScheme/A2Tests/A2tester.h"

TEST_CASE("A2 Inversion Test") {

	double K = 110;
	double T = 5;
	double theta = 2.0 / 3.0;
	double epsilon_t = .01;
	double kappa = .4;
	double mean_volatility = .02;
	double sigma = .2;
	

	heston_hull_white::finite_difference::examples::ADIDiscretizationInputs discretization_inputs{ \
	20 /* m1_ */, K / 20.0 /* d_1 */, std::max(.5, std::exp(-T / (4.0))) * K /* S_left */, K /* S_right */, 14 * K /* S_max*/, \
	20 /* m_2 */, .02 /*d_2 = V_max/500*/, 10 /* V_max */, \
	20 /* m_3 */, 0.0025 /* d_3 = R_max/400 */, .1 /* c */, 1 /* R_max */ };


	std::shared_ptr<heston_hull_white::finite_difference::examples::ADIDiscretization> discretization = \
		std::make_shared< heston_hull_white::finite_difference::examples::ADIDiscretization>(discretization_inputs);

	std::size_t mean_vol_index = discretization->GetIndexV(mean_volatility);

	heston_hull_white::finite_difference::examples::A2<heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> a2{ mean_volatility, kappa, sigma , mean_vol_index, discretization};

	heston_hull_white::finite_difference::examples::A2ImplicitInversion< heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> a2_implicit_inverter{ a2 , theta, epsilon_t };

	// here I use transform in place because that is what I used in the Douglas algorithm

	heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double> rhs = discretization->IntrinsicOptionValues<heston_hull_white::finite_difference::examples::VanillaCall>(heston_hull_white::finite_difference::examples::VanillaCall{ 80 });

	heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double> copy_of_rhs = rhs;

	a2_implicit_inverter.SolveSystemInPlace(rhs);

	// copy to change name (who cares about cost here?)
	heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double> lhs = rhs;

	// just make another copy, why not
	heston_hull_white::finite_difference::tests::A2Tester<heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> a_2_tester{ mean_volatility, kappa, sigma , mean_vol_index, discretization };



	heston_hull_white::finite_difference::examples::A2<heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> i_minus_theta_times_epsilon_times_a2 = a2;
	a_2_tester.IMinusThetaEpsilonA2(i_minus_theta_times_epsilon_times_a2, theta, epsilon_t);

	i_minus_theta_times_epsilon_times_a2.TransformGridInPlace(lhs);

	std::size_t num_ss = lhs.NumS();
	std::size_t num_vs = lhs.NumV();
	std::size_t num_rs = lhs.NumR();

	// Now require that lhs == rhs
	for (std::size_t i = 0; i < num_ss; ++i) {
		for (std::size_t j = 0; j < num_vs; ++j) {
			for (std::size_t k = 0; k < num_rs; ++k) {

				//std::cout << "i = " << i << " j = " << j << " k = " << k << "  " <<  lhs.EntryByIndex(i, j, k) << "     " << copy_of_rhs.EntryByIndex(i, j, k) << "\n";

				REQUIRE(lhs.EntryByIndex(i, j, k) == Approx(copy_of_rhs.EntryByIndex(i, j, k)).epsilon(.01));
			}
		}
	}

}