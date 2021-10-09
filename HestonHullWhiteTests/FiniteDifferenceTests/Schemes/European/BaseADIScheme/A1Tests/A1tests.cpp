

#include "HestonHullWhiteTests/catch.hpp"
#include <memory>
#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"
#include "HestonHullWhiteTests/FiniteDifferenceTests/Schemes/European/BaseADIScheme/A1Tests/A1tester.h"

TEST_CASE("A1 Accuracy Tests") {

	double K = 110;
	double T = 5;

	heston_hull_white::finite_difference::examples::ADIDiscretizationInputs discretization_inputs{ \
	20 /* m1_ */, K / 20.0 /* d_1 */, std::max(.5, std::exp(-T / (4.0))) * K /* S_left */, K /* S_right */, 14 * K /* S_max*/, \
	20 /* m_2 */, .02 /*d_2 = V_max/500*/, 10 /* V_max */, \
	20 /* m_3 */, 0.0025 /* d_3 = R_max/400 */, .1 /* c */, 1 /* R_max */ };


	std::shared_ptr<heston_hull_white::finite_difference::examples::ADIDiscretization> discretization = \
		std::make_shared< heston_hull_white::finite_difference::examples::ADIDiscretization>(discretization_inputs );

	heston_hull_white::finite_difference::tests::A1Tester<heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> a_1_tester{ discretization };

	//std::cout << "\n\n\n\n A_TESTER RESULT " << a_1_tester.GetA1MatrixEntry(0, 0, 0)[1];

}