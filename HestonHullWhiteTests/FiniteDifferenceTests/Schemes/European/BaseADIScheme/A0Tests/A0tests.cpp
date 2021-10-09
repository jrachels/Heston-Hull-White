
//
//#include "HestonHullWhiteTests/catch.hpp"
//#include <memory>
//#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"
////#include "HestonHullWhiteTests/FiniteDifferenceTests/Schemes/European/BaseADIScheme/A1Tests/A1tester.h"
//
//TEST_CASE("A0 Accuracy Tests") {
//
//	double K = 110;
//	double T = 5;
//	double rho_s_v = -.3; /*rho_s_v*/
//	double rho_s_r = 0;
//	double rho_v_r = 0;
//	double eta = .2;
//	double mean_volatility = .02;
//	double vol_of_vol = .2;
//	
//
//	heston_hull_white::finite_difference::examples::ADIDiscretizationInputs discretization_inputs{ \
//	20 /* m1_ */, K / 20.0 /* d_1 */, std::max(.5, std::exp(-T / (4.0))) * K /* S_left */, K /* S_right */, 14 * K /* S_max*/, \
//	20 /* m_2 */, .02 /*d_2 = V_max/500*/, 10 /* V_max */, \
//	20 /* m_3 */, 0.0025 /* d_3 = R_max/400 */, .1 /* c */, 1 /* R_max */ };
//
//
//	std::shared_ptr<heston_hull_white::finite_difference::examples::ADIDiscretization> discretization = \
//		std::make_shared< heston_hull_white::finite_difference::examples::ADIDiscretization>(discretization_inputs);
//
//	heston_hull_white::finite_difference::tests::A1Tester<heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> a_1_tester{ discretization };
//
//	std::size_t mean_vol_index = discretization->GetIndexV(mean_volatility);
//
//	heston_hull_white::finite_difference::examples::A0< heston_hull_white::finite_difference::examples::ThreeDimensionalGrid<double>> a0{ discretization, mean_vol_index, vol_of_vol, eta, rho_s_v, rho_s_r, rho_v_r };
//
//
//
//	
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	std::cout << "\n\n\n\n A_TESTER RESULT " << a_1_tester.GetA1MatrixEntry(0, 0, 0)[1];
//
//}