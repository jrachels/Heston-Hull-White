#include "HestonHullWhiteTests/catch.hpp"

#include <algorithm>

#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"

TEST_CASE("ADIDiscretizationAccuracyTests") {

	double K = 110;
	double T = 5;

	heston_hull_white::finite_difference::examples::ADIDiscretizationInputs discretization_inputs{ \
	20 /* m1_ */, K / 20.0 /* d_1 */, std::max(.5, std::exp(-T / (4.0))) * K /* S_left */, K /* S_right */, 14 * K /* S_max*/, \
	20 /* m_2 */, .02 /*d_2 = V_max/500*/, 10 /* V_max */, \
	20 /* m_3 */, 0.0025 /* d_3 = R_max/400 */, .1 /* c */, 1 /* R_max */ };


	heston_hull_white::finite_difference::examples::ADIDiscretization discretization{ discretization_inputs };

	// conditions on S_coordinates
	REQUIRE(discretization.SCoefficientsByIndex(0)[0] == Approx(34.3036).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(1)[0] == Approx(47.8999).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(2)[0] == Approx(54.392).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(3)[0] == Approx(59.688).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(4)[0] == Approx(64.98286).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(5)[0] == Approx(70.27748046).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(10)[0] == Approx(96.750569).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(13)[0] == Approx(112.73632).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(14)[0] == Approx(120.9755).epsilon(0.1));
	REQUIRE(discretization.SCoefficientsByIndex(19)[0] == Approx(1540).epsilon(0.1));


	// Conditions on V_coordinates
	REQUIRE(discretization.VCoefficientsByIndex(0)[0] == 0);
	REQUIRE(discretization.VCoefficientsByIndex(1)[0] == Approx(.00704592).epsilon(0.001));
	REQUIRE(discretization.VCoefficientsByIndex(2)[0] == Approx(.0149408).epsilon(0.001));
	REQUIRE(discretization.VCoefficientsByIndex(3)[0] == Approx(.0246357).epsilon(0.001));
	REQUIRE(discretization.VCoefficientsByIndex(19)[0] == Approx(7.07945).epsilon(0.001));

	// Conditions on R_coordinates
	REQUIRE(discretization.RCoefficientsByIndex(0)[0] == -1);
	REQUIRE(discretization.RCoefficientsByIndex(1)[0] == Approx(-0.464027).epsilon(0.001));
	REQUIRE(discretization.RCoefficientsByIndex(2)[0] == Approx(-.189204).epsilon(0.001));
	REQUIRE(discretization.RCoefficientsByIndex(3)[0] == Approx(-.0482825).epsilon(0.001));
	REQUIRE(discretization.RCoefficientsByIndex(13)[0] == Approx(.108201).epsilon(0.001));
	REQUIRE(discretization.RCoefficientsByIndex(20)[0] == Approx(1).epsilon(0.0001));

	// Checking Alphas, Betas, Gammas, and Deltas

	// S_index = 0
	REQUIRE(discretization.SCoefficientsByIndex(0)[1] == Approx(-.00844461).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(0)[2] == Approx(-.04292011).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(0)[3] == Approx(.0513159).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(0)[4] == Approx(.00121021).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(0)[5] == Approx(-.004288).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(0)[6] == Approx(.00299197).epsilon(0.05));

	// S_index = 1
	REQUIRE(discretization.SCoefficientsByIndex(1)[1] == Approx(-.02297476).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(1)[2] == Approx(-.0819489).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(1)[3] == Approx(.104923673).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(1)[4] == Approx(.00707721).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(1)[5] == Approx(-.0222014).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(1)[6] == Approx(.015124).epsilon(0.05));

	// S_index = 7
	REQUIRE(discretization.SCoefficientsByIndex(7)[1] == Approx(-.09444119).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(7)[2] == Approx(0).margin(.0001));
	REQUIRE(discretization.SCoefficientsByIndex(7)[3] == Approx(.09444119).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(7)[4] == Approx(.03567655).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(7)[5] == Approx(-.0713531).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(7)[6] == Approx(.03567655).epsilon(0.05));


	// S_index = 19
	REQUIRE(discretization.SCoefficientsByIndex(19)[1] == 0);
	REQUIRE(discretization.SCoefficientsByIndex(19)[2] == 0);
	REQUIRE(discretization.SCoefficientsByIndex(19)[3] == 0);
	REQUIRE(discretization.SCoefficientsByIndex(19)[4] == Approx(.0000012799579).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(19)[5] == Approx(-.0000025591587).epsilon(0.05));
	REQUIRE(discretization.SCoefficientsByIndex(19)[6] == Approx(.0000012799579).epsilon(0.05));

	// V_index = 0
	REQUIRE(discretization.VCoefficientsByIndex(0)[1] == Approx(-208.857).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[2] == Approx(268.591).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[3] == Approx(-59.7335495).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[4] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[5] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[6] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[7] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[8] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(0)[9] == Approx(0).epsilon(0.05));

	// V_index = 1
	REQUIRE(discretization.VCoefficientsByIndex(1)[1] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[2] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[3] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[4] == Approx(-74.9951).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[5] == Approx(15.26173777).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[6] == Approx(59.7335495).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[7] == Approx(18998.4616).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[8] == Approx(-35953.962).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(1)[9] == Approx(16955.50035).epsilon(0.05));

	// V_index = 2
	REQUIRE(discretization.VCoefficientsByIndex(2)[1] == Approx(74.99527).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[2] == Approx(-268.5905).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[3] == Approx(193.5952).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[4] == Approx(-69.8131755).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[5] == Approx(23.51735).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[6] == Approx(46.29582).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[7] == Approx(14402.04139).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[8] == Approx(-26130.103).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(2)[9] == Approx(11728.06203).epsilon(0.05));

	// V_index = 7
	REQUIRE(discretization.VCoefficientsByIndex(7)[1] == Approx(24.5735).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[2] == Approx(-72.33973).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[3] == Approx(47.766228).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[4] == Approx(-17.628085).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[5] == Approx(8.693948).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[6] == Approx(8.934137).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[7] == Approx(757.43326).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[8] == Approx(-1296.65588).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(7)[9] == Approx(539.2226).epsilon(0.05));

	// V_index = 19

	REQUIRE(discretization.VCoefficientsByIndex(19)[1] == Approx(0.39999815).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[2] == Approx(-1.16683166).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[3] == Approx(0.766833).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[4] == Approx(-0.283179).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[5] == Approx(0.141253).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[6] == Approx(0.1419257).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[7] == Approx(0.1939219).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[8] == Approx(-0.33120808).epsilon(0.05));
	REQUIRE(discretization.VCoefficientsByIndex(19)[9] == Approx(0.13728615).epsilon(0.05));

	// R_index = 0
	REQUIRE(discretization.RCoefficientsByIndex(0)[1] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(0)[2] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(0)[3] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(0)[4] == Approx(3.481094).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(0)[5] == Approx(-6.9621887).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(0)[6] == Approx(3.481094).epsilon(0.05));

	// R_index = 1
	REQUIRE(discretization.RCoefficientsByIndex(1)[1] == Approx(-.6324148).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(1)[2] == Approx(-1.77291).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(1)[3] == Approx(2.405325).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(1)[4] == Approx(4.60231).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(1)[5] == Approx(-13.57787).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(1)[6] == Approx(8.97556).epsilon(0.05));

	// R_index = 11
	REQUIRE(discretization.RCoefficientsByIndex(11)[1] == Approx(-333.55121).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(11)[2] == Approx(162.6434).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(11)[3] == Approx(170.90781).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(11)[4] == Approx(273290.6).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(11)[5] == Approx(-468915.7).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(11)[6] == Approx(195625).epsilon(0.05));

	// R_index = 19
	REQUIRE(discretization.RCoefficientsByIndex(19)[1] == Approx(-2.932347).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(19)[2] == Approx(2.158425).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(19)[3] == Approx(0.7739217).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(19)[4] == Approx(13.373713).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(19)[5] == Approx(-20.24428).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(19)[6] == Approx(6.870571).epsilon(0.05));

	// R_index = 20

	REQUIRE(discretization.RCoefficientsByIndex(20)[1] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(20)[2] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(20)[3] == Approx(0).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(20)[4] == Approx(5.2001188).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(20)[5] == Approx(-10.4002376).epsilon(0.05));
	REQUIRE(discretization.RCoefficientsByIndex(20)[6] == Approx(5.2001188).epsilon(0.05));

	//std::cout << "START OF S ENTRIES \n\n";
	//for (int i = 0; i < 20; ++i) {
	//	for (double d : discretization.SCoefficientsByIndex(i)) {
	//		std::cout << d << "   ";
	//	}
	//	std::cout << "\n";
	//}




	//std::cout << "\n\n START OF V ENTRIES \n\n";
	//for (int j = 0; j < 20; ++j) {
	//	for (double d : discretization.VCoefficientsByIndex(j)) {
	//		std::cout << d << "   ";
	//	}
	//	std::cout << "\n";
	//}
	//std::cout << "\n\n START OF R ENTRIES \n\n";
	//for (int k = 0; k < 21; ++k) {
	//	for (double d : discretization.RCoefficientsByIndex(k)) {
	//		std::cout << d << "   ";
	//	}
	//	std::cout << "\n";
	//}



}