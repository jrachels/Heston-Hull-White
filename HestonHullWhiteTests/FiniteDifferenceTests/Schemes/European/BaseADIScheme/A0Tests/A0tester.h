#ifndef HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A0TESTS_A0TESTER_H
#define HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A0TESTS_A0TESTER_H

#include <memory>

#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace tests {
			template<typename GridType>
			class A0Tester {
			public:
				A0Tester(std::shared_ptr<const examples::ADIDiscretization> discretization, std::size_t mean_vol_index /*included so I don't have to compute it twice*/, \
					const double vol_of_vol, const double eta, const double rho_s_v, const double rho_s_r, const double rho_v_r) : \
					a0_(discretization, mean_vol_index, vol_of_vol, eta, rho_s_v, rho_s_r, rho_v_r) {}

				std::array<double, 19> GetA0MatrixEntry(int i, int j, int k) {
					return a0_.sparse_A0_.EntryByIndex(i, j, k);
				}

				examples::ThreeDimensionalGrid<std::array<double, 19>> GetA0Matrix() {
					return a0_.sparse_A0_;
				}

				std::size_t GetMeanVolIndex() {
					return a2_.mean_vol_index_;
				}

			private:
				examples::A0<GridType> a0_;

			};
		}
	}
}

#endif