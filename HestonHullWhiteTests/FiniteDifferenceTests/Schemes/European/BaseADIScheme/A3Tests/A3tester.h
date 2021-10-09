#ifndef HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A3TESTS_A3TESTER_H
#define HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A3TESTS_A3TESTER_H

#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace tests {
			template<typename GridType>
			class A3Tester {
			public:
				A3Tester(examples::A3<GridType>& a3) : a3_(a3) {}

				std::array<double, 3> GetA3MatrixEntry(int k) {
					return a3_.condensed_a3_matrix_[k];
				}

				std::vector<std::array<double, 3>> GetA3Matrix() {
					return a3_.condensed_a3_matrix_;
				}


				void IMinusThetaEpsilonA3(examples::A3<examples::ThreeDimensionalGrid<double>>& a_3_prime, double theta, double epsilon_t) {

					double multiplicative_constant = -theta * epsilon_t;
					//std::cout << "\n\n\nIMinusThetaEpsilonA3 tests: \n\n\n";
					std::size_t grid_size = a_3_prime.condensed_a3_matrix_.size();
					for (std::size_t i = 0; i < grid_size; ++i) {
						auto [entry_1, entry_2, entry_3] = a_3_prime.condensed_a3_matrix_[i];
						a_3_prime.condensed_a3_matrix_[i] = { multiplicative_constant * entry_1, (multiplicative_constant * entry_2) + 1, multiplicative_constant * entry_3 };
						//std::cout << "i: " << i << "     " << (multiplicative_constant * entry_1) << "     " << (multiplicative_constant * entry_2 + 1) << "    " << (multiplicative_constant * entry_3) << "\n";
						if (std::abs(multiplicative_constant * entry_2 + 1) < std::abs(multiplicative_constant * entry_1) + std::abs(multiplicative_constant * entry_3)) {
							throw;
						}
					}
				}


			private:
				examples::A3<GridType> a3_;

			};
		}
	}
}

#endif