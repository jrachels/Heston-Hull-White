#ifndef HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A1TESTS_A1TESTER_H
#define HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A1TESTS_A1TESTER_H

#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace tests {
			template<typename GridType>
			class A1Tester {
			public:
				A1Tester(std::shared_ptr<const examples::ADIDiscretization> discretization) : a1_(discretization) {}


				std::array<double, 3> GetA1MatrixEntry(int i, int j, int k) {
					return a1_.tridiagonal_matrix_.EntryByIndex(i, j, k);
				}

				examples::ThreeDimensionalGrid<std::array<double, 3>> GetA1Matrix() {
					return a1_.tridiagonal_matrix_;
				}

				void IMinusThetaEpsilonA1(examples::A1<examples::ThreeDimensionalGrid<double>>& a_1_prime, double theta, double epsilon_t) {

					double multiplicative_constant = -theta * epsilon_t;
					//std::cout << "\n\n\nIMinusThetaEpsilonA1 tests: \n\n\n";
					std::size_t grid_size = a_1_prime.tridiagonal_matrix_.Size();
					for (std::size_t i = 0; i < grid_size; ++i) {
						auto [entry_1, entry_2, entry_3] = a_1_prime.tridiagonal_matrix_.EntryByVectorIndex(i);
						a_1_prime.tridiagonal_matrix_.RefToEntryByVectorIndex(i) = { multiplicative_constant * entry_1, multiplicative_constant * entry_2 + 1, multiplicative_constant * entry_3 };
						//std::cout << "i: " << i << "     " << (multiplicative_constant * entry_1) << "     " << (multiplicative_constant * entry_2 + 1) << "    " << (multiplicative_constant * entry_3) << "\n";
						if (std::abs(multiplicative_constant * entry_2 + 1) < std::abs(multiplicative_constant * entry_1) + std::abs(multiplicative_constant * entry_3)) {
							throw;
						}
					}

				}

				double GetG1Constant(int i, int j, int k) {
					return a1_.g_1_constants_.EntryByIndex(i, j, k);
				}

				examples::ThreeDimensionalGrid<double> GetG1() {
					return a1_.g_1_constants_;
				}

			private:
				examples::A1<GridType> a1_;

			};
		}
	}
}

#endif