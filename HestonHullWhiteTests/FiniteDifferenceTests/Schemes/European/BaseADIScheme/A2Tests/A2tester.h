#ifndef HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A2TESTS_A2TESTER_H
#define HESTONHULLWHITETESTS_FINITEDIFFERENCETESTS_SCHEMES_EUROPEAN_BASEADISCHEME_A2TESTS_A2TESTER_H

#include "HestonHullWhiteFiniteDifference/Include/hestonhullwhitefinitedifference.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace tests {
			template<typename GridType>
			class A2Tester {
			public:
				A2Tester(const double mean_volatility, const double kappa, const double sigma /*vol of vol*/, std::size_t mean_vol_index /*included so I don't have to compute it twice*/, std::shared_ptr<const examples::ADIDiscretization> discretization) : \
					a2_(mean_volatility, kappa, sigma, mean_vol_index, discretization) {}

				std::array<double, 5> GetA2MatrixEntry(int i, int j, int k) {
					return a2_.condensed_a2_matrix_.EntryByIndex(i, j, k);
				}

				examples::ThreeDimensionalGrid<std::array<double, 5>> GetA2Matrix() {
					return a2_.condensed_a2_matrix_;
				}

				double GetG2Constant(int j) {
					return a2_.g_2_constants_[j];
				}

				void IMinusThetaEpsilonA2(examples::A2<examples::ThreeDimensionalGrid<double>>& a_2_prime, double theta, double epsilon_t) {

					const double multiplicative_constant = -theta * epsilon_t;

					std::size_t grid_size = a_2_prime.condensed_a2_matrix_.Size();
					for (std::size_t i = 0; i < grid_size; ++i) {
						auto [entry_1, entry_2, entry_3, entry_4, entry_5] = a_2_prime.condensed_a2_matrix_.EntryByVectorIndex(i);
						a_2_prime.condensed_a2_matrix_.RefToEntryByVectorIndex(i) = { multiplicative_constant * entry_1, multiplicative_constant * entry_2, multiplicative_constant * entry_3 + 1, \
							multiplicative_constant * entry_4, multiplicative_constant * entry_5 };
					}

				}


				std::vector<double> GetG2() {
					return a2_.g_2_constants_;
				}

				std::size_t GetMeanVolIndex() {
					return a2_.mean_vol_index_;
				}

			private:
				examples::A2<GridType> a2_;

			};
		}
	}
}

#endif