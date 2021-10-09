#ifndef HESTONHULLWHITEFINITEDIFFERENCE_DISCRETIZATIONEXAMPLES_ADIDISCRETIZATION_ADIDISCRETIZATION_H
#define HESTONHULLWHITEFINITEDIFFERENCE_DISCRETIZATIONEXAMPLES_ADIDISCRETIZATION_ADIDISCRETIZATION_H


// The ADIDiscretization class contains caches values related to the chosen ADI discretization of the grid.


// If one discretization is going to be used in multiple schemes, precompute alphas, betas, gammas, and deltas

// TO DO: reduce multiplications with alphas, betas, gammas, and deltas by not calling the corresponding function
// and instead doing an online update of the values? May cause stability issues.

// TO DO: Decide on order of IntrinsicOptionValues

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			namespace internal {
				class CompareADICoordinates {
				public:
					bool operator()(const std::array<const double, 7>& lhs, const double rhs) const {
						return lhs[0] < rhs;
					}

					bool operator()(const std::array<const double, 10>& lhs, const double rhs) const {
						return lhs[0] < rhs;
					}
				};

			}


			struct ADIDiscretizationInputs {
				// constants for discretization in s direction
				std::size_t m_1;
				double d_1;
				double S_left;
				double S_right;
				double S_max;

				// constants for discretization in v direction
				std::size_t m_2;
				double d_2;
				double V_max;

				// constants for discretization in r direction
				std::size_t m_3;
				double d_3;
				double c;
				double R_max;
			};

			// In this example, ThreeDimensionalGrid is the model GridType.
			class ADIDiscretization {
			public:
				ADIDiscretization(const ADIDiscretizationInputs& inputs) : num_ss_(inputs.m_1), num_vs_(inputs.m_2), num_rs_(inputs.m_3 + 1), comparer{}{
					// I don't know if these are the absolute minimums. Using weird values (too small, or with mean_vol too close to zero) of num_vs_ and mean_vol will definitely break A2.
					assert(num_ss_ > 0);
					assert(num_vs_ > 2);
					assert(num_rs_ > 1);

					ConstructSDiscretization(inputs);

					ConstructVDiscretization(inputs);

					ConstructRDiscretization(inputs);
				
				}


				std::array<std::size_t, 3> GetIndex(double S, double V, double R) const {
					std::size_t S_index = std::distance(ss_.begin(), std::lower_bound(ss_.begin(), ss_.end(), S, comparer));
					std::size_t V_index = std::distance(vs_.begin(), std::lower_bound(vs_.begin(), vs_.end(), V, comparer));
					std::size_t R_index = std::distance(rs_.begin(), std::lower_bound(rs_.begin(), rs_.end(), R, comparer));
					return { S_index, V_index, R_index };
				}

				std::array<double, 3> ValueByIndex(int S_index, int V_index, int R_index) const {
					return { ss_[S_index][0], vs_[V_index][0], rs_[R_index][0] };
				}
			
				// order should be changed to improve memory
				template <typename IntrinsicValuePolicyType>
				ThreeDimensionalGrid<double> IntrinsicOptionValues(const IntrinsicValuePolicyType& intrinsic_pricer) const {

					ThreeDimensionalGrid<double> option_values{ num_ss_, num_vs_, num_rs_ };

					for (std::size_t i = 0; i < num_ss_; ++i) {
						for (std::size_t j = 0; j < num_vs_; ++j) {
							for (std::size_t k = 0; k < num_rs_; ++k) {
								option_values.EntryByIndex(i, j, k) = intrinsic_pricer(ss_[i][0], vs_[j][0], rs_[k][0]);
							}
						}
					}
					return option_values;
				}


				std::array<const double, 7> SCoefficientsByIndex(std::size_t index) const {
					return ss_[index];
				}

				std::array<const double, 10> VCoefficientsByIndex(std::size_t index) const {
					return vs_[index];
				}

				std::array<const double, 7> RCoefficientsByIndex(std::size_t index) const {
					return rs_[index];
				}

				std::size_t NumS() const {
					return num_ss_;
				}

				std::size_t NumV() const {
					return num_vs_;
				}

				std::size_t NumR() const {
					return num_rs_;
				}

				std::size_t NumTriples() const {
					return num_ss_ * num_vs_ * num_rs_;
				}

				size_t GetIndexV(double V) const {
					size_t V_index = std::distance(vs_.begin(), std::lower_bound(vs_.begin(), vs_.end(), V, comparer));
					return V_index;
				}

			private:
				std::size_t num_ss_;
				std::size_t num_vs_;
				std::size_t num_rs_;

				// 0 is boundary condition, so there is no entry
				// at 0<s<s_max, entries are s, central scheme coefficients for first partial, central scheme coefficients for second partial
				// at s=s_max, entries are s, 0, 0, 0, special central scheme coefficients for second partial
				std::vector<std::array<const double, 7>> ss_;


				// NOTE: The discretization is independent of eta, so I need to compute both alphas and betas for all v
				// at v =0, entries are v, forward scheme coefficients for first partial, 0, 0, 0, 0, 0, 0
				// // at v = v_1, entries are v, 0, 0, 0, central scheme coefficients for first partial, central scheme coefficients for second partial
				// at v_1 < v < v_max, entries are v, backward scheme coefficients for first partial, central scheme coefficients for first partial, central scheme coefficients for second partial
				// at v=V_max represents a boundary condition, and is not present in the following array
				std::vector<std::array<const double, 10>> vs_;

				// at r = -R_max, entires are r, 0, 0, 0, special central scheme coefficients for second partial
				// at -R_max<r<R_max, entires are r, central scheme coefficients for first partial, central scheme coefficients for second partial
				// at r = R_max, entires are r, 0, 0, 0, special central scheme coefficients for second partial
				std::vector<std::array<const double, 7>> rs_;

				internal::CompareADICoordinates comparer;

				void ConstructSDiscretization(const ADIDiscretizationInputs& inputs);

				void ConstructVDiscretization(const ADIDiscretizationInputs& inputs);

				void ConstructRDiscretization(const ADIDiscretizationInputs& inputs);

				double Phi(double xi, double S_left, double S_right, double d_1, double xi_int) {
					if (xi < 0) {
						return S_left + d_1 * std::sinh(xi);
					}
					else if (xi > xi_int) {
						return S_right + d_1* std::sinh(xi-xi_int);
					}
					else {
						return S_left + d_1 * xi;
					}
				}


				double AlphaMinus2(double delta_previous, double delta_current) {
					return delta_current / (delta_previous * (delta_previous + delta_current));
				}

				double AlphaMinus1(double delta_previous, double delta_current) {
					return (-delta_previous - delta_current) / (delta_previous * delta_current);
				}

				double Alpha0(double delta_previous, double delta_current) {
					return (delta_previous + 2.0 * delta_current) / (delta_current * (delta_previous + delta_current));
				}

				double BetaMinus1(double delta_current, double delta_next) {
					return (-delta_next) / (delta_current * (delta_current + delta_next));
				}

				double Beta0(double delta_current, double delta_next) {
					return (delta_next - delta_current) / (delta_next * delta_current);
				}

				double BetaPlus1(double delta_current, double delta_next) {
					return (delta_current / (delta_next * (delta_current + delta_next)));
				}

				double Gamma0(double delta_next, double delta_next_2) {
					return (-2.0 * delta_next - delta_next_2) / (delta_next * (delta_next + delta_next_2));
				}

				double GammaPlus1(double delta_next, double delta_next_2) {
					return (delta_next + delta_next_2) / (delta_next * delta_next_2);
				}

				double GammaPlus2(double delta_next, double delta_next_2) {
					return (-delta_next) / (delta_next_2 * (delta_next + delta_next_2));
				}

				double DeltaMinus1(double delta_current, double delta_next) {
					return 2.0 / (delta_current * (delta_current + delta_next));
				}

				double Delta0(double delta_current, double delta_next) {
					return -2.0 / (delta_current * delta_next);
				}

				double DeltaPlus1(double delta_current, double delta_next) {
					return 2.0 / (delta_next * (delta_current + delta_next));
				}

			};

		}
	}
}

#endif