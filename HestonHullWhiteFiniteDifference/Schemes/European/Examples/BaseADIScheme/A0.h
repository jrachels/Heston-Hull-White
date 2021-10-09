#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A0_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A0_H

// A0 implements the calculation of the mixed partial terms in the Heston-Hull-White PDE using an explicit scheme.

// TODO: The computation of A_0 can be simplified by factoring. This may slightly obfuscate the computation, however.
// Example: (s_beta_0 * r_beta_plus_1 * s_r_partial_coefficient) + (v_beta_0 * r_beta_plus_1 * v_r_partial_coefficient)
// could be rewritten as r_beta_plus_1 *((s_beta_0 * s_r_partial_coefficient) + (v_beta_0 *  v_r_partial_coefficient)).
// This removed a multiplication.

// TODO: Consider what order I should iterate over the variables in transform grid.

// I specifically skip transform in place because it is not as simple as all of the other matrices

#include <cmath>
#include <memory>

#include <HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h>
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"


namespace heston_hull_white {
	namespace finite_difference {
		namespace tests {
			template<typename GridType> class A0Tester;
		}

		namespace examples {

			template<typename GridType>
			class A0 {
			public:
				friend class tests::A0Tester<GridType>;

				// TODO: perhaps I should just take a reference to the discretization?
				A0(std::shared_ptr<const ADIDiscretization> discretization, std::size_t mean_vol_index /*included so I don't have to compute it twice*/, \
					const double vol_of_vol, const double eta, const double rho_s_v, const double rho_s_r, const double rho_v_r);

				// TransformGrid computes the sume of the mixed partials at each point of the grid

				// To help keep track of entries, use the keys:
				// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)
				// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)
				[[nodiscard]] GridType TransformGrid(const GridType& vec);

				// I specifically skip transform in place because it is not as simple as all of the other matrices


			private:
				std::size_t num_ss_;
				std::size_t num_vs_;
				std::size_t num_rs_;
				std::size_t mean_vol_index_;
				// Entries are coefficients to adjacent u's in computation of A_0 at u(s_i, v_j, r_k). Order of coefficients will be indicated by the shorthand ( a, b, c) for the coefficient
				// of u(s_{i+a}, v_{j+b}, r_{k+c}).
				// if v_index < mean_vol_index, entries follow this pattern
				// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1), (-1, -1, 0), (1, -1, 0), (-1, 1, 0), (1, 1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -1, -1), (0, 1, -1), (0, -1, 1), (0, 1, 1)
				// else entries follow this pattern
				// (0, 0, 0), (-1, 0, 0), (1, 0, 0), (0, -2, 0), (0, -1, 0), (0, 0, -1), (0, 0, 1), (-1, -2, 0), (1, -2, 0), (-1, -1, 0), (1, -1, 0), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (1, 0, 1), (0, -2, -1), (0, -1, -1), (0, -2, 1), (0, -1, 1)
				ThreeDimensionalGrid<std::array<double, 19>> sparse_A0_;
			};
		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A0.tpp"

#endif