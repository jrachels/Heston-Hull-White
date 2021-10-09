#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A2_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A2_H

// A2 implements the calculation of the partial terms in the v-direction the Heston-Hull-White PDE using an explicit scheme.

// TODO: I did not enforce the boundary condition u(s, v_max, r, t) = s because I forgot. I don't know if this leads to instability

// TODO: this has to be a bit slower to  include the -ru term, since it changes the coefficients, meaning I have to save values for each r and v rather than just each v.
// Even if I don't change this, condensed_a2_matrix_ should probably be a two dimensional grid.

// TODO: condensed_a2_matrix_ could do with a more compact representation. Each row only has 3 or 4 nonzero entries, but I store 5 entries. I do this because I want to make finding the LU
// decomposition easier. Condensing the matrix could make the LU decomposition more complicated.

// TODO: If I am going to have five entries per row anyway, wouldn't it be better to just use second order central approximations over backwards approximations?

#include <cassert>
#include <memory>
#include <vector>

#include <HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h>
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace tests {
			template<typename GridType> class A2Tester;
		}

		namespace examples {

			template<typename GridType>
			class A2 {
			public:
				friend class tests::A2Tester<GridType>;

				// TODO: perhaps I should just take a reference to the discretization?
				A2(const double mean_volatility, const double kappa, const double sigma /*vol of vol*/, std::size_t mean_vol_index /*included so I don't have to compute it twice*/, std::shared_ptr<const ADIDiscretization> discretization);

				// Apply G2 adds in the constants produced by the finite difference scheme at the boundary.

				void ApplyG2InPlace(GridType& grid);

				// this could possibly be done more efficiently copying and then 
				void ApplyG2(const GridType& grid);

				// TransformGrid computes the sume of the v-direction partials at each point of the grid

				GridType TransformGrid(const GridType& vec);


				void TransformGridInPlace(GridType& vec);

				std::size_t MeanVolIndex() const {
					return mean_vol_index_;
				}

				ThreeDimensionalGrid<std::array<double, 5>> GetMatrix() const {
					return condensed_a2_matrix_;
				}

			private:
				std::size_t num_ss_;
				std::size_t num_vs_;
				std::size_t num_rs_;
				std::size_t mean_vol_index_;
				ThreeDimensionalGrid<std::array<double, 5>> condensed_a2_matrix_;
				std::vector<double> g_2_constants_;
			};
		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A2.tpp"

#endif