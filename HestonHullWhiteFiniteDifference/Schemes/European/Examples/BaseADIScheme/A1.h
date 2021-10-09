#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A1_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A1_H

// A1 implements the calculation of the partial terms in the s-direction the Heston-Hull-White PDE using an explicit scheme.

// TODO: Make GridType a function template parameter instead?

#include <array>
#include <memory>
#include <vector>

#include <HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h>
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"


namespace heston_hull_white {
	namespace finite_difference {

		namespace tests {
			template<typename GridType> class A1Tester;
		}

		namespace examples {

			template<typename GridType>
			class A1 {

			public:
				friend class tests::A1Tester<GridType>;

				// TODO: perhaps I should just take a reference to the discretization?
				A1(std::shared_ptr<const ADIDiscretization> discretization);

				// Apply G1 adds in the constants produced by the finite difference scheme at the boundary.

				void ApplyG1InPlace(GridType& grid);

				// this could possibly be done more efficiently copying and then 
				GridType ApplyG1(const GridType& grid);

				// TransformGrid computes the sume of the s-direction partials at each point of the grid

				GridType TransformGrid(const GridType& vec);

				void TransformGridInPlace(GridType& vec);

				ThreeDimensionalGrid<std::array<double, 3>> GetMatrix() const {
					return tridiagonal_matrix_;
				}

			private:
				std::size_t num_ss_;
				std::size_t num_vs_;
				std::size_t num_rs_;
				ThreeDimensionalGrid<std::array<double, 3>> tridiagonal_matrix_;
				ThreeDimensionalGrid<double> g_1_constants_;
			};

		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A1.tpp"

#endif