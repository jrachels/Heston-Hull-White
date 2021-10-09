#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A3_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A3_H

// A3 implements the calculation of the partial terms in the r-direction the Heston-Hull-White PDE using an explicit scheme.

#include <memory>
#include <vector>

#include <HestonHullWhiteFiniteDifference/DiscretizationExamples/ADIDiscretization/ADIDiscretization.h>
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"

namespace heston_hull_white {
	namespace finite_difference {

		namespace tests {
			template<typename GridType> class A3Tester;
		}

		namespace examples {

			template<typename GridType>
			class A3 {
				friend class tests::A3Tester<GridType>;
			public:
				// TODO: perhaps I should just take a reference to the discretization?
				A3(const double lambda, const double eta, const double theta_t, std::shared_ptr<const ADIDiscretization> discretization);

				// TransformGrid computes the sume of the r-direction partials at each point of the grid

				// TODO: This should be faster if the loops were:
				// Outer: r
				// Inner: s
				// Inner Inner: v
				// possibly switching s and v based on memory considerations
				GridType TransformGrid(const GridType& vec);

				void TransformGridInPlace(GridType& vec);

				std::vector<std::array<double, 3 >> GetMatrix() const {
					return condensed_a3_matrix_;
				}

			private:

				std::size_t num_ss_;
				std::size_t num_vs_;
				std::size_t num_rs_;
				std::vector<std::array<double, 3>> condensed_a3_matrix_;

			};
		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A3.tpp"

#endif