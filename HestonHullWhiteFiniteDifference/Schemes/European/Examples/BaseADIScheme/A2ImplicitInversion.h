#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A2IMPLICITINVERSION_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A2IMPLICITINVERSION_H

// A2ImplicitInversion implements the calculation of the partial terms in the v-direction the Heston-Hull-White PDE using an implicit scheme.

// TODO: currently I copy a2 and modify it to get the LU decomposition. I should probably modify a2 as I copy it.

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A2.h"
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			template<typename GridType>
			class A2ImplicitInversion {

			public:
				// performs the first half of reduction similar to Thomas's Algorithm

				A2ImplicitInversion(const A2<GridType>& a_2, double theta, double epsilon_t);

				// SolveSystem performs the second half of solving the linear system similar to Thomas's Algorithm

				// Solves the system
				// (I-theta*epsilon_t*A_2)Y_2=vec
				// for Y_2.
				GridType SolveSystem(const GridType& vec);

				void SolveSystemInPlace(GridType& vec);

			private:
				std::size_t num_vs_;
				std::size_t num_rs_;
				std::size_t mean_vol_index_;
				ThreeDimensionalGrid<std::array<double, 5>> LU_decomposition_;

			};

		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A2ImplicitInversion.tpp"

#endif