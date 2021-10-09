#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A1IMPLICITINVERSION_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A1IMPLICITINVERSION_H

// A1ImplicitInversion implements the calculation of the partial terms in the s-direction the Heston-Hull-White PDE using an implicit scheme.

// TODO: currently I copy a1 and modify it to get the LU decomposition. I should probably modify a1 as I copy it.

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A1.h"
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			template<typename GridType>
			class A1ImplicitInversion {
			public:

				// performs the first half of Thomas's Algorithm
				A1ImplicitInversion(const A1<GridType>& a_1, double theta, double epsilon_t);

				// Solve System performs the second half of Thomas's Algorithm.

				// Solves the system
				// (I-theta*epsilon_t*A_1)Y_1=Y_0-theta*epsilon_t*U_(n-1) 
				// for Y_1.
				GridType SolveSystem(const GridType& vec);

				void SolveSystemInPlace(GridType& vec);

			private:
				ThreeDimensionalGrid<std::array<double, 3>> LU_decomposition_;
			};

		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A1ImplicitInversion.tpp"

#endif