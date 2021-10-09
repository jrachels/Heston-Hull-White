#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A3IMPLICITINVERSION_H
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_BASEADISCHEME_A3IMPLICITINVERSION_H

// A3ImplicitInversion implements the calculation of the partial terms in the r-direction the Heston-Hull-White PDE using an implicit scheme.

// TODO: currently I copy a3 and modify it to get the LU decomposition. I should probably modify a3 as I copy it.

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/A3.h"
#include "HestonHullWhiteFiniteDifference/GridExamples/ThreeDimensionalGrid/threedimensionalgrid.h"


namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			template<typename GridType>
			class A3ImplicitInversion {
			public:

				// performs the first half of reduction similar to Thomas's Algorithm

				A3ImplicitInversion(const A3<GridType>& a_3, double theta, double epsilon_t);

				// SolveSystem performs the second half of solving the linear system similar to Thomas's Algorithm

				// Solves the system
				// (I-theta*epsilon_t*A_3)Y_3=vec 
				// for Y_3.
				GridType SolveSystem(const GridType& vec);


				void SolveSystemInPlace(GridType& vec);

			private:
				std::vector <std::array<double, 3 >> condensed_LU_decomposition_;

			};

		}
	}
}

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/BaseADIScheme/src/A3ImplicitInversion.tpp"

#endif