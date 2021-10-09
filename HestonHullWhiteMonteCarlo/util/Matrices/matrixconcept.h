// TO DO: Matrix concept should allow you to specific the required member functions (accessor and size) and not apply specific names of them.
// even better if they can have common default names like size1 if the size or accessor functions aren't specified.

#ifndef HESTONHULLWHITEMONTECARLO_UTIL_MATRICES_MATRIXCONCEPT_H
#define HESTONHULLWHITEMONTECARLO_UTIL_MATRICES_MATRIXCONCEPT_H

#include <concepts>

namespace heston_hull_white {
	namespace util {

		template<typename T>
		concept MatrixConcept = requires (T a, int b, int c) { //intellisense issue
			{a(b, c)} -> std::convertible_to<double>;
			{a.rows()}->std::convertible_to<size_t>;
			{a.columns()}->std::convertible_to<size_t>;
		};

	}
}
#endif