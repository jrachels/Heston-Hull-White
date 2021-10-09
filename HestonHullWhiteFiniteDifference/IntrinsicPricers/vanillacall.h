#ifndef HESTONHULLWHITEFINITEDIFFERENCE_INTRINSICPRICERS_VANILLACALL_H
#define HESTONHULLWHITEFINITEDIFFERENCE_INTRINSICPRICERS_VANILLACALL_H

#include <algorithm>

namespace heston_hull_white {
	namespace finite_difference {
		namespace examples {

			class VanillaCall {
			public:
				VanillaCall(const double strike) : strike_(strike) {}

				double operator() (double s, double v, double r) const {
					return std::max(s - strike_, 0.0);
				}

			private:
				const double strike_;

			};

		}
	}
}

#endif