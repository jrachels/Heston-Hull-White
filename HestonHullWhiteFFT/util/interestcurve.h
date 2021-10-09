// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef HESTONHULLWHITEFFT_UTIL_INTERESTCURVE_INTERESTCURVE_H
#define HESTONHULLWHITEFFT_UTIL_INTERESTCURVE_INTERESTCURVE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <vector>

// TO DO: In the Interest curve constructor, Thomas's algorithm and inserting values into the matrices could be performed at the same time. This would reduce subtractions
// and memory accesses.

// some parts of this are borrowed from H. Kammeyer's InterestCurve implementation

namespace heston_hull_white {
	namespace fourier {
		namespace util {
			class InterestCurve {
			public:

				// assumes dates in discount curve are already ordered
				// assumes discount_curve[0][0] = 0 is the first starting date. 
				template<typename StdArrayLikeObject>
				InterestCurve(std::vector<StdArrayLikeObject> discount_curve);


				double ZeroBondCurve(double t) const;

				double Theta(double t, double lambda, double eta) const;


			private:
				// provides greatest index i such that dates_[i] <= date
				std::size_t IndexOfGreatestLowerBound_Date(double date) const;

				double ZeroBondToYield(double t, double ZeroBond) const;

				// discount_factor_[i] is the discount factor at date dates_[i]
				std::vector<double> dates_;
				//std::vector<double> discount_factor_;
				std::vector<std::array<double, 4>> discount_factor_spline_;
			};
		}
	}
}

template<typename StdArrayLikeObject>
heston_hull_white::fourier::util::InterestCurve::InterestCurve(std::vector<StdArrayLikeObject> discount_curve) : discount_factor_spline_(discount_curve.size()) {
	// This is probably the best place ot check to make sure dates_ are in order, etc. I will just assume they are.
		// double new_date = discount_curve[i][0]
		// if (new_date < previous date) throw 
		// dates_.push_back(new_date);
		// previous date = new_date;
	// I also assume that the discount factors are monotonic, I believe.
	size_t size = discount_curve.size();
	size_t size_minus_1 = size - 1;
	dates_.reserve(size);

	// constructing it with size is fine because std::array does not initialize values, I think.
	std::vector<std::array<double, 3>> intermediate_spline_matrix;
	intermediate_spline_matrix.reserve(size);
	intermediate_spline_matrix.push_back({ 0.0, 2.0, 0.0 });
	discount_factor_spline_[0][2] = 0.0;

	auto [current_date, current_discount_factor] = discount_curve[0];
	auto [next_date, next_discount_factor] = discount_curve[1];

	dates_.push_back(current_date);
	discount_factor_spline_[0][0] = current_discount_factor;


	for (std::size_t i = 1; i < size_minus_1; ++i) {

		double previous_date = current_date;
		double current_date = next_date;

		double previous_discount_factor = current_discount_factor;
		double current_discount_factor = next_discount_factor;
		next_date = discount_curve[i + 1][0];
		next_discount_factor = discount_curve[i + 1][1];


		intermediate_spline_matrix.push_back({ 1.0 / 3.0 * (current_date - previous_date), 2.0 / 3.0 * (next_date - previous_date), 1.0 / 3.0 * (next_date - current_date) });
		discount_factor_spline_[i][2] = ((next_discount_factor - current_discount_factor) / (next_date - current_date)) - ((current_discount_factor - previous_discount_factor) / (current_date - previous_date));


		dates_.push_back(current_date);
		discount_factor_spline_[i][0] = current_discount_factor;

		// cache the x's
		//intermediate_spline_matrix[i] = [1.0 / 3.0 * (x[i] - x[i - 1]), 2.0 / 3.0 * (x[i + 1] - x[i - 1]), 1.0 / 3.0 * (x[i + 1] - x[i])];

		//discount_factor_spline_[i][1] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);

	}

	dates_.push_back(next_date);
	discount_factor_spline_[size_minus_1][0] = next_discount_factor;

	intermediate_spline_matrix.push_back({ 0, 2.0, 0 });
	discount_factor_spline_[size_minus_1][2] = 0.0;

	// Thomas's Algorithm

	for (std::size_t i = 1; i < size; ++i) {
		double W = intermediate_spline_matrix[i][0] / intermediate_spline_matrix[i - 1][1];
		intermediate_spline_matrix[i][1] = intermediate_spline_matrix[i][1] - W * intermediate_spline_matrix[i - 1][2];
		discount_factor_spline_[i][2] = discount_factor_spline_[i][2] - W * discount_factor_spline_[i - 1][2];
	}

	discount_factor_spline_[size_minus_1][2] = discount_factor_spline_[size_minus_1][2] / intermediate_spline_matrix[size_minus_1][1];

	for (std::size_t i = size - 2; i > 0; --i) {
		discount_factor_spline_[i][2] = (discount_factor_spline_[i][2] - (intermediate_spline_matrix[i][2] * discount_factor_spline_[i + 1][2])) / (intermediate_spline_matrix[i][1]);
	}

	// end thomas's algorithm

	// set b and d

	for (std::size_t i = 0; i < size_minus_1; i++) {

		double epsilon_t = dates_[i + 1] - dates_[i];

		discount_factor_spline_[i][3] = (1.0 / 3.0) * ((discount_factor_spline_[i + 1][2] - discount_factor_spline_[i][2]) / (epsilon_t));
		discount_factor_spline_[i][1] = ((discount_factor_spline_[i + 1][0] - discount_factor_spline_[i][0]) / (epsilon_t)) - \
			((1.0 / 3.0) * (2.0 * discount_factor_spline_[i][2] + discount_factor_spline_[i + 1][2]) * (epsilon_t));

	}


	discount_factor_spline_[size_minus_1][3] = 0.0;

	double epsilon_t_last = dates_[size_minus_1] - dates_[size_minus_1 - 1];

	discount_factor_spline_[size_minus_1][1] = (3.0 * discount_factor_spline_[size_minus_1 - 1][3] * epsilon_t_last * epsilon_t_last) + \
		(2.0 * discount_factor_spline_[size_minus_1 - 1][2] * epsilon_t_last) + discount_factor_spline_[size_minus_1 - 1][1];

}




#endif