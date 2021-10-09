#ifndef HESTONHULLWHITEMONTECARLO_UTIL_INTERESTCURVE_INTERESTCURVE_H
#define HESTONHULLWHITEMONTECARLO_UTIL_INTERESTCURVE_INTERESTCURVE_H

#include <algorithm>
#include <cmath>
#include <vector>

// some parts of this are borrowed from H. Kammeyer's InterestCurve implementation

namespace heston_hull_white {
	namespace monte_carlo {
		namespace util {
			class InterestCurve {
			public:

				// assumes dates in discount curve are already ordered
				// assumes discount_curve[0][0] = 0 is the first starting date. 
				template<typename StdArrayLikeObject>
				InterestCurve(std::vector<StdArrayLikeObject> discount_curve) {
					//std::get<0, me>;
					size_t size = discount_curve.size();
					dates_.reserve(size);
					discount_factor_.reserve(size);
					//if (discount_curve[0][0] != 0) {

					//}

					// double previous date = 0;
					for (size_t i = 0; i < size; ++i) {
						// This is probably the best place ot check to make sure dates_ are in order, etc. I will just assume they are.
						// double new_date = discount_curve[i][0]
						// if (new_date < previous date) throw 
						// dates_.push_back(new_date);
						// previous date = new_date;
						dates_.push_back(discount_curve[i][0]);
						discount_factor_.push_back(discount_curve[i][1]);
					}
				}

				double ZeroBondCurve(double t) {
					const std::size_t num_dates_minus_one = dates_.size() - 1;
					size_t index_of_greatest_lower_bound = IndexOfGreatestLowerBound_Date(t);

					// extrapolate if past last date
					if (index_of_greatest_lower_bound == (num_dates_minus_one)) {
						return std::exp(-ZeroBondToYield(dates_[index_of_greatest_lower_bound], discount_factor_[index_of_greatest_lower_bound]) * t);
					}

					double X1 = dates_[index_of_greatest_lower_bound];
					double X2 = dates_[index_of_greatest_lower_bound + 1];
					double Y1 = discount_factor_[index_of_greatest_lower_bound];
					double Y2 = discount_factor_[index_of_greatest_lower_bound + 1];

					double ZBC = (Y2 - Y1) / (X2 - X1) * (t - X1) + Y1;
					return ZBC;
				}


				// copied from H. Kammeyer, removed an unnecessary negative
				double InstForward(double t) //instantaneous forward rates f(0,t) as defined in [BrMe], definition 1.4.2
				{
					const double epsilon = 1e-6; //"length" of numeric differential quotient
					double CV_1 = ZeroBondCurve(t);
					double CV_2 = ZeroBondCurve(t + epsilon);

					return (1.0 / CV_1) * ((CV_1 - CV_2) / epsilon); // use CV_1-CV_2 to remove the negative sign
				}

				// copied straight from H. Kammeyer
				double YieldCurve(double t) //zero coupon curve Y(0,t) as defined in [BrMe], definition 1.3.1, will be computed by linear interpolation
				{
					if (t < dates_[1]) //See documentation, section 4.4.1
					{
						//double X1 = itsDiscountCurve(2, 1);
						//double X2 = itsDiscountCurve(3, 1);
						//double Y1 = itsDiscountCurve(2, 2);
						//double Y2 = itsDiscountCurve(3, 2);
						double X1 = dates_[1];
						double X2 = dates_[2];
						double Y1 = discount_factor_[1];
						double Y2 = discount_factor_[2];

						double YieldY1 = ZeroBondToYield(X1, Y1);
						double YieldY2 = ZeroBondToYield(X2, Y2);

						return ((YieldY2 - YieldY1) / (X2 - X1)) * (t - X1) + YieldY1;
					}
					else //this case is not even needed at present
						return ZeroBondToYield(t, ZeroBondCurve(t));
				}

			private:
				// provides greatest index i such that dates_[i] <= date
				size_t IndexOfGreatestLowerBound_Date(double date) const {
					std::vector<double>::const_iterator first = dates_.begin();
					// I try to use standard library algorithms to be optimal. no alternative to upper_bound
					std::vector<double>::const_iterator iter_to_greatest_lower_bound = std::prev(std::upper_bound(first, dates_.end(), date));
					return std::distance(first, iter_to_greatest_lower_bound);
				}

				double ZeroBondToYield(double t, double ZeroBond) //conversion function
				{
					if (t < 1.0) //spot interest rates for maturities smaller than one year (e.g. LIBOR rates) are simply compounded rates.
						return (1.0 - ZeroBond) / (t * ZeroBond);
					else //annualy compounded rates
						return 1.0 / pow(ZeroBond, 1.0 / t) - 1.0;
				}

				// discount_factor_[i] is the discount factor at date dates_[i]
				std::vector<double> dates_;
				std::vector<double> discount_factor_;
			};
		}
	}
}

#endif