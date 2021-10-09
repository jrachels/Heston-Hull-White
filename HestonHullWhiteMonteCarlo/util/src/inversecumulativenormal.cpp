#include "HestonHullWhiteMonteCarlo/util/inversecumulativenormal.h"
#include <cmath>
#include <limits>

double internal::InverseCumulativeNormal::tail_value(double x) {
    // this part is unnecessary since general_canonical is guaranteed to work output a value between 0 and 1. - jeremy

    //if (x <= 0.0 || x >= 1.0) {
    //    // try to recover if due to numerical error
    //    if (close_enough(x, 1.0)) {
    //        return std::numeric_limits<double>::max(); // largest value available
    //    }
    //    else if (std::fabs(x) < std::numeric_limits<double>::epsilon()) {
    //        return std::numeric_limits<double>::min(); // largest negative value available
    //    }
    //    else {
    //        QL_FAIL("InverseCumulativeNormal(" << x
    //            << ") undefined: must be 0 < x < 1");
    //    }
    //}

    double z;
    if (x < x_low_) {
        // Rational approximation for the lower region 0<x<u_low
        z = std::sqrt(-2.0 * std::log(x));
        z = (((((c1_ * z + c2_) * z + c3_) * z + c4_) * z + c5_) * z + c6_) /
            ((((d1_ * z + d2_) * z + d3_) * z + d4_) * z + 1.0);
    }
    else {
        // Rational approximation for the upper region u_high<x<1
        z = std::sqrt(-2.0 * std::log(1.0 - x));
        z = -(((((c1_ * z + c2_) * z + c3_) * z + c4_) * z + c5_) * z + c6_) /
            ((((d1_ * z + d2_) * z + d3_) * z + d4_) * z + 1.0);
    }

    return z;
}


const double internal::InverseCumulativeNormal::a1_ = -3.969683028665376e+01;
const double internal::InverseCumulativeNormal::a2_ = 2.209460984245205e+02;
const double internal::InverseCumulativeNormal::a3_ = -2.759285104469687e+02;
const double internal::InverseCumulativeNormal::a4_ = 1.383577518672690e+02;
const double internal::InverseCumulativeNormal::a5_ = -3.066479806614716e+01;
const double internal::InverseCumulativeNormal::a6_ = 2.506628277459239e+00;

const double internal::InverseCumulativeNormal::b1_ = -5.447609879822406e+01;
const double internal::InverseCumulativeNormal::b2_ = 1.615858368580409e+02;
const double internal::InverseCumulativeNormal::b3_ = -1.556989798598866e+02;
const double internal::InverseCumulativeNormal::b4_ = 6.680131188771972e+01;
const double internal::InverseCumulativeNormal::b5_ = -1.328068155288572e+01;

const double internal::InverseCumulativeNormal::c1_ = -7.784894002430293e-03;
const double internal::InverseCumulativeNormal::c2_ = -3.223964580411365e-01;
const double internal::InverseCumulativeNormal::c3_ = -2.400758277161838e+00;
const double internal::InverseCumulativeNormal::c4_ = -2.549732539343734e+00;
const double internal::InverseCumulativeNormal::c5_ = 4.374664141464968e+00;
const double internal::InverseCumulativeNormal::c6_ = 2.938163982698783e+00;

const double internal::InverseCumulativeNormal::d1_ = 7.784695709041462e-03;
const double internal::InverseCumulativeNormal::d2_ = 3.224671290700398e-01;
const double internal::InverseCumulativeNormal::d3_ = 2.445134137142996e+00;
const double internal::InverseCumulativeNormal::d4_ = 3.754408661907416e+00;

// Limits of the approximation regions
const double internal::InverseCumulativeNormal::x_low_ = 0.02425;
const double internal::InverseCumulativeNormal::x_high_ = 1.0 - x_low_;