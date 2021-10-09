/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2010 Kakhkhor Abdijalilov

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef HESTONHULLWHITEMONTECARLO_UTIL_INVERSECUMULATIVENORMAL_H
#define HESTONHULLWHITEMONTECARLO_UTIL_INVERSECUMULATIVENORMAL_H


namespace internal {

    class InverseCumulativeNormal {
    public:
        InverseCumulativeNormal(double average = 0.0, double sigma = 1.0) : average_(average), sigma_(sigma) {};
        // function
        double operator()(double x) const {
            return average_ + sigma_ * standard_value(x);
        }
        // value for average=0, sigma=1
        /* Compared to operator(), this method avoids 2 floating point
           operations (we use average=0 and sigma=1 most of the
           time). The speed difference is noticeable.
        */
        static double standard_value(double x) {
            double z;
            if (x < x_low_ || x_high_ < x) {
                z = tail_value(x);
            }
            else {
                z = x - 0.5;
                double r = z * z;
                z = (((((a1_ * r + a2_) * r + a3_) * r + a4_) * r + a5_) * r + a6_) * z /
                    (((((b1_ * r + b2_) * r + b3_) * r + b4_) * r + b5_) * r + 1.0);
            }

            // The relative error of the approximation has absolute value less
            // than 1.15e-9.  One iteration of Halley's rational method (third
            // order) gives full machine precision.
            // #define REFINE_TO_FULL_MACHINE_PRECISION_USING_HALLEYS_METHOD
#ifdef REFINE_TO_FULL_MACHINE_PRECISION_USING_HALLEYS_METHOD
// error (f_(z) - x) divided by the cumulative's derivative
            const Real r = (f_(z) - x) * M_SQRT2 * M_SQRTPI * exp(0.5 * z * z);
            //  Halley's method
            z -= r / (1 + 0.5 * z * r);
#endif

            return z;
        }
    private:
        /* Handling tails moved into a separate method, which should
           make the inlining of operator() and standard_value method
           easier. tail_value is called rarely and doesn't need to be
           inlined.
        */
        static double tail_value(double x);
//#if defined(QL_PATCH_SOLARIS)
//        CumulativeNormalDistribution f_;
//#else
//        static const CumulativeNormalDistribution f_;
//#endif
        double average_, sigma_;
        static const double a1_;
        static const double a2_;
        static const double a3_;
        static const double a4_;
        static const double a5_;
        static const double a6_;
        static const double b1_;
        static const double b2_;
        static const double b3_;
        static const double b4_;
        static const double b5_;
        static const double c1_;
        static const double c2_;
        static const double c3_;
        static const double c4_;
        static const double c5_;
        static const double c6_;
        static const double d1_;
        static const double d2_;
        static const double d3_;
        static const double d4_;
        static const double x_low_;
        static const double x_high_;
    };
}

#endif