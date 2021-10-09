#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_DOUGLAS_SRC_DOUGLAS_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_DOUGLAS_SRC_DOUGLAS_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/Douglas/douglas.h"

template<typename GridType>
GridType heston_hull_white::finite_difference::examples::Douglas<GridType>::Next(const GridType& current_prices) {
	current_time_ -= epsilon_;

	// TODO: Fix duct tape solution with current_time. On last loop, issues occur if current_time < 0.
	assert(current_time_ >= -(epsilon_ / 7.0));
	if (current_time_ < 0) {
		current_time_ = 0;
	}

	// construct an a3_
	std::unique_ptr<A3<GridType>> a3 = std::make_unique<A3<GridType>>(a3_builder_.GetA3(current_time_));
	A3ImplicitInversion<GridType> a3_implicit_inversion_(*a3.get(), theta_, epsilon_);

	// line 1
	GridType Y_0 = a0_.TransformGrid(current_prices);
	GridType Y_1 = a1_.TransformGrid(current_prices);
	GridType Y_2 = a2_.TransformGrid(current_prices);
	GridType Y_3 = previous_a3_->TransformGrid(current_prices);


	Y_0.Add(Y_1);
	Y_0.Add(Y_2);
	Y_0.Add(Y_3);
	a1_.ApplyG1InPlace(Y_0);
	a2_.ApplyG2InPlace(Y_0);
	Y_0.MultiplyByConstant<double>(epsilon_);
	Y_0.Add(current_prices); // Now equal to Y_0 from paper

	// line 2
	Y_1.MultiplyByConstant<double>(-epsilon_ * theta_);
	Y_1.Add(Y_0);
	a1_implicit_inversion_.SolveSystemInPlace(Y_1); // Now equal to Y_1 from paper

	Y_2.MultiplyByConstant<double>(-epsilon_ * theta_);
	Y_2.Add(Y_1);
	a2_implicit_inversion_.SolveSystemInPlace(Y_2); // Now equal to Y_2 from paper

	// line 3
	Y_3.MultiplyByConstant<double>(-epsilon_ * theta_);
	Y_3.Add(Y_2);
	a3_implicit_inversion_.SolveSystemInPlace(Y_3); // Now equal to Y_3 from paper

	// adjust a3 
	previous_a3_ = std::move(a3);

	return Y_3;

}

#endif