#ifndef HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_HUNDSDORFERVERWER_SRC_HUNDSDORFERVERWER_TPP
#define HESTONHULLWHITEFINITEDIFFERENCE_SCHEMES_EUROPEAN_EXAMPLES_HUNDSDORFERVERWER_SRC_HUNDSDORFERVERWER_TPP

#include "HestonHullWhiteFiniteDifference/Schemes/European/Examples/HundsdorferVerwer/hundsdorferverwer.h"

template<typename GridType>
GridType heston_hull_white::finite_difference::examples::HundsdorferVerwer<GridType>::Next(const GridType& current_prices) {
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

	GridType Line_5_Intermediate = Y_0;

	a1_.ApplyG1InPlace(Y_0);
	a2_.ApplyG2InPlace(Y_0);
	Y_0.MultiplyByConstant<double>(epsilon_);
	Y_0.Add(current_prices); // Now equal to Y_0 from paper

	// line 2
	Y_1.MultiplyByConstant<double>(-epsilon_ * theta_);
	//GridType Y_1_tilde = Y_1;
	Y_1.Add(Y_0);
	a1_implicit_inversion_.SolveSystemInPlace(Y_1); // Now equal to Y_1 from paper

	Y_2.MultiplyByConstant<double>(-epsilon_ * theta_);
	//GridType Y_2_tilde = Y_2;
	Y_2.Add(Y_1);
	a2_implicit_inversion_.SolveSystemInPlace(Y_2); // Now equal to Y_2 from paper

	// line 3
	Y_3.MultiplyByConstant<double>(-epsilon_ * theta_);
	//GridType Y_3_tilde = Y_3;
	Y_3.Add(Y_2);
	a3_implicit_inversion_.SolveSystemInPlace(Y_3); // Now equal to Y_3 from paper

	// line 4
	GridType Y_0_tilde = a0_.TransformGrid(Y_3);
	GridType A_1_Y_3 = a1_.TransformGrid(Y_3);
	GridType A_2_Y_3 = a2_.TransformGrid(Y_3);
	GridType A_3_Y_3 = a3->TransformGrid(Y_3);
	Y_0_tilde.Add(A_1_Y_3);
	Y_0_tilde.Add(A_2_Y_3);
	Y_0_tilde.Add(A_3_Y_3);
	Y_0_tilde.Subtract(Line_5_Intermediate);
	Y_0_tilde.MultiplyByConstant<double>((.5) * (epsilon_));
	Y_0_tilde.Add(Y_0);

	// line 5
	A_1_Y_3.MultiplyByConstant<double>(-(epsilon_) * (theta_));
	Y_0_tilde.Add(A_1_Y_3);
	a1_implicit_inversion_.SolveSystemInPlace(Y_0_tilde);


	A_2_Y_3.MultiplyByConstant<double>(-(epsilon_) * (theta_));
	Y_0_tilde.Add(A_2_Y_3);
	a2_implicit_inversion_.SolveSystemInPlace(Y_0_tilde);

	// line 6
	A_3_Y_3.MultiplyByConstant<double>(-(epsilon_) * (theta_));
	Y_0_tilde.Add(A_3_Y_3);
	a3_implicit_inversion_.SolveSystemInPlace(Y_0_tilde);

	// adjust a3 
	previous_a3_ = std::move(a3);

	return Y_0_tilde;

}

#endif