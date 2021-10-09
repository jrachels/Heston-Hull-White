#include "HestonHullWhiteMonteCarlo/util/Matrices/matrix.h"

#include <cmath>

void heston_hull_white::util::Matrix::ScaleBy(double scale) {
	std::size_t size = mat_.size();
	for (std::size_t i = 0; i < size; i++) {
		double temp = mat_[i];
		mat_[i] = temp * scale;
	}
	return;
}

double heston_hull_white::util::Matrix::MaxAbsOfEntries() const {
	double max = 0;
	for (double entry : mat_) {
		if (std::abs(entry) > max) {
			max = entry;
		}
	}
	return max;
}

heston_hull_white::util::Matrix heston_hull_white::util::Matrix::Transpose() const {
	Matrix transpose(width_, height_);
	// should test if I should width or height first for speed
	for (int i = 0; i < width_; ++i) {
		for (int j = 0; j < height_; ++j) {
			transpose(i, j) = (this)->operator()(j, i);
		}
	}
	return transpose;
}


