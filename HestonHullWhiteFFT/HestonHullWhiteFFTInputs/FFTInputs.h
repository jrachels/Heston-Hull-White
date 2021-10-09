// @author Jeremy Rachels
// Contact: jeremyrachels@gmail.com

#ifndef HESTONHULLWHITEFFT_HESTONHULLWHITEFFTINPUTS_FFTINPUTS_H
#define HESTONHULLWHITEFFT_HESTONHULLWHITEFFTINPUTS_FFTINPUTS_H

#include <cstddef>

namespace heston_hull_white {
	namespace fourier {

		struct FFTInputs {
			std::size_t number_grid_points_ = 4096;	//number of points of the discrete grid for Fourier transform
			const double integration_epsilon_ = 4.0 * 1e-3;		//size of integration step
			//const double lambda_;	//fineness of Fourier transform grid, calculated from other variables
			//const double b_; // max and min, calculated from other variables


			/*unsigned int FFT::N = 4096;
			double FFT::eta = 4.0 * 1e-3;
			double FFT::lambda = (2.0 * M_PI) / (FFT::N * FFT::eta);
			double FFT::b = M_PI / FFT::eta;*/

		};
	}
}

#endif