# Heston-Hull-White
This repository contains three stand-alone projects implementing three different algorithms for pricing options under the Heston-Hull-White model. This code was written as a private exercise in numerical methods for option pricing. Some sections exist only to test design concepts and language features but otherwise provide only nominal reductions in computation time and memory usage. Some features you may expect to be present are excluded. For example, the finite difference class provides no method for extrapolating prices between grid points.

## FFT

The folder

> Heston-Hull-White/HestonHullWhiteFFT/ 

contains an implementation of the algorithm from this paper

https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1399389

for pricing European call options. An example of how to use this method can be found in 

> Heston-Hull-White/HestonHullWhiteTests/FFTTests/EuropeanCallPricer/europeancallpricertests.cpp

## Monte Carlo

The folder

> Heston-Hull-White/HestonHullWhiteMonteCarlo/

contains an implementation of the Quadratic Exponential algorithm from this paper

https://papers.ssrn.com/sol3/papers.cfm?abstract_id=946405

for pricing derivatives under the HHW model. I have slightly modified the algorithm to account for nonzero correlation between interest rates and each/either of s or v.

The Monte Carlo method has been written to be able to price multiple path-dependent and path-independent options simultaneously. This complicates the process of executing the Monte Carlo method. Examples are provided in 

> Heston-Hull-White/HestonHullWhiteTests/MonteCarloTests/European/EuropeanMonteCarloTests.cpp

Step 1:

Define a ```StateType``` that models the evolution of parameters affecting the price of the underlying. The example I provide is ```HestonHullWhite``` in

> Heston-Hull-White/HestonHullWhiteMonteCarlo/StateTypes/Examples/HestonHullWhite/

Step 2:

Define classes representing the option types you want to price. Examples can be found in

> Heston-Hull-White/HestonHullWhiteMonteCarlo/OptionTypes/Examples/

These classes require a member function ```Price(const StateType& state) const;``` that should return a ```double```, an ```std::array<double>```, or an ```std::vector<double>``` representing the price(s) of the derivative(s) at the given state.

Step 3: 

Construct an ```std::shared_ptr< heston_hull_white::monte_carlo::util::InterestCurve >``` by calling the constructor of ```InterestCurve``` on an ```std::vector<std::array<double, 2>>``` of (date, discount_factor) pairs.


Step 4:

Instantiate a ```heston_hull_white::monte_carlo::HestonHullWhiteInputParameters```.


Step 5:

Instantiate a ```heston_hull_white::monte_carlo::PathIndependents``` with the path-independent options you are pricing.

Step 6:

Instantiate a ```heston_hull_white::monte_carlo::PathDependents``` with the path-dependent options you are pricing.

Step 7: 

Construct a ```heston_hull_white::monte_carlo::european::MultipleDerivativesSimulator``` or similar that can sample option prices.

Step 8: 

Construct a ```heston_hull_white::monte_carlo::european::MonteCarloSimulation```.

Step 9:

Access the result with the member function ```MeanOutcomes()```.


## Finite Difference

The folder

> Heston-Hull-White/HestonHullWhiteFiniteDifference/

contains an implementation of several finite difference algorithms for pricing European options under the HHW model. These are the Douglas, Craig-Sneyd, Modified Craig-Sneyd, and Hundsdorfer ```Schemes```. The paper on which these methods are based is given here:

https://arxiv.org/pdf/1111.4087.pdf .

Examples of how to use these methods can be found in 

> Heston-Hull-White/HestonHullWhiteTests/FiniteDifferenceTests/Schemes/European/

Recommended input parameters for instantiating ```heston_hull_white::finite_difference::examples::ADIDiscretizationInputs``` are provided at the top of page 5 of the linked ADI paper.

As stated previously, I provide no method for extrapolating prices between grid points. You will need to decide how you want to accomplish this. 

## Testing

Each of the above three projects is stand-alone. They can be compiled individually with a C++20 compatible compiler. Testing is provided by a fourth project found in the folder

> Heston-Hull-White/HestonHullWhiteTests/

You can run these tests with 

> Heston-Hull-White/x64/Debug/HestonHullWhiteTests.exe

These tests take about 10 minutes on my personal computer.


## References
<a id="1">[1]</a> 
Kienitz, J., Kammeyer, H.: An implementation of the Hybrid-Heston-Hull-White model. Available at SSRN 1399389 (2009)

<a id="2">[2]</a>
L. Andersen, *Simple and efficient simulation of the Heston stochastic volatility model*, J. Comput. Finance, 11 (2008), pp. 1–42.

<a id="3">[3]</a>
T. Haentjens and K. J. In ’t Hout, *ADI finite difference schemes for the Heston-Hull-White PDE*, arXiv:1111.4087, (2011).
 
