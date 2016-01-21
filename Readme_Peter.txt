fixes to be merged into official QuantLib
- (in master) gsr refined observability (evaluation date, termstructure)
- (in master) type Volatility vs Real in swaptionvolstructure
- (in master) proxynonstandardswaptionengine: struct/class state halper, !(integrationPoints_ ==0), nextExOrigDate is not used
- (todo) CapHelper for normal and sln vols
- (todo) Swaption instrument, extend impliedVolatility method to normal vols
- (todo) several compiler warnings about unused fields
- (todo) fix this:
../ql/experimental/models/mcgaussian1dnonstandardswaptionengine.hpp:87:10: warning: 'QuantLib::McGaussian1dNonstandardSwaptionEngine<QuantLib::GenericPseudoRandomMultiThreaded<QuantLib::MersenneTwisterMultiThreaded, QuantLib::InverseCumulativeNormal>, QuantLib::GenericRiskStatistics<QuantLib::GenericGaussianStatistics<QuantLib::GeneralStatistics> > >::reset' hides overloaded virtual function [-Woverloaded-virtual]
    void reset() const;
         ^
../ql/experimental/models/mcgaussian1dnonstandardswaptionengine.hpp:401:13: note: in instantiation of template class 'QuantLib::McGaussian1dNonstandardSwaptionEngine<QuantLib::GenericPseudoRandomMultiThreaded<QuantLib::MersenneTwisterMultiThreaded, QuantLib::InverseCumulativeNormal>, QuantLib::GenericRiskStatistics<QuantLib::GenericGaussianStatistics<QuantLib::GeneralStatistics> > >' requested here
        new McGaussian1dNonstandardSwaptionEngine<RNG, S>(
            ^
montecarlo_multithreaded.cpp:312:9: note: in instantiation of member function 'QuantLib::MakeMcGaussian1dNonstandardSwaptionEngine<QuantLib::GenericPseudoRandomMultiThreaded<QuantLib::MersenneTwisterMultiThreaded, QuantLib::InverseCumulativeNormal>, QuantLib::GenericRiskStatistics<QuantLib::GenericGaussianStatistics<QuantLib::GeneralStatistics> > >::operator shared_ptr' requested here
        MakeMcGaussian1dNonstandardSwaptionEngine<PseudoRandomMultiThreaded>(
        ^
../ql/pricingengine.hpp:68:14: note: hidden overloaded virtual function 'QuantLib::GenericEngine<QuantLib::NonstandardSwaption::arguments, QuantLib::NonstandardSwaption::results>::reset' declared here: different qualifiers (none vs const)
        void reset() { results_.reset(); }
             ^
1 warning generated.
- (in master) fix messages in lgm tests
- (todo) make the zerobond option method virtual in Gaussian1dModel (leave the default implementation as is, this can be used e.g. by the Markov functional model),  override the zerobond option method in Gsr with a closed form solution, adapt the gaussian1d jamshidian swaption engine such that it incorporates the transformation to a one curve swaption
