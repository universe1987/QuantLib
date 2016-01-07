fixes to be merged into official QuantLib
- (in master) gsr refined observability (evaluation date, termstructure)
- (in master) type Volatility vs Real in swaptionvolstructure
- (in master) proxynonstandardswaptionengine: struct/class state halper, !(integrationPoints_ ==0), nextExOrigDate is not used
- (todo) CapHelper for normal and sln vols
- (todo) several compiler warnings about unused fields
- (in master) fix messages in lgm tests
- (todo) make the zerobond option method virtual in Gaussian1dModel (leave the default implementation as is, this can be used e.g. by the Markov functional model),  override the zerobond option method in Gsr with a closed form solution, adapt the gaussian1d jamshidian swaption engine such that it incorporates the transformation to a one curve swaption
