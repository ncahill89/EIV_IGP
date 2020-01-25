# EIV_IGP
The Errors in Variables Integrated Gaussian Process (EIV IGP) Model is used to perform Bayesian inference on historical rates of sea-level change. 

The input data to the model can be from tidegauge measurements and/or proxy reconstructions from cores of coastal sediment. These data are complicated by multiple sources of uncertainty, some of
which arise as part of the data collection exercise. Notably, the proxy reconstructions include temporal uncertainty from dating of the sediment core using
techniques such as radiocarbon. The EIV IGP model places a Gaussian process prior on the rate of sea-level change, which is then integrated to provide the mean of the likelihood for the observed data. The model is set
in an errors-in-variables framework to take account of age uncertainty. The resulting model captures the continuous and dynamic evolution of sea-level
change with full consideration of all available sources of uncertainty.
