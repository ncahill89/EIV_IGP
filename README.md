# EIV_IGP

⚠️ Known Bug Notice: A bug affecting temporal erros was identified in earlier versions of this repository.

✅ Fix: The issue has been resolved in the most recent version. Please download or clone the latest version from the repository.

ℹ️ Impact: Previous results should not be significantly affected, but I recommend updating to ensure correct behavior moving forward.


__About the Model__

The Errors in Variables Integrated Gaussian Process (EIV IGP) Model is used to perform Bayesian inference on historical rates of sea-level change. 

The input data to the model can be from tidegauge measurements and/or proxy reconstructions from cores of coastal sediment. These data are complicated by multiple sources of uncertainty, some of which arise as part of the data collection exercise. Notably, the proxy reconstructions include temporal uncertainty from dating of the sediment core using techniques such as radiocarbon. The EIV IGP model places a Gaussian process prior on the rate of sea-level change, which is then integrated to provide the mean of the likelihood for the observed data. The model is set in an errors-in-variables framework to take account of age uncertainty. The resulting model captures the continuous and dynamic evolution of sea-level change with full consideration of all available sources of uncertainty. For a more detailed model description check out the paper which you will find [here](https://www.jstor.org/stable/24522592?seq=4#metadata_info_tab_contents).

__Running the code__

Dowload the repo and click the `EIVIGP.Rproj` to open the project in Rstudio. Then you need to open the `runIGP.R` script to run the model. Example datasets are provided. 

__JAGS (Just Another Gibbs Sampler)__

Running these models requires JAGS (Just Another Gibbs Sampler) which you can install from [here](https://sourceforge.net/projects/mcmc-jags/).
