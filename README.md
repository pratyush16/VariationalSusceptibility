# VariationalSusceptibility

This repository consists of the data and results obtained from modifying the code by Gomes et al. from the paper [Herd immunity thresholds for SARS-CoV-2 estimated from unfolding epidemics](https://www.medrxiv.org/content/10.1101/2020.07.23.20160762v2.full.pdf) 
 
We tested changing the amount of time at maximum social distancing for 5,10,...,150 days. We also tested changing the social distancing curve to linearly decrease until it reaches 5,10,...,100% of the maximum social distancing at which it stays until the end of the simulation. Lastly, we did a simple test changing the rampdown to 1000 days and the total simulation to 1500 days. The way this was implemented was running a python script to generate a new file for each parameter change and running a bash script to run all the files on separate nodes of a supercomputer. This can be seen in the codegen files in the corresponding results folders.

All the inferred parameters for each is in a corresponding folder in the results section. Additionally, we added data up to August for Belgium and England and these datasets are in the data folder.
