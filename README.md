# Sensitivity Analyses of Variational Susceptibility COVID mMdel

This repository consists of the data and results obtained from modifying the code by Aguas et al. from the paper [Herd immunity thresholds for SARS-CoV-2 estimated from unfolding epidemics](https://www.medrxiv.org/content/10.1101/2020.07.23.20160762v2.full.pdf) 
 
Briefly, we tested changing the amount of time at maximum social distancing for 5,10,...,150 days. We also tested changing the social distancing curve to linearly decrease until it reaches 5,10,...,100% of the maximum social distancing at which it stays until the end of the simulation. Lastly, we did a sensitivty analysis on the rampdown of social distancing by testing 120,240,..., 1440 day long rampdowns. The way this was implemented was running a python script to generate a new file for each parameter change and running a bash script to run all the files on separate nodes of a supercomputer. This can be seen in the codegen files in the corresponding analysis folders.

All the inferred parameters for each test is in a corresponding folder in the results section. Additionally, we added data up to August for Belgium and England and these datasets are in the data folder.

To replicate our analyses, run the codegen file for the corresponding test you'd like to replicate. Then running all of the generated Parameter Estimation files will give you the parameters to put into the Epidemic matlab code and see Epidemic curves as well as Rc and Reff curves. Furthermore, plugging in the estimated R0 and CV into the formula for HIT in the paper above will give you the HIT for that parameter. Some graphs of the sensitivity of HIT to the parameters analyzed can be seen in the ResultsAnalysis python notebook in the results folder in each analysis folder.
