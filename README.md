A Bayes Sampler by JWG SIMONS

This repository includes:
- A latex file titled "Bayesian Assignment JWG Simons". This latex file includes .  
- An annotated R file titled "Bayesian Assignment JWG Simons - Code". This is the code with which the reader can reproduce the results in the report.  
- This README file, with instructions on which pieces of code to run to produce the results. 

Please select and run the following lines of code to reproduce the results in the report: 

Preparing the data: 
- Running lines 10 - 19 installs the package with the data if necessary, loads the data, and executes the necessary data operations. 
- Running lines 22 - 29 grand-mean centers the variables, and gives a summary of the grand-centered data. 
- Running lines 32 - 37 executes an MLE regression analysis of the second model (with all variables included), and stores coefficients 
                        and standard deviations for later use in the MH step.

Estimation - Gibbs sampler & Metropolis-Hastings algorithm: 
- Running lines 44 - 131 defines the gibbs sampler function.
- Running lines 134 - 275 defines the Metropolis-Hastings algorithm. 
- Running line 278 executes the Gibbs sampler and stores the results.  
- Running line 280 executes the MH algorithm and stores the results.  

Assessing model convergence:
- Running lines 286 - 306 installs the required packages when necessary, and defines a function for plotting history plots. 
- Running lines 308 - 311 prints the history plots for the first model. 
- Running lines 313 - 317 prints the history plots for the second model.
- Running lines 322 - 351 defines a function for plotting auto-correlation plots.  
- Running lines 354 - 357 prints the autocorrelation plots for the first model. 
- Running lines 359 - 363 prints the autocorrelation plots for the second model.  
- Running lines 367 - 368 calculates the MC error for the first model on line 367 and the second model on line 368.
- Running lines 370 - 371 provides TRUE/FALE output for whether the MC error is smaller than 5% of sample standard deviation. 
- Running lines 374 - 375 prints the acceptance ratio of the MH step for the parameters "ads" and "premium". 

Posterior predictive checks:
- Running lines 382 - 386 combines the chains in each of the models. 
- Running lines 388 - 440 defines a function for computing the PPC for both models. 
- Running line 443 computes and prints the PPC for the first model, running line 444 computes and prints the PPC for the second model. 

Model comparison and selection with the DIC: 
- Running lines 451 - 498 defines a function for computing the DIC in both models. 
- Running line 500 computes the DIC for the first model, running line 501 computes the DIC for the second model.  
 
Hypothesis evaluation with the Bayes factor:  
- Lines 507 - 535 sample from the posterior and prior of the second model (please read the annotation for the specific steps)
- Lines 539 - 543 calculate the Bayes factor for the first hypothesis. 
- Lines 547 - 551 calculate the Bayes factor for the second hypothesis. 
- Lines 554 - 562 calculate the Bayes factor for the third hypothesis. 

Parameter estimates and credible intervals: 
- Lines 569 and 571 provide the EAP estimates and associated CI's for the second model, respectively. 
- Lines 573 - 590 define a function for plotting the posteriors of the parameters, with EAP estimates and CI's on the x-axis. It also loads the gridExtra package for plotting if necessary.
- Lines 592 - 595 plot the posterior distributions for each of the parameters in the second model. 
