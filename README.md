# COVID19distancing

This R code accompanies a Conversation article on the effects of spatial and temporal physical distancing in reducing the rate of COVID-19 infection (authors: Mike Lee, Corey Bradshaw, Craig Dalton)

Contact: Professor Corey J. A. Bradshaw, Global Ecology, College of Science and Engineering, GPO Box 2100, Flinders University, Adelaide, South Australia 5001, Australia

corey.bradshaw@flinders.edu.au
+61 (0) 400 697 665

Cellular automaton model to investigate reductions in COVID-19 infection rate with spatial and temporal physical distancing

This R code produces an R×C matrix populated by N individuals.

For each cell, individual is resampled from a maximum of half the matrix dimension, with each movement conditioned a binomial probability set by the user.
 
The process can be run over t time steps and s iterations. 

Note, there is no recovery component included.
 
The user can define two additional parameters:
 
1. A reduction in movement probability (spatial distancing)
2. A reduction of in exposure (temporal distancing)

