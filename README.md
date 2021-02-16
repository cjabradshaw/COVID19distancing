# COVID19distancing

This <a href="https://github.com/cjabradshaw/COVID19distancing/blob/master/COVID19%20distancing.R">R code</a> accompanies a <em>The Conversation</em> article on the effects of spatial and temporal physical distancing in reducing the rate of COVID-19 infection: '<a href="https://theconversation.com/want-to-make-social-distancing-even-more-effective-its-about-time-as-well-as-space-134551">Want to make social distancing even more effective? It’s about time (as well as space</a>'

Authors:
- <a href="https://www.flinders.edu.au/people/mike.lee">Mike Lee</a>
- <a href="http://www.flinders.edu.au/people/corey.bradshaw">Corey Bradshaw<a/>
- <a href="https://www.newcastle.edu.au/profile/craig-dalton">Craig Dalton</a>

For code queries, contact: Professor Corey J. A. Bradshaw, <a href="https://globalecologyflinders.com">Global Ecology</a>, College of Science and Engineering, GPO Box 2100, Flinders University, Adelaide, South Australia 5001, Australia

corey.bradshaw@flinders.edu.au
+61 (0) 400 697 665

Cellular automaton model to investigate reductions in COVID-19 infection rate with spatial and temporal physical distancing

This R code produces an R×C matrix populated by <em>N</em> individuals.

For each cell, individual is resampled from a maximum of half the matrix dimension, with each movement conditioned a binomial probability set by the user.
 
The process can be run over <em>t</em> time steps and <em>s</em> iterations. 

Note, there is no recovery component included.
 
The user can define two additional parameters:
 
1. A reduction in movement probability (spatial distancing)
2. A reduction of in exposure (temporal distancing)

