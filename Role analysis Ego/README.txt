These files contain R code to re-create the examples in the paper 
``Role Analysis in Networks using Mixtures of Exponential Random Graph Models''.

The foo_run.R files call other files to run code for the 3 examples shown in the paper. 
All code files are individually commented but here is a brief description of each:

lawyers_run.R         runs the Lazega lawyers example.
prosper_run.R         runs the prosper.com 2010 loans network example.
mrsim_run.R           runs many simulated data examples to explore bias and error in estimation.
sim_run.R             calls sim_networks.R and then fits EM using mixtures_of_ergms.R
rsim_run.R            runs the ego-ERGM EM fit for existing simulated networks. 
sim_networks.R        simulates networks to be used as ego-networks. 
likelihoods.R         calculates the loglikelihood and associated functions of a network given a set of parameters.
start.R               performs mixture model initialisation. 
terms_ego_ergm.R      calculates summary statistics. 
initialise.R          fits the two-stage model which is also used to initialise the ego-ERGM EM algorithm. 
mixtures_of_ergms.R   function to perform EM algorithm to find ego-ERGM parameters and memberships.
plot_ego_ergm.R       plots a whole network with nodes coloured according to egoERGM results.
regular_equiv.R       fits a regular equivalence model to a network. 
