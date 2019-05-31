This code base will perform a series of learning experiments, and associated postprocessing, to produce the  results shown in figure 7 of  the paper [1]. The code supports other learning paradigms than those reported in the paper, but should be 'good to go', as supplied, with the published paradigm.

The code was originally written for batch processing on clusters using the Matlab Distributed Compute Engine (DCE) and associated toolbox; hence there are several levels of wrapping around the main learning code file - STDE_Shen_batch.m - to support this. The rationale here was to speed up the simulation  of a set of learning experiments, with different initial synaptic conductances,  and report averages (10 runs were used in Fig 7).

However, we assume most people won't be using the DCE...
To do a single experiment, it is best to proceed as follows:
(i) Run the top level script make_msns_in_channel.m with one of N_msns_D1, N_msns_D2 set to 1 (These are two parameters in the script, which is  currently set to use a single D1 MSN). Using this script ensures a correct call to the parameter construction function make_stde_pars.m, which returns a cell-array of parameter structures called all_pars. [The final line of make_msns_in_channel.m (currently commented out) calls the helper function do_batch.m which, in turn, calls STDE_Shen_batch.m.  To use the DCE, uncomment this line]

(ii) Call STDE_Shen_batch.m with the individual components (structures) of all_pars generated in (i). 
This will produce a results1.mat file which can be examined using the postprocessing tools described below.

Postprocessing
(a) To view the spike response profile (as in Figs 7A,D of [1]) run the function find_mean_spikes.m
Again, this was crafted with batch processing in mind and it assumes a set of results files results1.mat, results2.mat... For one run, put one of its paramters D1_xpts, D2_xpts to 1 , and the other to zero (according which type was run in (i) above).  The function parameter, markers,  is a vector of boundary trials between phases. It should reflect the values in the  array paramter trial_counts, defined in   make_msns_in_channel. Thus, using the values currently there, the markers vector would be set to [15, 55, 85, 115, 155].

(b) To see synaptic conductances in a view similar to that in Figs 7B,E, use find_mean_gs.m. For one experiment, xpt_nos will  be 1.
(c) To see synaptic conductances in a view similar to that in Figs 7C,F, use find_mean_gs_overall.m. Again, for one experiment, xpt_nos will just be 1.


[1] Gurney, K. N., Humphries, M. D., & Redgrave, P. (2015). A New Framework for Cortico-Striatal Plasticity: Behavioural Theory Meets In Vitro Data at the Reinforcement-Action Interface. PLoS Biology, 13(1), e1002034. http://doi.org/10.1371/journal.pbio.1002034

