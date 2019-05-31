% make MSNs in a single channel, both D1 and D2

% The results (one for each MSN) are stored in files name results#, where # is
% an expt number.
% These raw results files will then be used by code like make
% MAKE_MSN_RESPONSES

clear all;

% define the structure component designators of 
% all_pars so they can be altered in a sensible way

S_SIM_GENERAL_STRUCT    = 1;
S_XPT_STRUCTURE         = 2;
S_CTX_AND_STIM_STRUCT   = 3;
S_STDP_ELIGIBILITY      = 4;
S_DOPAMINE              = 5;
S_MSN                   = 6;
S_PATTERN               = 7;
              
% set one of the following to 1 for non-batched mode, and the
% other to zero.
N_msns_D1 = 1; % set to one for non-batched mode with D1 MSN
N_msns_D2 = 0;% set to one for non-batched mode with D2 MSN

% phases(1) = RANDOM_PATTERNS;
% phases(2) = PATTERN_DISCOVERY;
% phases(3) = PATTERN_MATCH_WITH_PHASIC_DA;
% phases(4) = RANDOM_PATTERNS;
% phases(5) = PATTERN_MATCH_WITH_DA_DIP;
% phases(6) = RANDOM_PATTERNS;

trial_counts(1) = 15;   
trial_counts(2) = 40;   
trial_counts(3) = 30;    
trial_counts(4) = 30;   
trial_counts(5) = 40; 
trial_counts(6) = 30; 


N = 200;                % number of synapses
rp = randperm(N);       % for determining teh set of strong afferents
% seed = 1;

% ================================================ %

%% set up parameters

xpt_no = 1;
neuron_type = 'D1';     
for i = 1:N_msns_D1
    seed = i;
    all_pars = make_stde_pars(neuron_type, trial_counts, N, rp, seed);
    all_pars{S_SIM_GENERAL_STRUCT}.xpt_no = xpt_no;
    pars{xpt_no} = all_pars;
    xpt_no = xpt_no + 1;
end

neuron_type = 'D2';     
for i = 1:N_msns_D2
    seed = i + N_msns_D2;
    all_pars = make_stde_pars(neuron_type, trial_counts, N, rp, seed);
    all_pars{S_SIM_GENERAL_STRUCT}.xpt_no = xpt_no;
    pars{xpt_no} = all_pars;
    xpt_no = xpt_no + 1;
end
No_xpts = xpt_no - 1;

% =============================================== %

%% batch mode stuff   
% do_batch(pars, 'STDE_Shen_batch', No_xpts);
% uncommnet if you are using the Matlab DCE with this helper function


