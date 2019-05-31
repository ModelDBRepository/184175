function all_pars = make_stde_pars(neuron_type, trial_counts, N, randp, seed)

% Set up STDE parameters for use by STDE_Shen
%
% neuron_type is 'D1' or 'D2'
% 
% trial_counts is a vector of trial counts in each phase of the experiment
% It should have 5 components: 
% #trials with random patterns
% #trials with DA burst
% #trials with random patterns
% #trials with DA dip
% #trials with random patterns
%
% N is the number of synapses
%
% randp is a random permutaion of the N synapses describing where the
% strong afferents oocur
%
% seed is random seed
%

% ===  experimental structure and timing ====

s_xpt_structure = struct();

%%% possible phases of a trial 

% for replicating in vivo stim xpts
NO_STIM_NO_PHASIC_DA = 1;               % Slow wave only

WITH_STIM_NO_PHASIC_DA = 2;             % Slow wave with electrical stim

WITH_STIM_WITH_PHASIC_DA = 3;           % no slow wave ('activataed state') with 
                                        % electrical stim and phasic DA
                                    
NO_STIM_WITH_PHASIC_DA = 4;             % no slow wave ('activataed state') with 
                                        % phasic DA

WITH_STIM_WITH_REDUCED_PHASIC_DA = 5;   % no slow wave ('activataed state') with 
                                        % electrical stim and phasic DA

% for pattern matching xpts
PATTERN_MATCH_WITH_PHASIC_DA = 6;       % background cortical ctrivity, then 
                                        % salience on strong afferents ,
                                        % then background. Phasic DA burst after
                                        % salience
                                        
PATTERN_MATCH_WITH_DA_DIP = 7;          % background cortical ctrivity, then 
                                        % salience on strong afferents ,
                                        % then background. Phasic DA dip after
                                        % salience
                                        
RANDOM_PATTERNS = 8;                    % background cortical ctrivity, then 
                                        % salience on random afferents ,
                                        % then background
PATTERN_DISCOVERY = 9;                  % same as #6, but strong afferents gradually
                                        % become more likely and phasic DA
                                        % occurs stchatically dependent on
                                        % number of strong afferents

% phases
% phases(1) = RANDOM_PATTERNS;
% phases(2) = PATTERN_DISCOVERY;
% phases(3) = PATTERN_MATCH_WITH_PHASIC_DA;
% phases(4) = RANDOM_PATTERNS;

phases(1) = RANDOM_PATTERNS;
phases(2) = PATTERN_MATCH_WITH_PHASIC_DA;
phases(3) = RANDOM_PATTERNS;
phases(4) = PATTERN_MATCH_WITH_DA_DIP;
phases(5) = RANDOM_PATTERNS;


% Number of trial in each interval
xpt_intervals(1) = trial_counts(1);  % 15
xpt_intervals(2) = trial_counts(2); % 20 match with s_pattern.tau_da_habituate
xpt_intervals(3) = trial_counts(3); % 10 
xpt_intervals(4) = trial_counts(4); % 30
xpt_intervals(5) = trial_counts(5); % 30


s_xpt_structure.phases = phases;
s_xpt_structure.intervals = xpt_intervals;
s_xpt_structure.half_period = 1.0;          %  pre- and post stim time
s_xpt_structure.DA_time = 1.55;             %  time of occurence of start of DA pulse -
                                            % needs to be greater than half_period
s_xpt_structure.stim_duration = 0.015;

s_xpt_structure.init_strong_frac_std = 0.25;  % initial std dev of fraction of
                                              % strong afferents for
                                              % PATTERN DISCOVERY
                                              
% ==============  end of experimental structure and timing ==============

% ==============  general simulation parameters ==============

s_sim_general_struct = struct();

s_xpt_structure.neuron_type = neuron_type;
s_sim_general_struct.dt = 0.0001; % 
s_sim_general_struct.rseed = seed;
s_sim_general_struct.xpt_no = 1;
s_sim_general_struct.spike_diag = 0;         % for recording afferent spikes
s_sim_general_struct.diag_elig_synapse = 0;  % set to zero to do no diagnostic here
s_sim_general_struct.diag_g_syn = 1;         % for recording g_syn

% ==============  end general simulation parameters ===========

% ============== cortex and stimulus parameters ==========================

s_ctx_and_stim_struct = struct();

s_ctx_and_stim_struct.background_rate = 8;  % background cortical firing rate when no slow wave

%% slow wave
s_ctx_and_stim_struct.freq_slow_wave = 0.8;
s_ctx_and_stim_struct.A_depth = 0.7;
s_ctx_and_stim_struct.sw_alpha_star = 0.05;
s_ctx_and_stim_struct.mean_sw_rate = 5;
s_ctx_and_stim_struct.sw_stim_reset_period_type = 1;
s_ctx_and_stim_struct.half_period_frac = 0.8;       % fraction of half period of slow wave 
                                                    % compelted at stim reset

s_ctx_and_stim_struct.stim_rate = 55;             % rate used to produce spikes during stim period
s_ctx_and_stim_struct.correlation = 0.15;

% ============== end of cortex and stimulus parameters ======================

% ============== STDP and eligibility ================================

s_stdp_elgibility = struct();

s_stdp_elgibility.tau_elig_pos = 0.3;
s_stdp_elgibility.tau_elig_neg = 0.3;
s_stdp_elgibility.tau_pos = 0.02;
s_stdp_elgibility.tau_neg = 0.02;
  
if neuron_type == 'D1'
    s_stdp_elgibility.k_hat_pos_hi = 1.3; %1.3
    s_stdp_elgibility.k_hat_pos_lo = -0.4; % -0.4
    s_stdp_elgibility.k_hat_neg_hi = 0.0;
    s_stdp_elgibility.k_hat_neg_lo = -0.5; % -0.475
    s_stdp_elgibility.LR = 0.65; % 0.65     
    
elseif neuron_type == 'D2'
    s_stdp_elgibility.k_hat_pos_hi = 0.35; % 0.3
    s_stdp_elgibility.k_hat_pos_lo = 0.3; % 0.3
    s_stdp_elgibility.k_hat_neg_hi = -0.85; % -0.8
    s_stdp_elgibility.k_hat_neg_lo = 0.3; % 0.3
    s_stdp_elgibility.LR = 0.65; % 0.65   
else
    error(make_stde_pars:Invalidneuron_type', 'Invalid neuron type for making paramters: must be D1 or D2');
end

% ============== end of STDP and eligibility ============================

% =============== dopamine ===============================

s_dopamine = struct();

LINEAR = 1;
NAKA_RUSHTON = 2;

s_dopamine.phasic = 20;
s_dopamine.tonic = 3.0;   % 3.0 
s_dopamine.tau_phi = 0.02;           % DA dynamics time constant
s_dopamine.DA_lo = 0.0;

s_dopamine.DA_std = 0.55; % 0.25

s_dopamine.reduced_phasic_factor = 0.1;

s_dopamine.DA_hi = 21;      % used in linear relation between alpha and DA only

if strcmp(neuron_type, 'D1');
    s_dopamine.NK_rho = 1.2;    
    s_dopamine.NK_theta = 6.0; 
    s_dopamine.NK_max = 1.2;
elseif strcmp(neuron_type, 'D2')
    s_dopamine.NK_rho = 1.4; %    
    s_dopamine.NK_theta = 1.8; 
    s_dopamine.NK_max = 1.0;    
else
    error(make_stde_pars:Invalidneuron_type', 'Invalid neuron type for making paramters: must be D1 or D2');
end

s_dopamine.DA_fn_type = NAKA_RUSHTON;

%% now do relation between DA and phi factor for full DA model
% assume Naka_Rusthon relation between phi and DA
s_dopamine.phi_max = 1.2;
s_dopamine.phi_rho = 1.8;
s_dopamine.phi_theta = 4.5;


% =============== end dopamine ===============================

% =============== msn ===============================

s_msn = struct();

s_msn.N         = N;    % Number of cortical afferents (usually 200)

s_msn.strong_perm  = randp;  % describes the selection of strong afferents

s_msn.syn_scaling = 0.4;        % scales all synaptic conductances from original DA_Iz model
s_msn.std_syn_frac = 0.1;       % std dev of synaptic strength at start id std_syn_frac * mean
s_msn.init_min_syn_frac = 0.5;  % minimum initial synaptic strength as fraction of mean
s_msn.max_syn_frac = 5;       % maximum synaptic weight as fraction of mean


% =============== end msn ===============================

% ================== pattern matching ===============

s_pattern = struct();

s_pattern.No_strong_affs  = 50;
s_pattern.rate_hi         = 25;      % rate of firing on the strong afferents
s_pattern.rate_lo         = 3.0;       % rate of firing on the weak afferents
s_pattern.alpha_hi        = 0.1;       % correlation between strong
                                     % afferents
s_pattern.alpha_lo        = 0;       % correlation between weak
                                     % afferents
s_pattern.salience_time   = 0.4;     % duration of salient input
s_pattern.tau_da_habituate = 20;     % (10) time constant for
                                      % habituation of phasic DA (in trial counts) 
% ================== end pattern matching =============


all_pars = {s_sim_general_struct, s_xpt_structure, s_ctx_and_stim_struct, ...
                        s_stdp_elgibility, s_dopamine, s_msn, s_pattern};

