function  STDE_Shen_batch(s_sim_general_struct, s_xpt_structure, ...
    s_ctx_and_stim_struct, s_stdp_elgibility, s_dopamine, s_msn, s_pattern)

% Batched response for many single MSNs to user defined sequence of inputs over a series of trials
%
% The structures in teh parameter list will usually be produced using
% MAKE_STDE_PARS
% More details of those parameters are given there.
%
% A cell array of spike times - one vector in the array
% for each trial is saved to the
% results file 'RESULTS#', where # is the value of the parameter 'XPT_NO'
% (one of the general_sim parameters)
%
% Also saved in RESULTS# is STRONG_AFF_INDS, a row vector of indices
% specifying the strong afferents to the MSN
%
% Further variables are recorded depending on setting of flags.
%
% Thus, if  SPIKE_DIAG is set then all the afferent
% spikes are saved in the array T_SPIKES_SS with
% dimesnions N_synapses x N_timesteps
%
% If DIAG_G_SYN is set,  arrays of synaptic conductances for the entire
% experiment are saved in G_SYN_SS (for ampa) and G_SYN_NMDA_SS, G_SYN_GABA_SS (nmda and gaba resp)
% with dimension N_trials x N_synapses. The
% g_syn are recorded at the ned of each trial
%
% If DIAG_ELIG_SYNAPSE is set to a non-zero value N, then the eligibility
% traces and kernel values for synapes N are recorded. if N = 0, then no
% recording takes place.
% Eligibility traces are recorded in DIAG_POS_ELIGS_SS and DIAG_NEG_ELIGS_SS
% kernesl are stored in DIAG_POS_PLAS_SS and DIAG_NEG_PLAS_SS. These all
% have dimension of N_trials x N_timesteps
%
% Note taht some diagnostics can have huge memory demands - use carefully!
%
% If no saving to a results file is required, tehn comment out thesave line
% at the end of the file
%

% make all_pars for saving later
all_pars = {s_sim_general_struct, s_xpt_structure, ...
    s_ctx_and_stim_struct, s_stdp_elgibility, s_dopamine, s_msn, s_pattern};

%% experimental structure
neuron_type     = s_xpt_structure.neuron_type;
phases          = s_xpt_structure.phases;
xpt_ints        = s_xpt_structure.intervals;
half_period     = s_xpt_structure.half_period;
stim_duration   = s_xpt_structure.stim_duration;
time_DA         = s_xpt_structure.DA_time;

%% general sim pars
dt                  = s_sim_general_struct.dt;
seed                = s_sim_general_struct.rseed;
xpt_no              = s_sim_general_struct.xpt_no;
spike_diag          = s_sim_general_struct.spike_diag;
diag_elig_synapse   = s_sim_general_struct.diag_elig_synapse;
diag_g_syn          = s_sim_general_struct.diag_g_syn;

% stimulus and cortex
background_rate             = s_ctx_and_stim_struct.background_rate;

% slow wave
freq_slow_wave              = s_ctx_and_stim_struct.freq_slow_wave;
A_depth                     = s_ctx_and_stim_struct.A_depth;
sw_alpha_star               = s_ctx_and_stim_struct.sw_alpha_star;
mean_sw_rate                = s_ctx_and_stim_struct.mean_sw_rate;
sw_stim_reset_period_type   = s_ctx_and_stim_struct.sw_stim_reset_period_type;
half_period_frac            = s_ctx_and_stim_struct.half_period_frac;

stim_rate                   = s_ctx_and_stim_struct.stim_rate;
correlation_in_stim         = s_ctx_and_stim_struct.correlation;

% STDP and eligibility
tau_elig_pos    = s_stdp_elgibility.tau_elig_pos;
tau_elig_neg    = s_stdp_elgibility.tau_elig_neg;
tau_pos         = s_stdp_elgibility.tau_pos;
tau_neg         = s_stdp_elgibility.tau_neg;
k_hat_pos_hi    = s_stdp_elgibility.k_hat_pos_hi;
k_hat_pos_lo    = s_stdp_elgibility.k_hat_pos_lo;
k_hat_neg_hi    = s_stdp_elgibility.k_hat_neg_hi;
k_hat_neg_lo    = s_stdp_elgibility.k_hat_neg_lo;
LR              = s_stdp_elgibility.LR;

%% dopamine
da_tonic        = s_dopamine.tonic;

da_phasic               = s_dopamine.phasic;
reduced_phasic_factor   = s_dopamine.reduced_phasic_factor;
tau_phi                 = s_dopamine.tau_phi;
DA_hi                   = s_dopamine.DA_hi;
DA_lo                   = s_dopamine.DA_lo;
DA_std                  = s_dopamine.DA_std;
DA_fn_type              = s_dopamine.DA_fn_type;
NK_rho                  = s_dopamine.NK_rho;
NK_theta                = s_dopamine.NK_theta;
NK_max                  = s_dopamine.NK_max;

LINEAR = 1;
NAKA_RUSHTON = 2;

phi_max = s_dopamine.phi_max;
phi_rho = s_dopamine.phi_rho;
phi_theta = s_dopamine.phi_theta;

%%% MSN
N                   = s_msn.N;                  % Number of cortical afferents
syn_scaling         = s_msn.syn_scaling;        % scales all synaptic conductances from original DA_Iz model
std_syn_frac        = s_msn.std_syn_frac;       % std dev of synaptic strength at start id std_syn_frac * mean
init_min_syn_frac   = s_msn.init_min_syn_frac;  % minimum initial synaptic strength as fraction of mean
max_syn_frac        = s_msn.max_syn_frac;       % maximum synaptic weight as fraction of mean

%tau_syn  = s_msn.tau_syn;               % synaptic time constant
% V_syn    = s_msn.V_syn;                 % synaptic reversal potential

randp    = s_msn.strong_perm;           % sets strong afferent set

%% pattern matching input
No_strong_affs  = s_pattern.No_strong_affs;
rate_hi         = s_pattern.rate_hi;    % rate of firing on the strong afferents
rate_lo         = s_pattern.rate_lo;    % rate of firing on the weak afferents
alpha_hi        = s_pattern.alpha_hi;   % correlation between strong
% afferents
alpha_lo        = s_pattern.alpha_lo;   % correlation between weak
% afferents
salience_time   = s_pattern.salience_time;      % duration of salientinput
tau_da_habituate = s_pattern.tau_da_habituate;  % time constant for
% habituation of phasic DA


%% uses seconds and Volts

rand('state', seed);
randn('state', seed);

%%% possible phases of a trial
NO_STIM_NO_PHASIC_DA = 1;
WITH_STIM_NO_PHASIC_DA = 2;
WITH_STIM_WITH_PHASIC_DA = 3;
NO_STIM_WITH_PHASIC_DA = 4;
WITH_STIM_WITH_REDUCED_PHASIC_DA = 5;
PATTERN_MATCH_WITH_PHASIC_DA = 6;
PATTERN_MATCH_WITH_DA_DIP = 7;
RANDOM_PATTERNS = 8;
PATTERN_DISCOVERY = 9;
PATTERN_MATCH_WITH_SALIENCE_DECAY = 10;

%%% Experiment structure
No_phases = length(phases);
for i = 1:No_phases
    xpt_phases(1,i) = phases(i);
    xpt_phases(2,i) = xpt_ints(i);
end
trial_counts = xpt_phases(2,:);


%%% New DA Iz model pars %%%%%%%%%%%%

%% temporary assignmnst while debugging code
% D1 = 0.8;
% D2 = 0.8;

%load new_Iz_DA_pars
% MS neuron parameters in saved file
%k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);

k0 = 1.0;
a = 0.01;
b = -20.0;
c = -55.0;
vr0 = -80.0;
vpeak = 40.0;

% found MS parameters: X = [C,vt,d]
%C = X(1); vt =X(2); d = X(3);
C = 15.2294;
vt = -29.7303;
d0 = 90.9096;


% extra DA model parameters in saved file
% KIR = XD1(1);    % KIR modifier
% LCA = XD1(2);    % LCA modifier
KIR = 0.0289;
LCA = 0.3308;

DA_conc = da_tonic;
phi_nk1 = DA_conc .^ phi_rho;
D1 =  phi_max .* phi_nk1 ./ (phi_nk1 + phi_theta .^ phi_rho);
D2 = D1;

vrD1 = vr0 .* (1 + D1 .* KIR);
dD1 = d0 .* (1 - D1 .* LCA);

% D2 - intrinsic
% alpha = XD2;

alpha = 0.032;
kD2 = k0 .* (1 - alpha .* D2);

% synaptic
% cD1 = Xd1all;
% cD2 = Xd2all;

cD1 = 6.3;
cD2 = 0.215;

% all PSP parameters in saved file
Egaba = -60;
Enmda = 0;
Eampa = 0;

% these should stay in the same ratio
% PSPampa = Xsyn; %%

PSPampa = 6.8687; %
ampa_nmda = 2.0;
ampa_gaba = 1.4;

ts_ampa = 6.0;
ts_nmda = 160.0;
ts_gaba = 4.0;

PSPnmda = PSPampa ./ ampa_nmda;
PSPgaba = PSPampa ./ ampa_gaba;

PSPampa = PSPampa ./ ts_ampa; % normalisation used inline in original code
PSPnmda = PSPnmda ./ ts_nmda; % ditto
PSPgaba = PSPgaba ./ ts_gaba; % ditto


D1_syn_fact = (1 + cD1 .* D1);
D2_syn_fact = (1 - cD2 .* D2);

Mg = 1.0;

ms_dt = 1000 .* dt; % dt im milliseconds
dt_over_C = ms_dt ./ C; % because membrame equation assumes ms

if strcmp(neuron_type, 'D1')
    vr = vrD1;
    d = dD1;
    k = k0;
elseif strcmp(neuron_type, 'D2')
    vr = vr0;
    d = d0;
    k = kD2;
else
    error('STDP_Shen_batch:neuron_type', 'Invalid neurontype (D1 or D2)')
end

SynExp_ampa = exp(-ms_dt / ts_ampa);
SynExp_nmda = exp(-ms_dt / ts_nmda);
SynExp_gaba = exp(-ms_dt / ts_gaba);


%%%%%%%%% general parameters


%%% spikes and stimulus

rates_stim = [background_rate stim_rate background_rate];          % mean firing rate during stim trials
rates_no_stim = [background_rate background_rate background_rate];        % mean firing rate during non-stim trials

periods = [half_period stim_duration half_period];            % phases of cortical activity within a trial
alphas_stim = [0 correlation_in_stim 0];         % correlation in each phase for stim trials
alphas_no_stim = [0 0.0 0];      % correlation in each phase for non-stim trials


%%% synaptic conductances

g_syn_ampa_mean = syn_scaling .* PSPampa;
std_ampa = std_syn_frac .* g_syn_ampa_mean;

g_syn_ampa = g_syn_ampa_mean .* ones(1, N) + std_ampa .* randn(1, N);

g_syn_ampa(g_syn_ampa <= 0) = init_min_syn_frac .* g_syn_ampa_mean;   % don't let any synapses be less than init_min_g_syn

live_syn_ampa = true(1,N);   % live synapses - synapses that are reduced below zero
% are set to sero and are killed off. They
% can no longer be altered

max_syn_ampa = g_syn_ampa_mean .* max_syn_frac;

g_syn_nmda_mean = syn_scaling .* PSPnmda;
std_nmda = std_syn_frac .* g_syn_nmda_mean;
g_syn_nmda = g_syn_nmda_mean .* ones(1, N) + std_nmda .* randn(1, N);
g_syn_nmda(g_syn_nmda <= 0) = init_min_syn_frac .* g_syn_nmda_mean;   % don't let any synapses be less than init_min_g_syn

g_syn_gaba_mean = syn_scaling .* PSPgaba;
std_gaba = std_syn_frac .* g_syn_gaba_mean;
g_syn_gaba = g_syn_gaba_mean .* ones(1, N) + std_gaba .* randn(1, N);
g_syn_gaba(g_syn_gaba <= 0) = init_min_syn_frac .* g_syn_gaba_mean;   % don't let any synapses be less than init_min_g_syn

%%%%% STDP

msn6 = dt ./ tau_neg;
msn7 = dt ./ tau_pos;

%%% eligibility traces

msn10 = dt ./ tau_elig_pos;
msn11 = dt ./ tau_elig_neg;

%%%% DA
% d[DA]/dt = -[DA] + m_0

da_time = round(time_DA ./ dt);     % time for phasic dopamine in time steps

msn5 = dt ./ tau_phi;               % optimisation

%%%%  LR = 20;                % "learn rate" parameter in the old code with single DA pulse

%%%% DA old formulation

% phi_0 = LR ./ tau_phi;  % to give same effective (integrated) effect
%
% m_0 =  phasic_tonic_ratio .* phi_0 .* (tau_phi ./ tau_elig_ltp);

%%%% DA new (simplified formulation)

msn8 = da_tonic .* (dt ./ tau_phi);

% DA_conc now defined above when defining factors for DA modulation

%%%% pattern matching
strong_aff_perm = randp;
strong_aff_inds = find(strong_aff_perm <= No_strong_affs);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% end parameters
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% initialise
g_syn_ss = []; % stores  over trials

diag_pos_eligs_ss = [];
diag_neg_eligs_ss = [];
diag_pos_plas_ss = [];
diag_neg_plas_ss = [];

V = vr;

u = 0;
S_neg = 0;
S_pos = 0;
pos_eligs = zeros(N,1);
neg_eligs = zeros(N,1);
total_trial_count = 1;

%% start

%% slow wave pars
period_slow_wave = 1 ./ freq_slow_wave;
N_half_period_sw = round(0.5 .* period_slow_wave ./ dt);

Ttrial_time = sum(periods);
if Ttrial_time <  0.5 .* period_slow_wave
    error('STDP_stim:trialtime', 'Trial time is less than half a slow wave')
end

sw_stim_reset_phase = round(half_period_frac .* (0.5 ./ freq_slow_wave) ./ dt);
fractional_phase_slow_wave = 0;     % integral number of increments of dt
% into a half cycle
slow_wave_half_cycle_type = 1;      % 'up' state is 1, down state is 0

p_up = (1 + A_depth) .* mean_sw_rate .* dt;
p_down = (1 - A_depth) .* mean_sw_rate .* dt;

sw_alpha = A_depth .* sw_alpha_star;

%%%%%%%%%%%%%%%%

T_spikes_ss = []; % all spikes in a run - only use if diagnsotic on

Gampa = 0;
Gnmda = 0;
Ggaba = 0;

for xpt_phase = 1:No_phases
    No_trials = trial_counts(xpt_phase);
    trial_type = xpt_phases(1, xpt_phase);
    switch trial_type
        case NO_STIM_NO_PHASIC_DA
            phi = 0;
        case WITH_STIM_NO_PHASIC_DA
            phi = 0;
        case WITH_STIM_WITH_PHASIC_DA
            phi = da_phasic;
        case NO_STIM_WITH_PHASIC_DA
            phi = da_phasic;
        case WITH_STIM_WITH_REDUCED_PHASIC_DA
            phi = reduced_phasic_factor .* da_phasic;
        case PATTERN_MATCH_WITH_PHASIC_DA
            phi = da_phasic;
        case PATTERN_MATCH_WITH_DA_DIP
            phi = 0;
        case RANDOM_PATTERNS
            phi = 0;
        case PATTERN_DISCOVERY
            phi = da_phasic;
        case PATTERN_MATCH_WITH_SALIENCE_DECAY
            phi = da_phasic;
        otherwise
            error('STDP_stim:XptPhasesInvalid', 'Invalid expt phase');
    end
    
    
    for trial = 1:No_trials
        
        % make spikes
        
        Tspikes = [];
        
        switch trial_type
            case NO_STIM_NO_PHASIC_DA       % cortical slow wave only
                
                tspan = 0:dt:Ttrial_time;
                lt = length(tspan);
                
                [Tspikes, fractional_phase_slow_wave, slow_wave_half_cycle_type] = ...
                    make_slow_wave_stim(lt, N_half_period_sw, ...
                    fractional_phase_slow_wave, slow_wave_half_cycle_type, ...
                    N, dt, sw_alpha, p_up, p_down);
                
            case WITH_STIM_NO_PHASIC_DA     % cortical slow wave and stim
                
                tspan = 0:dt:half_period;
                lt = length(tspan);
                
                [spikes, fractional_phase_slow_wave, slow_wave_half_cycle_type] = ...
                    make_slow_wave_stim(lt, N_half_period_sw, ...
                    fractional_phase_slow_wave, slow_wave_half_cycle_type, ...
                    N, dt, sw_alpha, p_up, p_down);
                
                Tspikes = [Tspikes spikes];
                
                %%%%%%%%%%%% stim %%%%%%%%%%%%%
                tspan = 0:dt:stim_duration;
                seg = length(tspan);
                p = stim_rate .* dt;
                
                T = make_spikes(N, seg, p, correlation_in_stim);
                Tspikes = [Tspikes T];
                
                %%%%%%%%%%% final period after reset %%%%%%%%%%%%
                slow_wave_half_cycle_type = sw_stim_reset_period_type;
                fractional_phase_slow_wave = sw_stim_reset_phase;
                
                tspan = 0:dt:half_period;
                lt = length(tspan);
                
                [spikes, fractional_phase_slow_wave, slow_wave_half_cycle_type] = ...
                    make_slow_wave_stim(lt, N_half_period_sw, ...
                    fractional_phase_slow_wave, slow_wave_half_cycle_type, ...
                    N, dt, sw_alpha, p_up, p_down);
                
                Tspikes = [Tspikes spikes];
                
            case {WITH_STIM_WITH_PHASIC_DA, WITH_STIM_WITH_REDUCED_PHASIC_DA}
                % no slow wave (activated state) and stim
                
                for i=1:length(periods)
                    
                    tspan = 0:dt:periods(i);
                    seg = length(tspan);
                    alpha = alphas_stim(i);
                    rate = rates_stim(i);
                    p = rate .* dt;
                    
                    T = make_spikes(N, seg, p, alpha);
                    Tspikes = [Tspikes T];
                    
                end
                
            case NO_STIM_WITH_PHASIC_DA     % no slow wave (activated state) no stim
                for i=1:length(periods)
                    
                    tspan = 0:dt:periods(i);
                    seg = length(tspan);
                    alpha = alphas_no_stim(i);
                    rate = rates_no_stim(i);
                    p = rate .* dt;
                    
                    T = make_spikes(N, seg, p, alpha);
                    Tspikes = [Tspikes T];
                    
                end
            case {PATTERN_MATCH_WITH_PHASIC_DA, PATTERN_MATCH_WITH_DA_DIP}
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
                %%% burst of salience
                tspan = 0:dt:salience_time;
                lt = length(tspan);
                p_hi = rate_hi .* dt;
                p_lo = rate_lo .* dt;
                strong_af_T = make_spikes(No_strong_affs, lt, p_hi, ...
                    alpha_hi);
                weak_af_T = make_spikes(N - No_strong_affs, lt, p_lo, ...
                    alpha_lo);
                all_T = [strong_af_T; weak_af_T];
                T = all_T(strong_aff_perm, :);
                Tspikes = [Tspikes T];
                
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
            case RANDOM_PATTERNS
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
                %%% burst of salience
                tspan = 0:dt:salience_time;
                lt = length(tspan);
                p_hi = rate_hi .* dt;
                p_lo = rate_lo .* dt;
                strong_af_T = make_spikes(No_strong_affs, lt, p_hi, ...
                    alpha_hi);
                weak_af_T = make_spikes(N - No_strong_affs, lt, p_lo, ...
                    alpha_lo);
                all_T = [strong_af_T; weak_af_T];
                aff_perm = randperm(N);
                T = all_T(aff_perm, :);
                Tspikes = [Tspikes T];
                
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
            case PATTERN_DISCOVERY
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
                %%% burst of salience
                tspan = 0:dt:salience_time;
                lt = length(tspan);
                p_hi = rate_hi .* dt;
                p_lo = rate_lo .* dt;
                
                % randomisation of strong afferents
                 
                % The following code keeps a fraction of the original
                % permuation fixed and randonly re-permutes the rest
                % Thus when frac_olap = 1, the originalk perm is kept and
                % when frac_olpa = 0, none of themn are preserved.
                % see test_overlap_st_af_code.m for a demo
                
                frac_olap = trial ./ No_trials; % No_trials is for this phase only
                rds = rand(1,N);                % N is total number of afferents
                id_olap = rds < frac_olap;
                s_fixed = strong_aff_perm .* id_olap; % strong_aff_perm is rand perm of all N synpses
                s_new_perm = s_fixed;
               
                non_fix = setdiff(strong_aff_perm, s_fixed);
                newp = randperm(length(non_fix));
                non_fix = non_fix(newp);
                id_nolap = ~id_olap;
                s_new_perm(find(id_nolap)) = non_fix;
                
                %%%%%%%%%%%%%%%%%%%
                
                strong_af_T = make_spikes(No_strong_affs, lt, p_hi, ...
                    alpha_hi);
                weak_af_T = make_spikes(N - No_strong_affs, lt, p_lo, ...
                    alpha_lo);
                all_T = [strong_af_T; weak_af_T];
                T = all_T(s_new_perm, :);
                Tspikes = [Tspikes T];
                
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
            case PATTERN_MATCH_WITH_SALIENCE_DECAY
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
                %%% burst of salience
                tspan = 0:dt:salience_time;
                lt = length(tspan);
                p_hi = rate_hi .* dt;
                p_lo = rate_lo .* dt;
                
                % randomisation of strong afferents
                
                % The following code keeps a fraction of the original
                % permuation fixed and randonly re-permutes the rest
                % Thus when frac_olap = 1, the originalk perm is kept and
                % when frac_olpa = 0, none of themn are preserved.
                % see test_overlap_st_af_code.m for a demo
                
                frac_olap = 1 - trial ./ No_trials; % No_trials is for this phase only
                rds = rand(1,N);                % N is total number of afferents
                id_olap = rds < frac_olap;
                s_fixed = strong_aff_perm .* id_olap; % strong_aff_perm is rand perm of all N synpses
                s_new_perm = s_fixed;
               
                non_fix = setdiff(strong_aff_perm, s_fixed);
                newp = randperm(length(non_fix));
                non_fix = non_fix(newp);
                id_nolap = ~id_olap;
                s_new_perm(find(id_nolap)) = non_fix;
                
                %%%%%%%%%%%%%%%%%%%
                
                strong_af_T = make_spikes(No_strong_affs, lt, p_hi, ...
                    alpha_hi);
                weak_af_T = make_spikes(N - No_strong_affs, lt, p_lo, ...
                    alpha_lo);
                all_T = [strong_af_T; weak_af_T];
                T = all_T(s_new_perm, :);
                Tspikes = [Tspikes T];
                
                %% background activity
                tspan = 0:dt:half_period;
                lt = length(tspan);
                alpha = alpha_lo;
                p = rate_lo .* dt;
                T = make_spikes(N, lt, p, alpha);
                Tspikes = [Tspikes T];
                
            otherwise
                error('STDP_multi_syn_DA_dyn_batch:XptPhasesInvalid', 'Invalid expt phase');
        end % generating spikes for the trial
        
        if spike_diag
            T_spikes_ss = [T_spikes_ss Tspikes];
        end
        
        % make gaba spikes
        gaba_spikes = [];
        tspan = 0:dt:half_period;
        lt = length(tspan);
        p = rate_lo .* dt;
        T = make_spikes(N, lt, p, 0);
        gaba_spikes = [gaba_spikes T];
        
        %%% burst of salience
        tspan = 0:dt:salience_time;
        lt = length(tspan);
        p_hi = rate_hi .* dt;
        p_lo = rate_lo .* dt;
        strong_af_T = make_spikes(No_strong_affs, lt, p_hi, 0);
        weak_af_T = make_spikes(N - No_strong_affs, lt, p_lo, 0);
        all_T = [strong_af_T; weak_af_T];
        aff_perm = randperm(N);
        T = all_T(aff_perm, :);
        gaba_spikes = [gaba_spikes T];
        
        
        %% background activity
        tspan = 0:dt:half_period;
        lt = length(tspan);
        p = rate_lo .* dt;
        T = make_spikes(N, lt, p, 0);
        gaba_spikes = [gaba_spikes T];
        
        %%%%%%%%%%%%% do MSN
        
        post_spikes_times = [];
        
        if  diag_elig_synapse
            diag_pos_eligs = zeros(1,trial_time);
            diag_neg_eligs = zeros(1,trial_time);
            diag_pos_plas = zeros(1,trial_time);
            diag_pos_plas = zeros(1,trial_time);
        end
        
        % make trial-based noise for DA conc
        DA_noise = DA_std .* randn;
        
        trial_time = size(Tspikes,2);
        for i = 1:trial_time
            time = dt .* i;
            
            %% set up any possible pos-timing plasticity
            S_pos = S_pos + Tspikes(:,i);  % set up potential pos-timing plasticity
            % based on pre-synaptic activity
            % waiting for a post-synaptic
            % spike. Implicit multiplicative
            % constant of 1, and scaling
            % occurs elsewhere.
            
            S_pos = S_pos - S_pos .* msn7; % decay the pos-timing potential
            
            %%% new neg-timing eligibilities
            neg_eligs = neg_eligs + S_neg .* Tspikes(:,i); % deploy the potential
            % negative timing
            % plasticity
            neg_eligs = neg_eligs - neg_eligs .* msn11; % decay the
            % neg-timing eligibilities
            
            % find synaptic current using first order ODE method
            
            % In the following the variables PSPX (X is ampa etc) are
            % conductances in spite of their name!
            
            Gampa = Gampa + (g_syn_ampa * Tspikes(:,i)); %normalisationwith ts_ampa done in init.
            Gampa = Gampa * SynExp_ampa;
            
            Gnmda = Gnmda + (g_syn_nmda * Tspikes(:,i));
            Gnmda = Gnmda * SynExp_nmda;
            
            Ggaba = Ggaba + (g_syn_gaba * gaba_spikes(:,i));
            Ggaba = Ggaba * SynExp_gaba;
            
            B_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-V .* 0.062));
            
            I_nmda = B_nmda .* Gnmda .* (Enmda - V);
            if strcmp(neuron_type, 'D1')
                I_nmda = D1_syn_fact .* I_nmda;
            end
            
            I_ampa = Gampa .* (Eampa - V);
            if strcmp(neuron_type, 'D2')
                I_ampa = D2_syn_fact .* I_ampa;
            end
            
            I_gaba = Ggaba .* (Egaba - V);
            
            I_syn = I_gaba + I_nmda + I_ampa;
            
            % the main membrane dynamics
            V = V + dt_over_C .* (k .* (V - vr) .* (V - vt) -u + I_syn);
            
            u = u + ms_dt .* a .* (b .* (V - vr) - u);
            
            % spikes?
            if V >= vpeak
                
                V = vpeak;
                V = c;
                u = u + d;
                
                post_spikes_times = [post_spikes_times time];
                
                %%% set up ltd
                S_neg = S_neg + 1;              % set up potential neg-timing plasticity
                % awaiting pre-synaptic
                % spikes. Use constant of 1
                % as scaling will be
                % absorbed in other
                % parameters (k_hats and LR)
                %% new pos-timing  eligibilities
                pos_eligs = pos_eligs + S_pos;  % deploy the potential
                % pos-timing eligibility
            end;
            
            S_neg = S_neg - S_neg .* msn6;      % decay the neg-timing potential
            pos_eligs = pos_eligs - pos_eligs .* msn10; % decay the
            % pos-timing eligibilities
            
            %% update the DA
            % This conditional code could be in the next outer loop (over
            % trials) but leaving for now...
            if trial_type == PATTERN_MATCH_WITH_PHASIC_DA
                phi = da_phasic .* exp(-trial ./ tau_da_habituate);
            end
            
            % In the next possibility, DA_conc has to be interpreted as
            % a 'virtual' concentration since we allow it go
            % negative. This s kludge to make the DA dip behave like zero
            % DA for an extended period. Things are fixed in the line
            % commneted a %% ***** %%% below, where all negative DA concs get
            % mapped into the same value of alpha
            
            % see remark above about placement of this code too
            if trial_type == PATTERN_MATCH_WITH_DA_DIP
                phi = -da_phasic  .* exp(-trial ./ tau_da_habituate);
            end
            
            % pattern discovery
            if trial_type == PATTERN_DISCOVERY
                do_phasic = rand < (trial ./ No_trials);
                if do_phasic
                    phi = da_phasic;
                else
                    phi = 0;
                end
            end
            if trial_type == PATTERN_MATCH_WITH_SALIENCE_DECAY
                do_phasic = rand < (1 - trial ./ No_trials);
                if do_phasic
                    phi = da_phasic;
                else
                    phi = 0;
                end
            end
            if i == da_time
                phi = phi + DA_noise;  % add noise to DA conc
                DA_conc = DA_conc + phi;
            end
            DA_conc = DA_conc - DA_conc .* msn5 + msn8;
            
            % find alpha_da (blend of two plasticity regimes)
            
            if DA_fn_type == LINEAR
                if DA_conc <= DA_lo %%% ***** %%%%
                    alpha_da = 0;
                elseif (DA_conc > DA_lo) && (DA_conc < DA_hi)
                    alpha_da = (DA_conc - DA_lo) ./ (DA_hi - DA_lo);
                else
                    alpha_da = 1;
                end
            elseif DA_fn_type == NAKA_RUSHTON
                if DA_conc <= DA_lo %%% ***** %%%%
                    alpha_da = 0;
                else
                    NK1 = (DA_conc - DA_lo) .^ NK_rho;
                    alpha_da = NK_max .* NK1 ./ (NK1 + NK_theta .^ NK_rho);
                end
            else
                error('STDE_Shen:Invalidalphafn', 'Invalid alpha definition function');
            end
            alpha_da(alpha_da > 1) = 1;
            
            % update the phi factors in the full DA model
            if DA_conc <= 0
                D1 = 0;
            else
                phi_nk1 = DA_conc .^ phi_rho;
                D1 =  phi_max .* phi_nk1 ./ (phi_nk1 + phi_theta .^ phi_rho);
            end
            D2 = D1;
            
            if strcmp(neuron_type, 'D1')
                vr = vr0 .* (1 + D1 .* KIR);
                d = d0 .* (1 - D1 .* LCA);
                D1_syn_fact = (1 + cD1 .* D1);
            elseif strcmp(neuron_type, 'D2')
                k = k0 .* (1 - alpha .* D2);
                D2_syn_fact = (1 - cD2 .* D2);
            else
                error('STDP_Shen_batch:neuron_type', 'Invalid neurontype (D1 or D2)')
            end
            
            % The following assumes that he sign of plastic change is
            % built into the k_hat values
            pos_plasticity = alpha_da .* k_hat_pos_hi + ...
                (1 - alpha_da) .* k_hat_pos_lo;
            neg_plasticity = alpha_da .* k_hat_neg_hi + ...
                (1 - alpha_da) .* k_hat_neg_lo;
            
            elig = pos_eligs(live_syn_ampa) .* pos_plasticity + neg_eligs(live_syn_ampa) .* neg_plasticity;
            
            g_syn_ampa(live_syn_ampa) = g_syn_ampa(live_syn_ampa) + LR .* elig' .* dt;
            zs = find(g_syn_ampa < 0);
            g_syn_ampa(zs) = 0;
            % live_syn_ampa(zs) = false;
            g_syn_ampa(g_syn_ampa > max_syn_ampa) = max_syn_ampa;
            
            if diag_elig_synapse
                diag_pos_eligs(i) = pos_eligs(diag_elig_synapse);
                diag_neg_eligs(i) = neg_eligs(diag_elig_synapse);
                diag_pos_plas(i) = pos_plasticity;
                diag_neg_plas(i) = neg_plasticity;
            end
        end
        if diag_g_syn
            g_syn_ss = [g_syn_ss; g_syn_ampa];
        end
        if diag_elig_synapse
            diag_pos_eligs_ss = [diag_pos_eligs_ss diag_pos_eligs];
            diag_neg_eligs_ss = [diag_neg_eligs_ss diag_neg_eligs];
            diag_pos_plas_ss = [diag_pos_plas_ss diag_pos_plas];
            diag_neg_plas_ss = [diag_neg_plas_ss diag_neg_plas];
        end
        
        post_spikes_ss{total_trial_count} = post_spikes_times;
        total_trial_count = total_trial_count + 1;
    end
end

basename = 'results';
fname = [basename num2str(xpt_no)];

eval(['save ' fname ' all_pars strong_aff_inds  post_spikes_ss g_syn_ss g_syn_nmda g_syn_gaba T_spikes_ss diag_pos_eligs_ss diag_neg_eligs_ss diag_pos_plas_ss diag_neg_plas_ss']);

return

