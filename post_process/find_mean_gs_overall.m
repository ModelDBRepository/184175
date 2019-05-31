function [means_strong, means_weak, iprods] = find_mean_gs_overall(xpt_nos)

% find mean (over xpts) weak and strong synaptic conductances g_syn, over all trials
%
% XPT_NOS is range of  experiment number to retrieve results file RESULTS#
% mean is taken over these
% 
% MEAN_STRONG is a vector of  mean of strong synspses ofve experiments with
% components for each trial
%
% MEAN_WEAK is similar for weak synpses
%
% IPRODS is inner product of synapse vector (all synapses) with itself, for each trial

%% assume same number synapse AND same striong afferent set
fname = ['results' num2str(xpt_nos(1))];
load(fname, 'g_syn_ss', 'strong_aff_inds');
N_synapses = size(g_syn_ss, 2);
weak_aff_indices = setdiff([1:N_synapses], strong_aff_inds);

No_xpts = length(xpt_nos);
No_trials = size(g_syn_ss, 1);
gs_all = zeros(No_xpts, No_trials, N_synapses);

for i = 1:No_xpts
    fname = ['results' num2str(xpt_nos(i))];
    load(fname, 'g_syn_ss', 'strong_aff_inds');
    
    gs_strong = g_syn_ss(:,strong_aff_inds); % trials x synapses
    mean_gs_strong(i,:) = mean(gs_strong,2); % 1 x trials
    
    gs_non_strong = g_syn_ss(:,weak_aff_indices);
    mean_gs_nonstrong(i,:) = mean(gs_non_strong, 2); % 1 x trials
    
    gs_all(i,:,:) = g_syn_ss;
   
end

means_strong = mean(mean_gs_strong);
means_weak = mean(mean_gs_nonstrong);

mean_gs_all = mean(gs_all, 1);
sumsq = mean_gs_all .* mean_gs_all;
iprods = squeeze(sum(sumsq, 3));
iprods = sqrt(iprods);
iprods = iprods';

figure(1)
plot(means_strong);
hold on
plot(means_weak, 'r');
hold off

set(gcf, 'PaperOri', 'portrait')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPos', [0 0 20 14])
% (20 and 14 are thus measured in cm)
fnme = ['mean_gs_overall.png'];
print(gcf, '-dpng', fnme, '-r100')

figure(2)
plot(iprods, 'k');
set(gcf, 'PaperOri', 'portrait')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPos', [0 0 20 14])
% (20 and 14 are thus measured in cm)
fnme = ['mean_iprods_overall.png'];
print(gcf, '-dpng', fnme, '-r100')

end



