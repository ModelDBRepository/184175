function [data_hist, data_scatter] = find_mean_gs(max_g, No_bins, xpt_nos, plot_trial)

% find mean synaptic conductances g_syn (over series of xpts), and histogram bins for selected trial
%
% MAX_G is maximum synaptic conductance to show in bar graphs -
% probably needs to be around 2)
%
% NO_BINS number of bins in a bar graph (try 20)
%
% XPT_NOS is range of  experiment number to retrieve results file RESULTS#
% mean is taken over this range
%
% PLOT_TRIAL - trial number to deal with
%
% DATA_HIST is an array of No_bins x 2, with columns for histograms for
% each of strong and weak synapses 
%
% DATA_SCATTER is the mean gs_at the given plot trial for each synapse
%


start_bin = 0;
end_bin = max_g;
bins = linspace(start_bin, end_bin, No_bins);

gs_sort_ss = [];
for i = 1:length(xpt_nos)
    fname = ['results' num2str(xpt_nos(i))];
    load(fname, 'g_syn_ss', 'strong_aff_inds');

    gs_select = g_syn_ss(plot_trial, :);

    N_synapses = size(g_syn_ss, 2);
    weak_aff_indices = setdiff([1:N_synapses], strong_aff_inds);

    %n(i,:) = histc(gs_select, bins);
    
    gs_strong = gs_select(strong_aff_inds);
    nstrong(i,:) = histc(gs_strong, bins);
    
    gs_non_strong = gs_select(weak_aff_indices);
    n_nonstrong(i,:) = histc(gs_non_strong, bins);
    
    gs_sort = [gs_strong gs_non_strong];
    gs_sort_ss = [gs_sort_ss; gs_sort];
   
end
gs_mean = mean(gs_sort_ss);
pcent_hs_strong = (100 ./ length(strong_aff_inds)) .* mean(nstrong);
pcent_hs_nonstrong = (100 ./ length(weak_aff_indices)) .*mean(n_nonstrong);

data_hist = [bins;pcent_hs_strong; pcent_hs_nonstrong];
data_hist = data_hist';

data_scatter = gs_mean';

figure(1)
plot(data_scatter, '*');
figure(2)
bar(data_hist(:,1), data_hist(:,2));
bar(data_hist(:,1), data_hist(:,3), 'r')




