function [mean_D1, mean_D2, new_marks] = find_mean_spikes(markers, D1_xpts, D2_xpts, trial_window)
% find moving average (windowed) mean, spike rates across all MSNs in 
% defined groups (D1 or D2)
% by 'windowed' is meant a convolution over individual trials but segmented
% across phase marker boundaries (so no 'blurring' across boundaries)
%
% MARKERS is a set of boundary trials between phases; don't include the
% last trial
% D1_xpts is a vector of xpt numbers for D1 neurons (oftne 1:10)
% D2_xpts is a vector of xpt numbers for D2 neurons (oftne 11:20)
% TRIAL_WINDOW is window size over which an average is taken 
% returns results in MEAN_PS_D1 and MEAN_PS_D2
% NEW_MARKS is a set of new phase boundaries based on teh segmentation used
%

xpt_nos = [D1_xpts D2_xpts]; 

No_xpts = length(xpt_nos);

ps_ss = [];
ps_spikes = [];
for i=1:No_xpts
    fname = ['results' num2str(xpt_nos(i))];
    load(fname, 'post_spikes_ss');
    No_trials = length(post_spikes_ss);
    for j = 1:No_trials
        spike_times = post_spikes_ss{j};
        ps_spikes(j) = length(spike_times);
    end
    
    ps_ss = [ps_ss; ps_spikes];
end
  
mean_ps_D1 =  mean(ps_ss(D1_xpts, :));
mean_ps_D2 =  mean(ps_ss(D2_xpts, :));

mask = ones(1, trial_window);
mean_D1 = [];
mean_D2 = [];
start = 1;
markers = [markers length(mean_ps_D1)]; %include final trial

for j = 1:length(markers)
    last = markers(j);
    D1_seg = mean_ps_D1(start: last);
    mu_D1_seg = conv(D1_seg, mask) ./ trial_window;
    mu_D1_seg = mu_D1_seg(trial_window: end - (trial_window - 1));
    D2_seg = mean_ps_D2(start: last);
    mu_D2_seg = conv(D2_seg, mask) ./ trial_window;
    mu_D2_seg = mu_D2_seg(trial_window: end - (trial_window - 1));
    mean_D1 = [mean_D1 mu_D1_seg];
    mean_D2 = [mean_D2 mu_D2_seg];
    marks(j) = length(mu_D1_seg);
    start = last + 1;
end

new_marks = cumsum(marks);    
figure(1)
plot(mean_D1);
hold on
my =max(mean_D1);
for k=1:length(new_marks)
    x = new_marks(k);
    plot([x x], [0 my], 'r-.');
end
set(gcf, 'PaperOri', 'portrait')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPos', [0 0 20 14])
% (20 and 14 are thus measured in cm)
fnme = ['mean_D1_spikes.png'];
print(gcf, '-dpng', fnme, '-r100')
hold off

figure(2)
plot(mean_D2);
hold on
my =max(mean_D2);
for k=1:length(new_marks)
    x = new_marks(k);
    plot([x x], [0 my], 'r-.');
end
set(gcf, 'PaperOri', 'portrait')
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPos', [0 0 20 14])
% (20 and 14 are thus measured in cm)
fnme = ['mean_D2_spikes.png'];
print(gcf, '-dpng', fnme, '-r100')
hold off

new_marks = [1 new_marks];

