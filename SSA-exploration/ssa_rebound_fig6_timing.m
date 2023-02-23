%% timing experiment results
% distributions showing rebound strength, shading significant units (for
%   each song)
% summary data showing rebound strength doesn't correlate with surprise
%   (for each song)
% sorted imagescs showing that rebounds do get more common as timescale
%   increases

%% prep
clear
close all
pipeline_out = 'D:\timing_expmts\multibird_summary_timing.mat';

%% load multibird summary
load(pipeline_out)

% motif and bout lengths

%% rbd strength per unit
% for each unit, for each song, rebound is significant if it is significantly greater than 0
% ie, if AVG - 1.65(STD) of rebound is greater than 0, rbd is significant
%
tmr = total_motif_resp;
tmr_norm = tmr ./ max(tmr, [], 1);

% habtest1 = 125.5;
% testhab2 = 870.5;
% habtest2 = 995.5;
% testhab3 = 1735.5;
% habtest3 = 1860.5;


key_motifs = [1 125;
    871 995;
    1736 1860];

rbds_by_song = zeros(sum(num_goods),3);
sig_by_song = zeros(sum(num_goods),3);

rbd_frac = zeros(sum(num_goods),3);
sig_frac = zeros(sum(num_goods),3);

for song_idx = 1:3
    this_m1_ids = key_motifs(song_idx,1):5:(key_motifs(song_idx,2));
    this_m5_ids = key_motifs(song_idx,1)+4:5:(key_motifs(song_idx,2));
    
    this_rbds = tmr(this_m1_ids(2:end),:) - ...
        tmr(this_m5_ids(1:end-1),:); % if rebound exists, this number will be positive
    
    avg_rbds = mean(this_rbds,1); % for each unit
    std_rbds = std(this_rbds,0,1); % for each unit
    sig_rbds = find((avg_rbds - 1.65*(std_rbds)) > 0);
       
    
    this_rbds_norm = tmr_norm(this_m1_ids(2:end),:) - ...
        tmr_norm(this_m5_ids(1:end-1),:); % if rebound exists, this number will be positive
    
    avg_rbds_frac = mean(this_rbds_norm,1); % for each unit
    std_rbds = std(this_rbds_norm,0,1); % for each unit
    sig_rbds_norm = find((avg_rbds_frac - 1.65*(std_rbds)) > 0);
    
    
    rbds_by_song(:, song_idx) = avg_rbds;
    sig_by_song(sig_rbds, song_idx) = 1;
    
    rbd_frac(:, song_idx) = avg_rbds_frac;
    sig_frac(sig_rbds_norm, song_idx) = 1;
end
   
%% plot and shade significant rebounds
% bin_edges = -0.5:0.025:0.5;
bin_edges = -18:3:30; % happen to know that these are bounds

for song_idx = 1:3    
    this_sig = find(sig_by_song(:,song_idx));
    
    figure(8181)
    subplot(1,3,song_idx)
    histogram(rbds_by_song(:,song_idx),'BinEdges', bin_edges)
    hold on;
    histogram(rbds_by_song(this_sig,song_idx),'BinEdges', bin_edges)
    hold off;
    xlabel('Rebound (Hz)')
    ylabel('Number of neurons')
    
    title(['Rebounds, Song ' num2str(song_idx)])
end

h5 = figure(8181);
print(['figure_pieces/fig6_timing_rbdDist_raw'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_timing_rbdDist_raw.fig'])

%% plot and shade significant rebounds
% bin_edges = -0.5:0.025:0.5;
bin_edges = -0.5:0.05:0.5; % happen to know that these are bounds

for song_idx = 1:3    
    this_sig = find(sig_frac(:,song_idx));
    
    figure(8182)
    subplot(1,3,song_idx)
    histogram(rbd_frac(:,song_idx),'BinEdges', bin_edges)
    hold on;
    histogram(rbd_frac(this_sig,song_idx),'BinEdges', bin_edges)
    hold off;
    xlabel('Rebound (prop. of max.)')
    ylabel('Number of neurons')
    
    title(['Rebounds, Song ' num2str(song_idx)])
end

h5 = figure(8182);
print(['figure_pieces/fig6_timing_rbdDist_frac'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_timing_rbdDist_frac.fig'])

%% summary showing rebounds aren't function of surprise
% for each song (just 2 and 3), for each "sig rbd unit", plot normalized
% rebound (maybe plot each unit)

% need IBIs - can get without going back bc we have motifInfo

% habtest1 = 125.5;
% testhab2 = 870.5;
% habtest2 = 995.5;
% testhab3 = 1735.5;
% habtest3 = 1860.5;        

key_motifs = [126 870;
    996 1735;
    1861 2595];

% exclude standard IBIs to clean up plot (they're very common)
std_ibis = [0.098 0.102; % song 1 doesn't matter, no sig rbds
    0.995 1.005;
    9.99 10.01];
        
for song_idx = 1:3
    m1s = key_motifs(song_idx,1):5:key_motifs(song_idx,2);
    m5s = key_motifs(song_idx,1)+4:5:key_motifs(song_idx,2);
    bout_starts = cell2mat(motifInfo(m1s,1));
    bout_ends = cell2mat(motifInfo(m5s,3));
    
    these_IBIs = bout_starts(2:end) - bout_ends(1:end-1);
    these_IBIs = these_IBIs - 0.03; % final correction due to IMI
    
    this_rbds = tmr_norm(m1s(2:end),:) - ...
        tmr_norm(m5s(1:end-1),:); % if rebound exists, this number will be positive
    
    if song_idx == 1 % no sig rebounds during song 1, so plot all
        sig_rbds = this_rbds;
        this_sig = length(this_rbds);
    else
        this_sig = find(sig_by_song(:,song_idx));
        sig_rbds = this_rbds(:,this_sig);
    end
    
    this_y = mean(sig_rbds,2); % average across units
    this_err = std(sig_rbds, [], 2) ./ sqrt(length(this_sig)); % avging across neurons
    
    
    this_std_low = std_ibis(song_idx,1);
    this_std_high = std_ibis(song_idx,2);
    
    
    [std_IBIs] = find(...
        (these_IBIs > this_std_low) & (these_IBIs < this_std_high));
     dev_IBIs = setdiff(1:length(these_IBIs),std_IBIs);
    
    true_std = mean(these_IBIs(std_IBIs));
    std_y = mean(this_y(std_IBIs));
    std_err = std(this_y(std_IBIs)) ./ sqrt(length(std_IBIs));
     
    figure(12223)
    subplot(3,1,song_idx)
    errorbar(these_IBIs(dev_IBIs), this_y(dev_IBIs), 2*this_err(dev_IBIs), 'o')
    hold on;
    errorbar(true_std, std_y, 2*std_err, 'ro')
    hold off;
    title(['Dev Rebounds by IBI, Song ' num2str(song_idx)])
    xlabel('Inter-bout Interval')
    ylabel('Rebound (prop. of max.)')
    set(gca, 'FontSize', 12)
end

h5 = figure(12223);
h5.Position = [100 100 600 800];
print(['figure_pieces/fig6_timing_summary_frac'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_timing_summary_frac.fig'])

%% imagescs
% just plot all motif response data sorted by rebound strength

key_motifs = [1 125;
    871 995;
    1736 1860];


tmr_norm = total_motif_resp ./ max(tmr, [], 1);

for song_idx = 1:3
    this_rbds = rbds_by_song(:,song_idx);
    
    [~,rbd_sort] = sort(this_rbds,'descend');
    
    song_start = key_motifs(song_idx,1);
    song_end = key_motifs(song_idx,2);
    
    figure(8383)
    subplot(3,1,song_idx)
    imagesc(tmr_norm(song_start:song_end,rbd_sort)')
    title(['Motif Responses, Song ' num2str(song_idx)])
    xlabel('Motif Number')
    ylabel('Unit ID')
end
    
h5 = figure(8383);
print(['figure_pieces/fig6__timing_rbd_imagesc'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_timing_rbd_imagesc.fig'])