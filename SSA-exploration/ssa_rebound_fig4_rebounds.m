% example raster of rebounding neuron - call out both types (across bouts
%   and across phases)
% raster of same unit, showing what's going on
% trajectory of same unit, showing that cross-phase rebound is basically 
%   just a return to earlier stage in habituation
% imagescs showing (both kinds of) rebound across datasets
% distributions of rebounds (both kinds, all 4 songs)

%% prep
clear
close all
pipeline_out = 'D:\ssa_expmts\multibird_summary_ssa.mat';

% bird 2, su 006, song1
% response curve looks great, raster looks great (but key motifs should be
% first 10 of p1 and first 5 of p2), response curve is pretty good too
example_bird = [2];
example_unit = [6];
example_song = [1];

%% load multibird summary
load(pipeline_out)

%% response summary example - can we focus on song 1?

psth_bird_id = example_bird;
psth_unit_type = 'su';
psth_unit_id = example_unit;
es = example_song;

% transition times for the stim
s1s2 = 125;
s2s3 = 250;
s3s4 = 375;
s4s1b = 500;
s1bs2b = 550;
s2bs3b = 600;
s3bs4b = 650;
std_end = 700;

km = [0 125 250 375 500 550 600 650 700];
% s1 s2 s3 s4
if strcmp(psth_unit_type, 'su')
    if psth_bird_id > 1
        base = sum(num_goods(1:psth_bird_id-1));
    else
        base = 0;
    end
    actual_id = base + psth_unit_id;
    this_row = total_motif_resp(1:std_end,actual_id);
else
    this_row = total_motif_resp_mu(1:std_end,actual_id);
end

h2=figure(1006);
subplot(3,1,[1 2])
plot(km(es)+1:km(es+1), this_row(km(es)+1:km(es+1)), 'r-')
hold on;
plot(126:175, this_row(km(es+4)+1:km(es+5)), 'r.')
hold off;
ylabel('resp strength (hz)')
title(['Std Motif Responses'])
xlim([1 175])

subplot(313)
plot(this_row(km(es)+1:km(es)+25), 'r*-')
xlim([1 25])
title('First Responses to Std 1')
ylim([min(0, min(this_row)) max(this_row)+10])

print('figure_pieces/fig4_ssaExample', '-dsvg', '-r300')
saveas(h2, 'figure_pieces/fig4_ssaExample.fig')

%% raster example - looks okay, but maybe exclude
for ex_idx = 1:length(example_bird)
    ex_bird = example_bird(ex_idx);
    ex_unit = example_unit(ex_idx);
    ex_song = example_song(ex_idx);
    
    switch ex_bird
        case 1
            load('D:\ssa_expmts\pipelines\pipeline_9059_210811_LH_NCM_g0.mat')
        case 2
            load('D:\ssa_expmts\pipelines\pipeline_9072_210810_LH_NCM_g0_v2.mat')
        case 3
            load('D:\ssa_expmts\pipelines\pipeline_9074_210811_LH_NCM_g0.mat')
    end
    
    raster_win = [-0.025 0.2];
    
    impt_trials = [1 501;...
        126 551;...
        251 601;...
        376 561];
    
    songIdx = ex_song;
    this_mot_len = (motifInfo{impt_trials(songIdx,1),3} - motifInfo{impt_trials(songIdx,1),1}) + 0.2;
    this_rast_win = raster_win + [0 this_mot_len];
    ats_1 = motifInfo(impt_trials(songIdx,1):impt_trials(songIdx,1)+9);
    ats_2 = motifInfo(impt_trials(songIdx,2):impt_trials(songIdx,2)+4);
    
    
    
    these_ats = cell2mat([ats_1'; ats_2']);
    these_gu = sp.cids(find(sp.cgs==uType));
    
    these_spikes = sp.st(sp.clu == these_gu(ex_unit)); 
    
    ssa_rebound_rasters_v2(these_spikes, these_ats, ex_unit, this_rast_win, micData, dataset_name, 'su')
    
    h5 = gcf;
    print(['figure_pieces/fig4_rbdRasterEx' '_' num2str(ex_idx)], '-dsvg', '-r300')
    saveas(h5, ['figure_pieces/fig4_rbdRasterEx'  '_' num2str(ex_idx) '.fig'])
end

%% trajectory example
% first 5 responses p1, last 5 responses p1, first 5 responses p2
for ex_idx = 1:length(example_bird)
    ex_bird = example_bird(ex_idx);
    ex_unit = example_unit(ex_idx) -3;
    ex_song = example_song(ex_idx);
    
    switch ex_bird
        case 1
            these_units = 1:num_goods(1);
        case 2
            these_units = num_goods(1)+1:num_goods(1)+num_goods(2);
        case 3
            these_units = num_goods(2)+1:sum(num_goods);
    end
    this_unit = these_units(ex_unit);
    
    relevant_motif_chunks = [1 5; 121 125; 126 130];
    
    line_colors = {'r-', 'k-', 'r--'};
    
    eval(['this_resp_vector = big_rv' num2str(ex_song) ';'])
    this_resp_piece = squeeze(this_resp_vector(this_unit,:,:));
%     this_resp_piece = this_resp_piece ./ max(this_resp_piece); % causes issues...
    
    for lineIdx = 1:3
        impt_motifs = relevant_motif_chunks(lineIdx,:);
        good_motifs = this_resp_piece(impt_motifs(1):impt_motifs(2),:);
        this_line = mean(good_motifs,1);
        figure(101212)
        plot(this_line, line_colors{lineIdx})
        if lineIdx == 1
            hold on;
        elseif lineIdx == 3
            hold off;
        end
    end
    
    title(['Bird ' num2str(ex_bird) ', Unit ' num2str(ex_unit) ', Song' num2str(songIdx)])
    ylabel('Resp. Str.')
    xlabel('Time')
    
    h5 = gcf;
    print(['figure_pieces/fig4_rbdCurveEx' '_' num2str(ex_idx)], '-dsvg', '-r300')
    saveas(h5, ['figure_pieces/fig4_rbdCurveEx'  '_' num2str(ex_idx) '.fig'])
end

%% distributions showing rebounds across experiment
% (we'll sort imagesc by the intra-phase rbd later)

tmr_norm = total_motif_resp ./ max(total_motif_resp, [], 1);

key_motifs = [1 126 251 376; 
    501 551 601 651];

for song_idx = 1:4
    this_m1_ids = key_motifs(1,song_idx):5:key_motifs(1,song_idx)+124;
    this_m5_ids = key_motifs(1,song_idx)+4:5:key_motifs(1,song_idx)+124;
    
    this_p1b25 = key_motifs(1,song_idx)+120:key_motifs(1,song_idx)+124;
    this_p2b1 = key_motifs(2,song_idx):key_motifs(2,song_idx)+4;
    
    
    this_p1_rbds = tmr_norm(this_m1_ids(2:end),:) - ...
        tmr_norm(this_m5_ids(1:end-1),:); % if rebound exists, this number will be positive
    avg_p1_rbds = mean(this_p1_rbds,1); % for each unit
    
    this_p2_rbds = mean(tmr_norm(this_p2b1,:)) - ...
        mean(tmr_norm(this_p1b25,:));
    % no averaging needed because just one number per unit
    
    eval(['rbds_p1_s' num2str(song_idx) ' = avg_p1_rbds;'])
    eval(['rbds_p2_s' num2str(song_idx) ' = this_p2_rbds;'])
    
    figure(8181)
    subplot(2,2,song_idx)
    histogram(avg_p1_rbds)
    title(['Inter-bout Rebound, Song ' num2str(song_idx)])
    
    figure(8282)
    subplot(2,2,song_idx)
    histogram(this_p2_rbds)
    title(['Inter-phase Rebound, Song ' num2str(song_idx)])
end

h5 = figure(8181);
print(['figure_pieces/fig4_ib_rbds'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig4_ib_rbds.fig'])

h5 = figure(8282);
print(['figure_pieces/fig4_ip_rbds'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig4_ip_rbds.fig'])

%% imagescs showing p1 rebounds

key_motifs = [1 125;
    126 250;
    251 375;
    376 500];

for song_idx = 1:4
    this_song_response =...
        total_motif_resp(key_motifs(song_idx,1):key_motifs(song_idx,2),:);
    tsr_norm = this_song_response ./ max(this_song_response,[],1);
    
    eval(['this_rbds = rbds_p1_s' num2str(song_idx) ';'])
    
    [~,rbd_sort] = sort(this_rbds,'descend');
    
    figure(8383)
    subplot(2,2,song_idx)
    imagesc(tsr_norm(:,rbd_sort)')
    title(['P1 Motif Responses, Song ' num2str(song_idx)])
    xlabel('Motif Number')
    ylabel('Unit ID')
end
    
h5 = figure(8383);
print(['figure_pieces/fig4_rbd_imagesc'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig4_rbd_imagesc.fig'])
