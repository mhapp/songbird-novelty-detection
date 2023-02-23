% figure 3 - habituation 

%% prep and set variables

clear
close all
pipeline_out = 'D:\ssa_expmts\multibird_summary_ssa.mat';

raster_bird = [2, 3];
raster_unit = [8, 22];
raster_song = [3, 2];

% 9072 su008, song 3 (raster ok)
% 9072 su039, song 1 (raster good)
% 9059 su070, song 4 (raster just fine)
% 9074 su022, song 2 (raster good)
curve_bird = [1, 2, 2];
curve_unit = [65, 4, 22];
curve_song = [4, 2, 1];
%% example rasters
for ex_idx = 1:length(raster_bird)
    ex_bird = raster_bird(ex_idx);
    ex_unit = raster_unit(ex_idx);
    ex_song = raster_song(ex_idx);
    
    switch ex_bird
        case 1
            load('D:\ssa_expmts\pipelines\pipeline_9059_210811_LH_NCM_g0.mat')
        case 2
            load('D:\ssa_expmts\pipelines\pipeline_9072_210810_LH_NCM_g0_v2.mat')
        case 3
            load('D:\ssa_expmts\pipelines\pipeline_9074_210811_LH_NCM_g0.mat')
    end
    
    raster_win = [-0.025 0.2];
    
    impt_trials = [1 116;...
        126 241;...
        251 366;...
        376 491];
    
    songIdx = ex_song;
    this_mot_len = (motifInfo{impt_trials(songIdx,1),3} - motifInfo{impt_trials(songIdx,1),1}) + 0.2;
    this_rast_win = raster_win + [0 this_mot_len];
    ats_1 = motifInfo(impt_trials(songIdx,1):impt_trials(songIdx,1)+9);
    ats_2 = motifInfo(impt_trials(songIdx,2):impt_trials(songIdx,2)+9);
    
    
    
    these_ats = cell2mat([ats_1'; ats_2']);
    these_gu = sp.cids(find(sp.cgs==uType));
    
    these_spikes = sp.st(sp.clu == these_gu(ex_unit)); 
    
    ssa_rebound_rasters_v2(these_spikes, these_ats, ex_unit, this_rast_win, micData, dataset_name, 'su')
    
    h5 = gcf;
    print(['figure_pieces/fig3_habRasterEx' '_' num2str(ex_idx)], '-dsvg', '-r300')
    saveas(h5, ['figure_pieces/fig3_habRasterEx'  '_' num2str(ex_idx) '.fig'])
end

%% example curve trajectories

load(pipeline_out)

for ex_idx = 1:length(curve_bird)
    ex_bird = curve_bird(ex_idx);
    ex_unit = curve_unit(ex_idx);
    ex_song = curve_song(ex_idx);
    
    switch ex_bird
        case 1
            these_units = 1:num_goods(1);
        case 2
            these_units = num_goods(1)+1:num_goods(1)+num_goods(2);
        case 3
            these_units = num_goods(2)+1:sum(num_goods);
    end
    this_unit = these_units(ex_unit);
    
    relevant_motif_chunks = [1 5; 6 25; 26 120; 121 125;];
    
    line_colors = {'r-', 'm-', 'b-', 'k-'};
    
    eval(['this_resp_vector = big_rv' num2str(ex_song) ';'])
    this_resp_piece = squeeze(this_resp_vector(this_unit,:,:));
%     this_resp_piece = this_resp_piece ./ max(this_resp_piece);
    
    for lineIdx = 1:4
        impt_motifs = relevant_motif_chunks(lineIdx,:);
        good_motifs = this_resp_piece(impt_motifs(1):impt_motifs(2),:);
        this_line = mean(good_motifs,1);
        max_line_val = max(max_line_val, max(this_line));
        
        plot(this_line, line_colors{lineIdx})
        if lineIdx == 1
            hold on;
        elseif lineIdx == 4
            hold off;
        end
    end
    
    title(['Bird ' num2str(ex_bird) ', Unit ' num2str(ex_unit) ', Song' num2str(songIdx)])
    ylabel('Resp. Str.')
    xlabel('Time')
    
    h5 = gcf;
    print(['figure_pieces/fig3_habCurveEx' '_' num2str(ex_idx)], '-dsvg', '-r300')
    saveas(h5, ['figure_pieces/fig3_habCurveEx'  '_' num2str(ex_idx) '.fig'])
end
        
%% gain model calculation and plotting

% goodness of fit, calculated as the max of the xcorr of the two sequences.
% then, also take peak lag (and track those just in case)
%
% Still tracking best habituation? 

%%
gain_vals = {};
err_vals = {};
bof_vals = {};
lag_vals = {};
for song_idx = 1:4
    eval(['this_rv = big_rv' num2str(song_idx) ';'])
    rv_inits = squeeze(mean(this_rv(:,1:5,:),2));
    rv_ends = squeeze(mean(this_rv(:,121:125,:),2));
    
    rv_long = [rv_inits rv_ends];
    max_r = max(rv_long,[],2);
    rv_inits = rv_inits ./ max_r;
    rv_ends = rv_ends ./ max_r;
    
    gain_checks = 0:0.01:3;
    
    this_gain_vals = [];
    this_err_vals = [];
    
    this_bof_vals = [];
    this_lag_vals = [];
    
    for unit_idx = 1:length(big_rv1(:,1,1))
        c1 = rv_inits(unit_idx,:);
        c2 = rv_ends(unit_idx,:);

        zero_lag = length(c1);%+1;
        
        auto1 = xcorr(c1);
        auto2 = xcorr(c2);
        auto_prod = auto1(zero_lag)*auto2(zero_lag);
        norm_factor = 1/sqrt(auto_prod);
        xc = xcorr(c1, c2) .* norm_factor;
        bof = xc(zero_lag);
        [bof,lag] = max(xc);
        lag = lag - zero_lag;
        
        
        best_err = 100; % maximum conceivable is just number of timesteps, <100
        best_gain = 0;
        for gain_idx = 1:length(gain_checks)
            c1x = c1*gain_checks(gain_idx);
            this_err = sum((c1x - c2).^2);
            
            if this_err < best_err
                best_err = this_err;
                best_gain = gain_checks(gain_idx);
            end
        end
        
        % great for visualizing individual examples
%         figure(9999);
%         plot(c1, 'r-')
%         hold on;
%         plot(c2, 'k-')
%         plot(c1*best_gain, 'r--')
%         hold off;
        
        if best_err<100
            this_gain_vals = [this_gain_vals; best_gain];
            this_err_vals = [this_err_vals; best_err];
            
            this_bof_vals = [this_bof_vals; bof];
            this_lag_vals = [this_lag_vals; lag];
        end

%         this_bof_vals = [this_bof_vals; bof];
%         this_lag_vals = [this_lag_vals; lag];

    end
    gain_vals{end+1} = this_gain_vals;
    err_vals{end+1} = this_err_vals;
    
    
    bof_vals{end+1} = this_bof_vals;
    lag_vals{end+1} = this_lag_vals;

    figure(91919)
    subplot(2,2,song_idx)
    histogram(this_bof_vals)%, 'BinEdges', 0:0.25:15)
    title('X-Corr Values')
    
    figure(91920)
    subplot(2,2,song_idx)
    histogram(this_gain_vals)
    title('Gain Values')
    
    % just for visualization - not included in dissertation
    figure(91921+song_idx)
    scatterhist(this_bof_vals, this_gain_vals,'Marker','.')
    title('Comparisons')
    xlabel('Xcorr Values')
    ylabel('Optimal Gain')
end

figure(91919)
h5 = gcf;
print(['figure_pieces/fig3_xcorrVals'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig3_xcorrVals.fig'])

figure(91920)
h5 = gcf;
print(['figure_pieces/fig3_gainVals'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig3_gainVals.fig'])
