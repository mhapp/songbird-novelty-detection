%% IBI experiment results
% distributions again
% summaries showing significance

%% prep
clear
close all
pipeline_out = 'D:\ibi_expmts\multibird_summary_ibi.mat';

%% load multibird summary
load(pipeline_out)

%% summary
% final shape of big_ibi_tensor is (song id X unit num X rep id X ibi id)

rebound_units = cell(2,1);
avg_rebound_strength = cell(2,1);

key_motifs = [1 2751; 1251 2801];
xvals = [.05, .1, .150, .2, .250, .3, .4, .5, .75, 1];



this_uType_str = 'su';
these_good_units = sum(num_goods);
these_motif_resp = total_motif_resp;
this_tensor_raw = big_ibi_tensor;
this_tensor_frac = big_ibi_frac;

for songIdx = 1:2 % only first two songs are important
    tkm = key_motifs(songIdx,:);
    
    these_frac_inc_values = zeros(length(these_good_units), 250); % if error, numBouts is not 25
    these_raw_inc_values = zeros(length(these_good_units), 250);
    these_dev_values = zeros(length(these_good_units), 250);
    
    bout_nums = [ceil(tkm(1)/5) ceil(tkm(1)/5)+249];
    
    for unitIdx = 1:sum(num_goods)
%         disp([this_uType_str num2str(unitIdx) ' of '
%         num2str(sum(num_goods)) ', Song ' num2str(songIdx)]) % helpful
%         for testing
        
        this_resp_p1 = total_motif_resp(tkm(1):tkm(1)+1249, unitIdx);
        this_resp_p2 = total_motif_resp(tkm(2):tkm(2)+49, unitIdx);
        
        max_p1 = max(this_resp_p1);
        
        stacked_p1 = reshape(this_resp_p1, [5, 250]);
        stacked_p2 = reshape(this_resp_p2, [5, 10]);
        
        zero_starts_p1 = find(stacked_p1(1,:) == 0);
        zero_starts_p2 = find(stacked_p2(1,:) == 0);
        
        
        if unitIdx == 1
            eval(['resp1_tensor_s' num2str(songIdx) ' = ' ...
                ' zeros(5, 250, length(these_good_units));']);
            eval(['resp2_tensor_s' num2str(songIdx) ' = ' ...
                ' zeros(5, 10, length(these_good_units));']);
        end
        eval(['resp1_tensor_s' num2str(songIdx) '(:,:,unitIdx) = ' ...
            ' stacked_p1;']);
        eval(['resp2_tensor_s' num2str(songIdx) this_uType_str '(:,:,unitIdx) = ' ...
            ' stacked_p2;']);
        
        mean_ib_p1s = mean(stacked_p1, 2);
        std_ib_p1s = std(stacked_p1, [], 2);
        
        mean_ib_p2s = mean(stacked_p2, 2);
        std_ib_p2s = std(stacked_p2, [], 2);
        
        num_bouts = length(this_resp_p1)/5;
        these_raw_ib_rbds = zeros(num_bouts,1);
        these_frac_ib_rbds = zeros(num_bouts,1);
        for bIdx = 1:num_bouts
            this_first_motif = (bIdx-1)*5 + 1;
            if bIdx > 1
                this_resp = this_resp_p1(this_first_motif:this_first_motif+4);
                last_resp = this_resp_p1(this_first_motif-5:this_first_motif-1);
                
                raw_ib_rbd = this_resp(1) - last_resp(end);
                frac_ib_rbd = raw_ib_rbd / max_p1;
                
                these_raw_ib_rbds(bIdx) = raw_ib_rbd;
                these_frac_ib_rbds(bIdx) = frac_ib_rbd;
                
                
                these_frac_inc_values(unitIdx, bIdx) = frac_ib_rbd;
                these_raw_inc_values(unitIdx, bIdx) = raw_ib_rbd;
            end
        end
        
        
        avg_raw_tensor = mean(squeeze(this_tensor_raw(songIdx,unitIdx,:,:)),1);
        std_raw_tensor = std(squeeze(this_tensor_raw(songIdx,unitIdx,:,:)),0,1) ./ sqrt(size(this_tensor_raw,3));

        end_dist = squeeze(this_tensor_frac(songIdx,unitIdx,:,end));
        rbd_sig = ttest(end_dist,0,'Tail','right');
        rbd_sig(isnan(rbd_sig)) = 0;
        if rbd_sig
%             disp('rbd sig!') % helpful for testing
            rebound_units{songIdx}(end+1) = unitIdx;
            avg_rebound_strength{songIdx} = [avg_rebound_strength{songIdx}; avg_raw_tensor]; % rows are units
        end
        
    end
    
    these_IBIs = xvals;
    
    %
    rbd_avg_across_trials = mean(squeeze(this_tensor_frac(songIdx,rebound_units{songIdx},:,:)),2);
    this_y = squeeze(mean(rbd_avg_across_trials,1));
    this_err = std(rbd_avg_across_trials, [], 1) ./ sqrt(length(rebound_units{songIdx}));
        
    XX = [ones(length(these_IBIs),1) these_IBIs'];
    b = XX\this_y;
    yCalc = XX*b;
    
    rsq = 1 - sum((this_y - yCalc).^2) / sum((this_y - mean(this_y)).^2);
    
    figure(2345)
    subplot(1,2,songIdx)
    errorbar(these_IBIs, this_y, 2*squeeze(this_err), 'o')
    hold on;
    plot(these_IBIs, yCalc, 'r-')
    hold off;
    title(['Rebounds by IBI, Song ' num2str(songIdx) ', ' this_uType_str])
    xlabel('Inter-bout Interval')
    ylabel('Rebound (prop. of max.)')
    set(gca, 'FontSize', 24)
    xlim([0 1.1])
    
    disp(['slope is ' num2str(b(2)) ' of max per sec.'])
    disp(['rsq is ' num2str(rsq) '.'])
end

h5 = figure(2345);
print(['figure_pieces/fig7_ibi_summary_frac_v2'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig7_ibi_summary_frac_v2.fig'])


%% distributions of strengths at 1s ibi

bin_edges = -20:4:60;

for song_idx = 1:2
    these_long_rbds = squeeze(big_ibi_tensor(song_idx,:,:,10));
    these_avg_long_rbds = mean(these_long_rbds, 2);
    
    this_sig = rebound_units{song_idx};
    
    figure(999)
    subplot(1,2,song_idx)
    histogram(these_avg_long_rbds,'BinEdges', bin_edges)
    hold on;
    histogram(these_avg_long_rbds(this_sig),'BinEdges', bin_edges)
    hold off;
    xlabel('Rebound (Hz)')
    ylabel('Number of neurons')
    
    title(['Rebounds, Song ' num2str(song_idx)])
end

h5 = figure(999);
print(['figure_pieces/fig7_ibi_rbd_dist_raw'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig7_ibi_rbd_dist_raw.fig'])

%% distributions of strengths at 1s ibi

bin_edges = -0.8:0.1:0.8;

for song_idx = 1:2
    these_long_rbds = squeeze(big_ibi_frac(song_idx,:,:,10));
    these_avg_long_rbds = mean(these_long_rbds, 2);
    
    this_sig = rebound_units{song_idx};
    
    figure(999)
    subplot(1,2,song_idx)
    histogram(these_avg_long_rbds,'BinEdges', bin_edges)
    hold on;
    histogram(these_avg_long_rbds(this_sig),'BinEdges', bin_edges)
    hold off;
    xlabel('Rebound (frac. of max.)')
    ylabel('Number of neurons')
    
    title(['Rebounds, Song ' num2str(song_idx)])
end

h5 = figure(999);
print(['figure_pieces/fig7_ibi_rbd_dist_frac'], '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig7_ibi_rbd_dist_frac.fig'])