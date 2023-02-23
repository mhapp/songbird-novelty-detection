% 1. load waveform data (possibly from separate pipeline)
% 2. perform 1d separation of fs v rs neurons using end slope
% 3. distribution of error sensitivity, fs v rs (learned & oddball)
% 4. oddball index, fs v rs
% 5. learned error index, fs v rs
%
% message is: "this distinction doesn't really matter," which further
% supports central claim that NCM is so densely interconnected that
% everything is multiplexed, pointing towards a dense model for error
% detection, rather than a sparse one

% at the very least, this shows that circuit-breaking approaches used when
% studying mammalian systems won't really apply to NCM

%% setup
clear
close all

the_good_birds = {'D:\error_expmts\pipeline_outputs\pipeline_7635_210723_LH_NCM_g0.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_9007_210722_LH_NCM_g0.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_8975_210610_LH_NCM_g1.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_9047_210722_LH_NCM_g0.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_9049_210724_LH_NCM_g0.mat'
    };

error_detection_pipeline = 'D:\error_expmts\multibird_summary_ncm_errors.mat';

rng(23) % just affects the sorting of summary fig

%% get waveform data across birds and process it somewhat
big_end_slopes = [];
big_wfs = [];

for bird_id = 1:length(the_good_birds)
    this_bird = the_good_birds{bird_id};
    load(this_bird)
    load(['D:/' bird_id '_wfData_' rec_date '.mat'])
    
    end_slopes = [];
    scaled_wfs = [];
    for wf_idx = 1:size(best_wfs,1) % only for good units - no problem
        twf = best_wfs(wf_idx,:);
        [~,this_trough] = min(twf(26:end-25)); % should be in middle...
        this_trough = this_trough + 25;
        this_chopped = twf(this_trough-25:this_trough+25);
        dmWF = this_chopped - mean(this_chopped(1:10));
        dmWF = dmWF ./ abs(min(dmWF));
        dwf = diff(dmWF);
        
        scaled_wfs = [scaled_wfs; dmWF];
        
        this_end = (dmWF(31) - dmWF(26)) / 5; 
        end_slopes = [end_slopes; this_end];
    end
    
    % then we can look at patterns outside of loop
    big_wfs = [big_wfs; scaled_wfs];
    big_end_slopes = [big_end_slopes; end_slopes];
end

% also load error pipeline output
load(error_detection_pipeline)

disp('data is loaded!')

%% perform 1d separation

wf_types = kmeans(big_end_slopes,2);

fs_ids = find(wf_types==2); % larger slopes
rs_ids = find(wf_types==1); % smaller slopes

edges = 0:0.006125:.5;
figure(12);
subplot(211)
histogram(big_end_slopes,edges)
title('End Slopes of Waveform')

subplot(212)
histogram(big_end_slopes(rs_ids),edges)
hold on;
histogram(big_end_slopes(fs_ids),edges)
hold off;
title('K-means Clustering')
legend('Putative RS Cells', 'Putative FS Cells', 'Location', 'Best')

h5 = gcf;
print('figure_pieces/fig6_fs_rs_split', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_fs_rs_split' '.fig'])

%% processing for error types, etc

% get magnitude for each pert, learned and oddball
a1_total = mean(big_resp_mat(:, 1:25),2);
odd_total = mean(big_resp_mat(:,26:50),2);
a2_total = mean(big_resp_mat(:,51:75),2);

err_idx = a2_total - ((a1_total +  odd_total)/2);
odd_idx = odd_total - ((a1_total + a2_total)/2);

ei_sig = err_idx(big_fli_ids);
oi_sig = odd_idx(big_rdi_ids);

num_neurons = length(big_resp_mat(:, 1))/18;

resp_per_nrn_fl = [];
sig_nrn_ids_fl = [];
resp_per_nrn_rd = [];
sig_nrn_ids_rd = [];
for nIdx = 1:num_neurons
    start_idx = ((nIdx-1)*18)+1;
    end_idx = start_idx + 17;
    
    % fl processing
    num_resp_fl = length(intersect(find(big_fli_ids >= start_idx), find(big_fli_ids <= end_idx)));
    sig_nrn_bool = ~isempty(intersect(find(big_fli_ids >= start_idx), find(big_fli_ids <= end_idx)));
    if sig_nrn_bool
        sig_nrn_ids_fl = [sig_nrn_ids_fl; nIdx];
    end
    resp_per_nrn_fl = [resp_per_nrn_fl; num_resp_fl];
    ei2(nIdx) = max(err_idx(start_idx:end_idx));
    
    % rd processing
    num_resp_rd = length(intersect(find(big_rdi_ids >= start_idx), find(big_rdi_ids <= end_idx)));
    sig_nrn_bool = ~isempty(intersect(find(big_rdi_ids >= start_idx), find(big_rdi_ids <= end_idx)));
    if sig_nrn_bool
        sig_nrn_ids_rd = [sig_nrn_ids_rd; nIdx];
    end
    resp_per_nrn_rd = [resp_per_nrn_rd; num_resp_rd];
    oi2(nIdx) = max(odd_idx(start_idx:end_idx));
    
end

%% probability of (learned, oddball) responses given (fs, rs)

% what is probability that cell has error given that it is fs v rs?
% learned errors
fs_learned_prop = length(intersect(sig_nrn_ids_fl, fs_ids))/length(fs_ids);
rs_learned_prop = length(intersect(sig_nrn_ids_fl, rs_ids))/length(rs_ids);
total_learned_prop = length(sig_nrn_ids_fl) / num_neurons;

% oddball errors
fs_oddball_prop = length(intersect(sig_nrn_ids_rd, fs_ids))/length(fs_ids);
rs_oddball_prop = length(intersect(sig_nrn_ids_rd, rs_ids))/length(rs_ids);
total_oddball_prop = length(sig_nrn_ids_rd) / num_neurons;


% what is probability that cell is fs v rs given that is has error?
learned_fs_prop = length(intersect(sig_nrn_ids_fl, fs_ids))/length(sig_nrn_ids_fl);
learned_rs_prop = length(intersect(sig_nrn_ids_fl, rs_ids))/length(sig_nrn_ids_fl);
total_rs_prop = length(rs_ids) / num_neurons;

% oddball errors
oddball_fs_prop = length(intersect(sig_nrn_ids_rd, fs_ids))/length(sig_nrn_ids_rd);
oddball_rs_prop = length(intersect(sig_nrn_ids_rd, rs_ids))/length(sig_nrn_ids_rd);
total_fs_prop = length(fs_ids) / num_neurons;




figure(14);
subplot(2,2,3)
bar([total_learned_prop fs_learned_prop rs_learned_prop])
title('Prop. of Neurons with Learned Errors (total, fs, rs)')

subplot(2,2,4)
bar([total_oddball_prop fs_oddball_prop rs_oddball_prop])
title('Prop. of Neurons with Oddball Errors (total, fs, rs)')


subplot(2,2,2)
bar([total_fs_prop learned_fs_prop oddball_fs_prop])
title('Prop. of Neurons that are FS (total, learned, error)')

subplot(2,2,1)
bar([total_rs_prop learned_rs_prop oddball_rs_prop])
title('Prop. of Neurons that are RS (total, learned, error)')

h5 = gcf;
print('figure_pieces/fig6_fs_rs_errors', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_fs_rs_errors' '.fig'])

%% response magnitudes, fs v rs

% first get sig ids for fs, rs
sig_err_rs = intersect(sig_nrn_ids_fl,rs_ids);
sig_odd_rs = intersect(sig_nrn_ids_rd,rs_ids);
sig_err_fs = intersect(sig_nrn_ids_fl,fs_ids);
sig_odd_fs = intersect(sig_nrn_ids_rd,fs_ids);

% learned error indices
[cts, edges] = histcounts(ei2(rs_ids));
[cts2, edges2] = histcounts(ei2(fs_ids));
figure(7);
subplot(211)
histogram(ei2(rs_ids), 'BinEdges', edges)
hold on;
histogram(ei2(sig_err_rs),'BinEdges',edges)
hold off;
title('Learned Error Indices - Max per Neuron - RS')

subplot(212)
histogram(ei2(fs_ids), 'BinEdges', edges2)
hold on;
histogram(ei2(sig_err_fs),'BinEdges',edges2)
hold off;
title('Learned Error Indices - Max per Neuron - FS')

h5 = gcf;
print('figure_pieces/fig6_fs_rs_flMags', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_fs_rs_flMags' '.fig'])


% oddball error indices
[cts, edges] = histcounts(oi2(rs_ids));
[cts2, edges2] = histcounts(oi2(fs_ids));
figure(8);
subplot(211)
histogram(oi2(rs_ids), 'BinEdges', edges)
hold on;
histogram(oi2(sig_odd_rs),'BinEdges',edges)
hold off;
title('Oddball Error Indices - Max per Neuron - RS')

subplot(212)
histogram(oi2(fs_ids), 'BinEdges', edges2)
hold on;
histogram(oi2(sig_odd_fs),'BinEdges',edges2)
hold off;
title('Oddball Error Indices - Max per Neuron - FS')

h5 = gcf;
print('figure_pieces/fig6_fs_rs_rdMags', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig6_fs_rs_rdMags' '.fig'])

%%
disp('All done')
