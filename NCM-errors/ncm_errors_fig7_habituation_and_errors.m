% 1. load big error pipeline output
% 2. cartoon of complementary habituation and error emergence curves v actual 
% 3. Scatter plots of error v habituation idxs, oddball v habituation idxs
% 4. network-wide habituation during each of 4 phases
%
% nb: habituation index can be less than 1 bc big_resp_mat values are
% between -1 and +1 (since all responses are corrected by iti response rate)
%
% point is to illustrate insufficiency of simple matched model of error
% emergence, while maintaining the possibility that the two processes are
% related, due to habituation over the course of the experiment (and of
% course it's also possible that habituation triggers some un-observable
% long-term process that only completes after habituation phase).

% add 1 more plot - imagesc of activity during the first period, then
% during the second, third and fourth periods (standards only)

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
load(error_detection_pipeline)

%% model cartoon (nb: these don't necessarily reflect real data)
x = 0:100;
decay_model = exp(-x/16);
decay_actual = exp(-x/4);
growth_model = 1 ./ (1 + exp(-0.2*(x-25)));
growth_actual = zeros(101,1);
growth_actual(52:101) = 1;

figure(1); 
subplot(211)
plot(x, decay_model./max(decay_model), 'linewidth', 4); 
hold on; 
plot(x, growth_model, 'linewidth', 4); hold off;
xlabel('Time')
ylabel('Response Strength')
legend('SSA', 'Error Sensitivity', 'Location', 'Best')
title('Matched Model of Error Emergence')
xticks([])

figure(1);
subplot(212)
plot(x, decay_actual./max(decay_actual), 'linewidth', 4); 
hold on; 
plot(x, growth_actual, 'linewidth', 4); hold off;
xlabel('Time')
ylabel('Response Strength')
legend('SSA', 'Error Sensitivity', 'Location', 'Best')
title('Experimental Data')
xticks([])

h5 = gcf;
print('figure_pieces/fig7_timescale_cartoon', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig7_timescale_cartoon' '.fig'])

%% pre-processing to id significant oddball and error responses
% get magnitude for each pert, learned and oddball
a1_start = mean(big_resp_mat(:,1:5),2);
a1_end = mean(big_resp_mat(:,21:25),2);
a2_end = mean(big_resp_mat(:,71:75),2);

a1_total = mean(big_resp_mat(:, 1:25),2);
odd_total = mean(big_resp_mat(:,26:50),2);
a2_total = mean(big_resp_mat(:,51:75),2);

hab_idx = a1_start - a2_end;
odd_idx = odd_total - ((a1_total + a2_total)/2);
err_idx = a2_total - ((a1_total +  odd_total)/2);

ei_sig = err_idx(big_fli_ids);
oi_sig = odd_idx(big_rdi_ids);

num_neurons = length(big_resp_mat(:, 1))/18;

resp_per_nrn_fl = [];
sig_nrn_ids_fl = [];
resp_per_nrn_rd = [];
sig_nrn_ids_rd = [];
ei2 = []; oi2 = []; hi2 = [];
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
    
    hi2(nIdx) = max(hab_idx(start_idx:end_idx));
end

%% scatter plots of learned and oddball error idx v habituation 

figure(5)
subplot(2,1,1)
plot(hi2, ei2, 'k.', 'MarkerSize', 12)
hold on;
plot(hi2(sig_nrn_ids_fl), ei2(sig_nrn_ids_fl), 'r.', 'MarkerSize', 12)
hold off;
xlabel('Hab. Idx.')
ylabel('Err. Idx.')
axis([-0.5 1.5 -0.5 0.75])
title('Learned Error Idx versus Habituation Index - max per neuron')

subplot(2,1,2)
plot(hi2, oi2, 'k.', 'MarkerSize', 12)
hold on;
plot(hi2(sig_nrn_ids_rd), oi2(sig_nrn_ids_rd), 'r.', 'MarkerSize', 12)
hold off;
xlabel('Hab. Idx.')
ylabel('Odd. Idx.')
axis([-0.5 1.5 -0.5 0.75])
title('Oddball Error Idx versus Habituation Index - max per neuron')


h5 = gcf;
print('figure_pieces/fig7_hab_scatters', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig7_hab_scatters' '.fig'])
%% habituation pre-processing - need to get all standard responses
% load individual bird

std_ids_bl1 = 1:19:475;
std_ids_pert = find(strcmp(motifInfo(476:2495,2),'std_A')); %motifInfo always same
std_ids_hab = 2496:2895;
std_ids_bl2 = 2896:19:3370;

big_std_bl1 = [];
big_std_pert = [];
big_std_hab = [];
big_std_bl2 = [];
for bird_id = 1:length(the_good_birds)
    this_bird = the_good_birds{bird_id};
    load(this_bird)
    
    % are motif_resp elements normalized to the ITI at all? NO (fixable if issue)
    big_std_bl1 = [big_std_bl1; motif_resp(std_ids_bl1,good_unit_ids)'];
    big_std_pert = [big_std_pert; motif_resp(std_ids_pert,good_unit_ids)'];
    big_std_hab = [big_std_hab; motif_resp(std_ids_hab,good_unit_ids)'];
    big_std_bl2 = [big_std_bl2; motif_resp(std_ids_bl2,good_unit_ids)'];
end


%% habituation across all phases

% for every good neuron across all datasets, we want the response to every
% single standard motif. we'll then calculate a habituation index for bl1,
% pert, hab, and bl2 phases by calculating 
% (mean_init - mean_final) ./ (mean_init + mean_final)

hab_bl1 = (mean(big_std_bl1(:,1:5),2) - mean(big_std_bl1(:,end-4:end),2)) ./ ...
    (mean(big_std_bl1(:,1:5),2) + mean(big_std_bl1(:,end-4:end),2));

hab_pert = (mean(big_std_pert(:,1:5),2) - mean(big_std_pert(:,end-4:end),2)) ./ ...
    (mean(big_std_pert(:,1:5),2) + mean(big_std_pert(:,end-4:end),2));

hab_hab = (mean(big_std_hab(:,1:5),2) - mean(big_std_hab(:,end-4:end),2)) ./ ...
    (mean(big_std_hab(:,1:5),2) + mean(big_std_hab(:,end-4:end),2));

hab_bl2 = (mean(big_std_bl2(:,1:5),2) - mean(big_std_bl2(:,end-4:end),2)) ./ ...
    (mean(big_std_bl2(:,1:5),2) + mean(big_std_bl2(:,end-4:end),2));


figure(245);
subplot(141)
histogram(hab_bl1)
title('Baseline Phase 1')

subplot(142)
histogram(hab_pert)
title('Pert. Phase')

subplot(143)
histogram(hab_hab)
title('Habituation Phase')

subplot(144)
histogram(hab_bl2)
title('Baseline Phase 2')

h5 = gcf;
print('figure_pieces/fig7_hab_across_phases', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig7_hab_across_phases' '.fig'])

%%
disp('all done')
