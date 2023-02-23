% 1. load multibird summary data
% 2. plot distributions of # of error responses per neuron
% 3. plot likelihood of sensitivity for syll types and pert types
% 4. plot indices (effect mag) for different types (syll and pert)

% NB: missing stat tests, which would be nice because some trends seem
% significant - probably just binomial?
% could possibly show that mag distributions are not significantly
% different

clear
close all
error_detection_pipeline = 'D:\error_expmts\multibird_summary_ncm_errors.mat';

rng(23) % important for simulating number of pert responses per neuron -NOT USED
%% setup 

load(error_detection_pipeline)

% these permutations will be useful
% sort by ptype and then by intra v whole syll pert
pType_perm = [1 7 13 2 8 14 3 9 15 4 10 16 5 11 17 6 12 18];

% sort by ptype and then by H v N v S
pKind_perm = [1 7 13 4 10 16 2 8 14 5 11 17 3 9 15 6 12 18];

% get magnitude for each pert, learned and oddball
a1_total = mean(big_resp_mat(:, 1:25),2);
odd_total = mean(big_resp_mat(:,26:50),2);
a2_total = mean(big_resp_mat(:,51:75),2);

err_idx = a2_total - ((a1_total +  odd_total)/2);
odd_idx = odd_total - ((a1_total + a2_total)/2);

ei_sig = err_idx(big_fli_ids);
oi_sig = odd_idx(big_rdi_ids);

%% get sim numbers - what if error responses uniformly distributed? - NOT USED!
nrns_per_bird = [173; 106; 120; 87; 85];
bStart = [1; 174; 280; 400; 487];
bEnd = [173; 279; 399; 486; 571];

total_nrns = sum(nrns_per_bird);

brdi_blank = rand(100,size(big_rare_dev_iti,1), size(big_rare_dev_iti,2));
bfli_blank = rand(100,size(big_fam_late_iti,1), size(big_fam_late_iti,2));

brdi_sim = brdi_blank;
bfli_sim = bfli_blank;
for bird_idx = 1:5
    these_ids = bStart(bird_idx):bEnd(bird_idx);
    brdi_sim(:,:,these_ids) = brdi_blank(:,:,these_ids) < ...
        sum(sum(big_rare_dev_iti(:,these_ids)) / ...
        numel(big_rare_dev_iti(:,these_ids)));
    
    bfli_sim(:,:,these_ids) = bfli_blank(:,:,these_ids) < ...
        sum(sum(big_fam_late_iti(:,these_ids)) / ...
        numel(big_fam_late_iti(:,these_ids)));
end

% brdi_sim = brdi_blank < (sum(sum(big_rare_dev_iti))/numel(big_rare_dev_iti));
% bfli_sim = bfli_blank < (sum(sum(big_fam_late_iti))/numel(big_fam_late_iti));

sim_edges = -0.5:1:18.5;

rpn_rd = zeros(100,19);
rpn_fl = zeros(100,19);
for sim_idx = 1:100
    this_brdi_sim = squeeze(brdi_sim(sim_idx,:,:));
    this_bfli_sim = squeeze(bfli_sim(sim_idx,:,:));
    
    this_rpn_rd = sum(this_brdi_sim,1);
    this_rpn_fl = sum(this_bfli_sim,1);
    
    this_hc_rd = histcounts(this_rpn_rd,sim_edges);
    this_hc_fl = histcounts(this_rpn_fl,sim_edges);
    
    rpn_rd(sim_idx,:) = this_hc_rd;
    rpn_fl(sim_idx,:) = this_hc_fl;
end
%% distributions of # of learned error responses 

num_neurons = length(big_resp_mat(:, 1))/18;


resp_per_nrn = [];
sig_nrn_ids = [];
for nIdx = 1:num_neurons
    start_idx = ((nIdx-1)*18)+1;
    end_idx = start_idx + 17;
    
    num_resp = length(intersect(find(big_fli_ids >= start_idx), find(big_fli_ids <= end_idx)));
    
    sig_nrn_bool = ~isempty(intersect(find(big_fli_ids >= start_idx), find(big_fli_ids <= end_idx)));
    if sig_nrn_bool
        sig_nrn_ids = [sig_nrn_ids; nIdx];
    end
    
    resp_per_nrn = [resp_per_nrn; num_resp];
end

% sim results
mean_sim = mean(rpn_fl,1); % average across sims
sd_sim = std(rpn_fl,0,1); % stds

figure(81);
histogram(resp_per_nrn(sig_nrn_ids))
% hold on;
% plot(mean_sim(2:19));
% hold off;
title('Learned Error Responses per Neuron')

h5 = gcf;
print('figure_pieces/fig5_errorsPerNeuron', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_errorsPerNeuron' '.fig'])

%% distributions of # of oddball error responses 

num_neurons = length(big_resp_mat(:, 1))/18;

resp_per_nrn = [];
sig_nrn_ids = [];
for nIdx = 1:num_neurons
    start_idx = ((nIdx-1)*18)+1;
    end_idx = start_idx + 17;
    
    num_resp = length(intersect(find(big_rdi_ids >= start_idx), find(big_rdi_ids <= end_idx)));
    
    sig_nrn_bool = ~isempty(intersect(find(big_rdi_ids >= start_idx), find(big_rdi_ids <= end_idx)));
    if sig_nrn_bool
        sig_nrn_ids = [sig_nrn_ids; nIdx];
    end
    
    resp_per_nrn = [resp_per_nrn; num_resp];
end

% sim results
mean_sim = mean(rpn_rd,1); % average across sims
sd_sim = std(rpn_rd,0,1); % stds

figure(82);
histogram(resp_per_nrn(sig_nrn_ids))
title('Oddball Error Responses per Neuron')
% hold on;
% plot(mean_sim(2:19));
% hold off;

h5 = gcf;
print('figure_pieces/fig5_oddErrorsPerNeuron', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_oddErrorsPerNeuron' '.fig'])


%% distributions of sensitivity by syllable type of learned error responses 

bfli_pTypes = sum(big_fam_late_iti,2);
bfli_neurons = big_fam_late_iti(:,familiar_late_iti);

pfli_bars = [sum(bfli_pTypes(1:6)), sum(bfli_pTypes(7:12)), sum(bfli_pTypes(13:18))];
figure(83);
subplot(311)
bar(pfli_bars)
title('Fam. Late - Syllable Type Sorting (b, d, f)')

subplot(312)
bar([sum(bfli_pTypes(pType_perm(1:9))), sum(bfli_pTypes(pType_perm(10:18)))]);
title('Fam. Late - Pert. Length Sorting (intra v whole)')

subplot(313)
bar([sum(bfli_pTypes(pKind_perm(1:6))), sum(bfli_pTypes(pKind_perm(7:12))),...
    sum(bfli_pTypes(pKind_perm(13:18)))])
title('Fam. Late - Pert. Type Sorting (HS, N, Sil)')

h5 = gcf;
print('figure_pieces/fig5_flErrorTypes', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_flErrorTypes' '.fig'])


%% distributions of sensitivity by syllable type of oddball error responses 

brdi_pTypes = sum(big_rare_dev_iti,2);
brdi_neurons = big_rare_dev_iti(:,rare_dev_iti);

prdi_bars = ([sum(brdi_pTypes(1:6)), sum(brdi_pTypes(7:12)), sum(brdi_pTypes(13:18))]);
figure(84);
subplot(311)
bar(prdi_bars)
title('Rare Dev. - Syllable Type Sorting (b, d, f)')

subplot(312)
bar([sum(brdi_pTypes(pType_perm(1:9))), sum(brdi_pTypes(pType_perm(10:18)))]);
title('Rare Dev. - Pert. Length Sorting (intra v whole)')

subplot(313)
bar([sum(brdi_pTypes(pKind_perm(1:6))), sum(brdi_pTypes(pKind_perm(7:12))),...
    sum(brdi_pTypes(pKind_perm(13:18)))])
title('Rare Dev. - Pert. Type Sorting (HS, N, Sil)')

h5 = gcf;
print('figure_pieces/fig5_rdErrorTypes', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_rdErrorTypes' '.fig'])


%% get effect magnitudes for learned and oddball error responses
% first get effect magnitudes by perturbation

ei_mags_by_p = {};
oi_mags_by_p = {};

fl_idx = 0;
rd_idx = 0;
for p_idx = 1:18
    this_row_fl = big_fam_late_iti(p_idx,:);
    this_row_rd = big_rare_dev_iti(p_idx,:);
    
    num_fls = length(find(this_row_fl==1));
    fl_new = fl_idx + num_fls;
    if num_fls > 0
        ei_mags_by_p{p_idx} = ei_sig(fl_idx+1:fl_new);
    end
    fl_idx = fl_new;
    
    num_rds = length(find(this_row_rd==1));
    rd_new = rd_idx + num_rds;
    if num_rds > 0
        oi_mags_by_p{p_idx} = oi_sig(rd_idx+1:rd_new);
    end
    rd_idx = rd_new;
end
    
%% plot learned error effect magnitudes (syllable)

% syll type
syll_rows = [1:6 ; 7:12 ; 13:18];
fl_mag_syll_dist = cell(3,1);
for syll_idx = 1:3
    for row_idx = syll_rows(syll_idx,:)
        fl_mag_syll_dist{syll_idx} =...
            [fl_mag_syll_dist{syll_idx}; ei_mags_by_p{row_idx}];
    end
end


edges = [0:0.05:1];
for syll_idx = 1:3
    figure(86);
    subplot(3,1,syll_idx)
    histogram(fl_mag_syll_dist{syll_idx},edges)
    if syll_idx == 1
        title('Fam. Late - Syllable Type Sorting (b, d, f)')
    end
end

h5 = gcf;
print('figure_pieces/fig5_flErrorMags', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_flErrorMags' '.fig'])


%% plot learned error effect magnitudes (error type (intra v whole))

% pert type
pType_rows = [1 7 13 2 8 14 3 9 15; 4 10 16 5 11 17 6 12 18];
fl_mag_pType_dist = cell(2,1);
for pType_idx = 1:2
    for row_idx = pType_rows(pType_idx,:)
        fl_mag_pType_dist{pType_idx} =...
            [fl_mag_pType_dist{pType_idx}; ei_mags_by_p{row_idx}];
    end
end


edges = [0:0.05:1];
for pType_idx = 1:2
    figure(87);
    subplot(2,1,pType_idx)
    histogram(fl_mag_pType_dist{pType_idx},edges)
    if pType_idx == 1
        title('Fam. Late - Pert. Type Sorting (intra v whole)')
    end
end

h5 = gcf;
print('figure_pieces/fig5_flErrorMags_type', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_flErrorMags_type' '.fig'])


%% plot learned error effect magnitudes (pert kind - harm stack v noise v silence)

% sort by ptype and then by H v N v S
pKind_rows = [1 7 13 4 10 16; 2 8 14 5 11 17; 3 9 15 6 12 18];
fl_mag_pKind_dist = cell(3,1);
for pKind_idx = 1:3
    for row_idx = pKind_rows(pKind_idx,:)
        fl_mag_pKind_dist{pKind_idx} =...
            [fl_mag_pKind_dist{pKind_idx}; ei_mags_by_p{row_idx}];
    end
end


edges = [0:0.05:1];
for pKind_idx = 1:3
    figure(88);
    subplot(3,1,pKind_idx)
    histogram(fl_mag_pKind_dist{pKind_idx},edges)
    if pKind_idx == 1
        title('Fam. Late - Pert. Kind Sorting (H v N v S)')
    end
end

h5 = gcf;
print('figure_pieces/fig5_flErrorMags_kind', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_flErrorMags_kind' '.fig'])

%% plot oddball error effect magnitudes (syllable)

% syll type
syll_rows = [1:6 ; 7:12 ; 13:18];
rd_mag_syll_dist = cell(3,1);
for syll_idx = 1:3
    for row_idx = syll_rows(syll_idx,:)
        rd_mag_syll_dist{syll_idx} =...
            [rd_mag_syll_dist{syll_idx}; oi_mags_by_p{row_idx}];
    end
end


edges = [0:0.05:1];
for syll_idx = 1:3
    figure(89);
    subplot(3,1,syll_idx)
    histogram(rd_mag_syll_dist{syll_idx},edges)
    if syll_idx == 1
        title('Rare Dev. - Syllable Type Sorting (b, d, f)')
    end
end

h5 = gcf;
print('figure_pieces/fig5_rdErrorMags', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_rdErrorMags' '.fig'])


%% plot oddball error effect magnitudes (error type (intra v whole))

% pert type
pType_rows = [1 7 13 2 8 14 3 9 15; 4 10 16 5 11 17 6 12 18];
rd_mag_pType_dist = cell(2,1);
for pType_idx = 1:2
    for row_idx = pType_rows(pType_idx,:)
        rd_mag_pType_dist{pType_idx} =...
            [rd_mag_pType_dist{pType_idx}; oi_mags_by_p{row_idx}];
    end
end


edges = [0:0.05:1];
for pType_idx = 1:2
    figure(90);
    subplot(2,1,pType_idx)
    histogram(rd_mag_pType_dist{pType_idx},edges)
    if pType_idx == 1
        title('Rare Dev. - Pert. Type Sorting (intra v whole)')
    end
end

h5 = gcf;
print('figure_pieces/fig5_rdErrorMags_type', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_rdErrorMags_type' '.fig'])

%% plot oddball error effect magnitudes (pert kind - h stck v noise v silence)

% sort by ptype and then by H v N v S
pKind_rows = [1 7 13 4 10 16; 2 8 14 5 11 17; 3 9 15 6 12 18];
rd_mag_pKind_dist = cell(3,1);
for pKind_idx = 1:3
    for row_idx = pKind_rows(pKind_idx,:)
        rd_mag_pKind_dist{pKind_idx} =...
            [rd_mag_pKind_dist{pKind_idx}; oi_mags_by_p{row_idx}];
    end
end


edges = [0:0.05:1];
for pKind_idx = 1:3
    figure(91);
    subplot(3,1,pKind_idx)
    histogram(rd_mag_pKind_dist{pKind_idx},edges)
    if pKind_idx == 1
        title('Rare Dev. - Pert. Kind Sorting (H v N v S)')
    end
end

h5 = gcf;
print('figure_pieces/fig5_rdErrorMags_kind', '-dsvg', '-r300')
saveas(h5, ['figure_pieces/fig5_rdErrorMags_kind' '.fig'])

%%
disp('all done')
