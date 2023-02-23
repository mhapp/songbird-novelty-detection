% it is very time-consuming to perform the analyses that determine whether
% neurons are responsive to perturbations of different types (as this
% involves over 50 distinct rank-sum calculations per neuron).
%
% to save time in analysis, we run this code one time after running the
% prelim pipeline on each bird. the result is large multi-bird data
% structures summarizing error-responsiveness data

the_good_birds = {'D:\error_expmts\pipeline_outputs\pipeline_7635_210723_LH_NCM_g0.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_9007_210722_LH_NCM_g0.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_8975_210610_LH_NCM_g1.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_9047_210722_LH_NCM_g0.mat',...
    'D:\error_expmts\pipeline_outputs\pipeline_9049_210724_LH_NCM_g0.mat'
    };

pipeline_output_save_target = 'D:\error_expmts\multibird_summary_ncm_errors.mat';

big_aud_bcs = [];
big_rare_dev_iti = [];
big_rare_dev_std = [];
big_fam_late_iti = [];
big_fam_late_std = [];

big_rdi_mat = [];
big_fli_mat = [];
big_std_mat = [];
big_rdi_ids =[];
big_fli_ids = [];
big_std_ids = [];
big_resp_mat = []; % we'll use this to calculate indices for errors (each row is neuron x pert combo)
brm_mean_norm = [];

big_adapt_resp_raw = [];
big_adapt_resp_bcs = [];

big_hab_adapt_raw = [];
big_odd_adapt_raw = [];

big_aud_resp_traj = [];
all_ids = 0;
total_ids = [];
for big_bird_num = 1:length(the_good_birds)
    clearvars -except big_bird_num the_good_birds big_aud_bcs ...
        big_rare_dev_iti big_fam_late_iti all_ids total_ids...
        big_rare_dev_std big_fam_late_std big_adapt_resp_bcs...
        big_adapt_resp_raw big_rdi_mat big_fli_mat big_std_mat...
        big_hab_adapt_raw big_odd_adapt_raw big_aud_resp_traj big_resp_mat...
        big_rdi_ids big_fli_ids big_std_ids pipeline_output_save_target ...
        brm_mean_norm
    close all
    
    this_bird = the_good_birds{big_bird_num};
    load(this_bird)
        
    disp(['starting bird ' num2str(big_bird_num)])
    
    
    %% sound feature sensitivity
    
    iti_resp_aud = iti_resp(1:5,good_unit_ids);
    iti_std = std(iti_resp_aud, [], 1);
    iti_mean = mean(iti_resp_aud, 1);
    
    iti_ceil_aud = iti_mean + 1.96*iti_std;
    iti_floor_aud = iti_mean - 1.96*iti_std;
    
    aud_resp_bc = zeros(length(unique_sylls),length(good_unit_ids));
    adapt_resp_bc = zeros(length(unique_sylls),length(good_unit_ids));
    for sIdx = 1:length(unique_sylls)
        this_syll = unique_sylls{sIdx};
        
        syll_resp = reshape(all_trial_resp,size(all_trial_resp,1)*size(all_trial_resp,2),size(all_trial_resp,3));
        syll_resp = syll_resp(:, good_unit_ids);
        
        these_trials_x = find(strcmp(syllableInfo(1:2375,3),this_syll));
        these_trials = these_trials_x(1:5);
        these_final_trials = these_trials_x(end-4:end);
        
        if sIdx >= 19
            if strcmp(stim_date, '210610')
                hab_trials_x = find(strcmp(syllableInfo(12476:14475,3),this_syll));
                odd_trials_x = find(strcmp(syllableInfo(2376:12475,3),this_syll));
                last_pre_hab = 12475;
                last_pre_odd = 2375;
            elseif strcmp(stim_date, '210613')
                hab_trials_x = find(strcmp(syllableInfo(12436:14435,3),this_syll));
                odd_trials_x = find(strcmp(syllableInfo(2376:12435,3),this_syll));
                last_pre_hab = 12435;
                last_pre_odd = 2375;
            end
            hab_trials = hab_trials_x(1:5);
            hab_final_trials = hab_trials_x(end-4:end);
            
            odd_trials = odd_trials_x(1:5);
            
            
            final_hab = syll_resp(last_pre_hab+hab_final_trials,:);
            mean_hab_final = mean(final_hab,1);
            these_hab = syll_resp(last_pre_hab+hab_trials,:);
            mean_hab_init = mean(these_hab,1);
            
            odd_init = syll_resp(last_pre_odd+odd_trials,:);
            mean_odd_init = mean(odd_init,1);
            
            hab_adapt_raw(sIdx,:) = (mean_hab_init - mean_hab_final) ./ (mean_hab_init + mean_hab_final);
            odd_adapt_raw(sIdx,:) = (mean_odd_init - mean_hab_final) ./ (mean_odd_init + mean_hab_final);
        end
        
        final_resp = syll_resp(these_final_trials,:);
        mean_final = mean(final_resp,1);
        ceil_final = mean_final + 1.96*std(final_resp,[],1);
        floor_final = mean_final - 1.96*std(final_resp,[],1);
        
        
        these_resp = syll_resp(these_trials,:);
        mean_resp = mean(these_resp, 1);
        this_std_resp = std(these_resp, [], 1);
        ceil_resp = mean_resp + 1.96*this_std_resp;
        floor_resp = mean_resp - 1.96*this_std_resp;
        
        this_resp_cond = (floor_resp > iti_ceil_aud) | (ceil_resp < iti_floor_aud);
        
        aud_resp_bc(sIdx,:) = this_resp_cond;
        
        adapt_resp_raw(sIdx,:) = (mean_resp - mean_final) ./ (mean_resp + mean_final);
        adapt_resp_bc(sIdx,:) = (floor_resp > ceil_final); % 1 if end and start don't overlap
        
        if sIdx >= 19
            big_aud_resp_traj = [big_aud_resp_traj; ...
                mean_resp' mean_final' mean_odd_init' mean_hab_init' mean_hab_final'];
        end
            
    end
    
    big_hab_adapt_raw = [big_hab_adapt_raw hab_adapt_raw];
    big_odd_adapt_raw = [big_odd_adapt_raw odd_adapt_raw];


    
    big_aud_bcs = [big_aud_bcs aud_resp_bc];
    
    big_adapt_resp_raw = [big_adapt_resp_raw adapt_resp_raw];
    big_adapt_resp_bcs = [big_adapt_resp_bcs adapt_resp_bc];
    
    big_rare_dev_iti = [big_rare_dev_iti rare_dev_bc_iti];
    big_rare_dev_std = [big_rare_dev_std rare_dev_bc_std];
    
    big_fam_late_iti = [big_fam_late_iti fam_late_bc_iti];
    big_fam_late_std = [big_fam_late_std fam_late_bc_std];
    
    
    std_aud_resp = find(sum(aud_resp_bc(19:23,:),1) > 0);
    
    %% big matrices of error responses
    
    sig_thresh = .05/18; % bonferroni-corrected sig val for detection
       
    for unitIdx = 1:length(good_unit_ids)
        disp(['pert ids: unit ' num2str(unitIdx) ' of ' num2str(length(good_unit_ids)) '.'])
        this_id = good_unit_ids(unitIdx);
        this_prepert_iti = squeeze(arbt_prepert_iti(:,:,this_id));
        this_pert_iti = squeeze(arbt_pert_iti(:,:,this_id));
        this_postpert_iti = squeeze(arbt_postpert_iti(:,:,this_id));
        
        
        for pIdx = 1:18
            mw_iti = ranksum(this_prepert_iti(pIdx,:), this_pert_iti(pIdx,:), 'method', 'exact', 'tail', 'left');
            this_erumd = mw_iti < sig_thresh; % 
            
            mw_iti = ranksum(this_prepert_iti(pIdx,1:10), this_prepert_iti(pIdx,end-9:end), 'method', 'exact');
            this_intra_er = mw_iti < sig_thresh;
            
            mw_iti = ranksum(this_postpert_iti(pIdx,:), this_pert_iti(pIdx,:), 'method', 'exact', 'tail', 'left');
            this_ltumd = mw_iti < sig_thresh;
            
            mw_iti = ranksum(this_postpert_iti(pIdx,1:10), this_postpert_iti(pIdx,end-9:end), 'method', 'exact');
            this_intra_lt = mw_iti < sig_thresh;
            
            mw_iti = ranksum(this_prepert_iti(pIdx,:), this_postpert_iti(pIdx,:), 'method', 'exact', 'tail', 'left');
            this_erult = mw_iti < sig_thresh;
            
            mw_iti = ranksum(this_pert_iti(pIdx,:), this_postpert_iti(pIdx,:), 'method', 'exact', 'tail', 'left');
            this_mdult = mw_iti < sig_thresh;
            
            mw_iti = ranksum(this_pert_iti(pIdx,1:10), this_pert_iti(pIdx,end-9:end), 'method', 'exact');
            this_intra_md = mw_iti < sig_thresh;
            
            
            this_resp_patt = [this_prepert_iti(pIdx,:) this_pert_iti(pIdx,:) ...
                this_postpert_iti(pIdx,:)];
            
            this_resp_patt = this_resp_patt ./ max(abs(this_resp_patt));
            trp_mean_norm = this_resp_patt ./ mean(this_prepert_iti(pIdx,:));
            
            all_ids = all_ids + 1;
    
            if (this_mdult && ~this_intra_md && this_erult && ~this_intra_er)
                big_fli_mat = [big_fli_mat; this_resp_patt];
                big_fli_ids = [big_fli_ids; all_ids];
            elseif (this_erumd && ~this_intra_er && this_ltumd && ~this_intra_lt)
                big_rdi_mat = [big_rdi_mat; this_resp_patt];
                big_rdi_ids = [big_rdi_ids; all_ids];
            elseif (~this_intra_er && ~this_intra_md && ~this_intra_lt)
                big_std_mat = [big_std_mat; this_resp_patt];
                big_std_ids = [big_std_ids; all_ids];
            end

            big_resp_mat = [big_resp_mat; this_resp_patt];
            total_ids = [total_ids; all_ids];
            
            brm_mean_norm = [brm_mean_norm; trp_mean_norm];

        end
    end

end

save(pipeline_output_save_target)