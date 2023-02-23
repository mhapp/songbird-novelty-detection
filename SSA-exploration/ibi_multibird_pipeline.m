%%% ibi multibird summary

% to save time in analysis, we run this code one time after running the
% prelim pipeline on each bird. the result is large multi-bird data
% structures summarizing timing response data

the_good_birds = {'D:\ibi_expmts\pipelines\pipeline_9131_220227_RH_NCM_g0.mat',
    'D:\ibi_expmts\pipelines\pipeline_9138_220228_RH_NCM_g0.mat',
    'D:\ibi_expmts\pipelines\pipeline_9329_220224_RH_NCM_g0.mat'};

pipeline_output_save_target = 'D:\ibi_expmts\multibird_summary_ibi.mat';

total_motif_resp = [];
total_motif_resp_mu = [];
total_wfs = [];

num_goods = [];
num_goods_mu = [];

big_ibi_tensor = [];
big_ibi_frac = [];
for big_bird_num = 1:length(the_good_birds)
    clearvars -except big_bird_num the_good_birds ...
        pipeline_output_save_target total_motif_resp total_motif_resp_mu...
        total_wfs num_goods num_goods_mu big_ibi_tensor big_ibi_frac
    close all
    
    this_bird = the_good_birds{big_bird_num};
    load(this_bird)
        
    disp(['starting bird ' num2str(big_bird_num)])
    
    total_motif_resp = [total_motif_resp all_motif_resp(:,good_unit_ids)];
    total_motif_resp_mu...
        = [total_motif_resp_mu all_motif_resp_mu(:,good_unit_ids_mu)];
    total_wfs = [total_wfs; best_wfs(good_unit_ids, 1:81)];
    
    num_goods = [num_goods; length(good_unit_ids)];
    num_goods_mu = [num_goods_mu; length(good_unit_ids_mu)];
    
    
    % IBI tensors for summary plotting
    motif_ids = stim_timeline.motifs{1,:};

    ibi_ids = find(contains(motif_ids,'IBI'));
    key_motifs = [1 2751; 1251 2801];

    ibi_tags = timeline.motifs(ibi_ids);
    ibi_tags_s1 = ibi_tags(1:250);
    ibi_tags_s2 = ibi_tags(251:500);
    ibi_types = unique(ibi_tags);
    
    % songs 1 and 2
    ibi_tensor_raw_su = zeros(2,length(good_unit_ids),25,10);
    ibi_tensor_raw_mu = zeros(2,length(good_unit_ids_mu),25,10);
    ibi_tensor_frac_su = zeros(2,length(good_unit_ids),25,10);
    ibi_tensor_frac_mu = zeros(2,length(good_unit_ids_mu),25,10);
    
    ctrl_tensor_raw_su = ibi_tensor_raw_su;
    ctrl_tensor_raw_mu = ibi_tensor_raw_su;
    ctrl_tensor_frac_su = ibi_tensor_raw_su;
    ctrl_tensor_frac_mu = ibi_tensor_raw_su;

    for uType = [1,2]
        if uType == 1
            this_uType_str = 'mu';
            these_good_units = good_unit_ids_mu;
            these_motif_resp = all_motif_resp_mu;
            this_tensor_raw = ibi_tensor_raw_mu;
            this_tensor_frac = ibi_tensor_frac_mu;
            this_ctrl_raw = ctrl_tensor_raw_mu;
            this_ctrl_frac = ctrl_tensor_frac_mu;
        elseif uType == 2
            this_uType_str = 'su';
            these_good_units = good_unit_ids;
            these_motif_resp = all_motif_resp;
            this_tensor_raw = ibi_tensor_raw_su;
            this_tensor_frac = ibi_tensor_frac_su;
            this_ctrl_raw = ctrl_tensor_raw_su;
            this_ctrl_frac = ctrl_tensor_frac_su;
        end
    
        for songIdx = 1:2 % only first two songs are important
            tkm = key_motifs(songIdx,:);
            
            these_frac_inc_values = zeros(length(these_good_units), 250); % if error, numBouts is not 25
            these_raw_inc_values = zeros(length(these_good_units), 250);
            these_dev_values = zeros(length(these_good_units), 250);
            
            bout_nums = [ceil(tkm(1)/5) ceil(tkm(1)/5)+249];
            these_IBIs = diff(cell2mat(bout_starts(bout_nums(1):bout_nums(2)))); % make sure these line up!
            % first IBI follows first bout
        
        
            for unitIdx = 1:length(these_good_units)
                disp([this_uType_str ', song' num2str(songIdx) ', unit' num2str(unitIdx)])
                ibi_counts = ones(length(ibi_types),1);
                disp([this_uType_str num2str(unitIdx) ' of ' num2str(length(these_good_units)) ', Song ' num2str(songIdx)])
                this_id = these_good_units(unitIdx);
                this_resp_p1 = these_motif_resp(tkm(1):tkm(1)+1249, this_id);
                this_resp_p2 = these_motif_resp(tkm(2):tkm(2)+49, this_id);
                
                this_max_p1 = max(this_resp_p1);
                
                num_bouts = length(this_resp_p1)/5;
                
                for bIdx = 1:num_bouts
                    this_first_motif = (bIdx-1)*5 + 1;
                    if bIdx > 1
                        this_ibi = ibi_tags(((songIdx-1)*250)+bIdx-1);
                        this_ibi_idx = find(strcmp(ibi_types,this_ibi));
                        
                        ibi_tape_num = ibi_counts(this_ibi_idx);
                        
                        this_resp = this_resp_p1(this_first_motif:this_first_motif+4);
                        last_resp = this_resp_p1(this_first_motif-5:this_first_motif-1);
                        
                        raw_ib_rbd = this_resp(1) - last_resp(end);
                        frac_ib_rbd = raw_ib_rbd / this_max_p1;
                        
                        ctrl_raw_rbd = this_resp(2) - last_resp(end);
                        ctrl_frac_rbd = this_resp(2) - last_resp(end);
                        
                        this_tensor_raw(songIdx,unitIdx,ibi_tape_num,this_ibi_idx) = raw_ib_rbd;
                        this_tensor_frac(songIdx,unitIdx,ibi_tape_num,this_ibi_idx) = frac_ib_rbd;
                        
                        this_ctrl_raw(songIdx,unitIdx,ibi_tape_num,this_ibi_idx) = ctrl_raw_rbd;
                        this_ctrl_frac(songIdx,unitIdx,ibi_tape_num,this_ibi_idx) = ctrl_frac_rbd;
                        
                        ibi_counts(this_ibi_idx) = ibi_counts(this_ibi_idx) + 1;
                    end
                end
            end
        end
        
        eval(['ibi_tensor_raw_' this_uType_str '= this_tensor_raw;'])
        eval(['ibi_tensor_frac_' this_uType_str '= this_tensor_frac;'])
        
        eval(['ctrl_tensor_raw_' this_uType_str '= this_ctrl_raw;'])
        eval(['ctrl_tensor_frac_' this_uType_str '= this_ctrl_frac;'])
    end
    
    big_ibi_tensor = cat(2, big_ibi_tensor, ibi_tensor_raw_su);
    big_ibi_frac = cat(2, big_ibi_frac, ibi_tensor_frac_su);
    % final shape is (song id X unit num X rep id X ibi id)

    disp('bird done')
end

clearvars -except pipeline_output_save_target total_motif_resp...
    total_motif_resp_mu total_wfs motifInfo num_goods num_goods_mu ...
    big_ibi_tensor big_ibi_frac

save(pipeline_output_save_target)
disp('ibi data across birds complete!')