%%% timing multibird summary

% to save time in analysis, we run this code one time after running the
% prelim pipeline on each bird. the result is large multi-bird data
% structures summarizing timing response data

the_good_birds = {'D:\timing_expmts\pipelines\pipeline_9206_222201_LH_NCM_g0.mat',
    'D:\timing_expmts\pipelines\pipeline_9299_222301_LH_NCM_g0.mat'};

pipeline_output_save_target = 'D:\timing_expmts\multibird_summary_timing.mat';

total_motif_resp = [];
total_motif_resp_mu = [];
total_wfs = [];

num_goods = [];
num_goods_mu = [];

for big_bird_num = 1:length(the_good_birds)
    clearvars -except big_bird_num the_good_birds ...
        pipeline_output_save_target total_motif_resp total_motif_resp_mu...
        total_wfs num_goods num_goods_mu big_rv1 big_rv2 big_rv3 big_rv4
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
    

    disp('bird done')
end

clearvars -except pipeline_output_save_target total_motif_resp...
    total_motif_resp_mu total_wfs motifInfo num_goods num_goods_mu 

save(pipeline_output_save_target)
disp('timing data across birds complete!')