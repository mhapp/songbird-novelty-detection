%%% SSA multibird summary

% to save time in analysis, we run this code one time after running the
% prelim pipeline on each bird. the result is large multi-bird data
% structures summarizing SSA data
%
% nb, this is now very slow because we grab response vectors for each unit

the_good_birds = {'D:\ssa_expmts\pipelines\pipeline_9059_210811_LH_NCM_g0.mat',
    'D:\ssa_expmts\pipelines\pipeline_9072_210810_LH_NCM_g0_v2.mat',
    'D:\ssa_expmts\pipelines\pipeline_9074_210811_LH_NCM_g0.mat'};

pipeline_output_save_target = 'D:\ssa_expmts\multibird_summary_ssa.mat';

total_motif_resp = [];
total_motif_resp_mu = [];
total_wfs = [];

num_goods = [];
num_goods_mu = [];

big_rv1 = [];
big_rv2 = [];
big_rv3 = [];
big_rv4 = [];


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
    
    % now get response vectors for single units
    resp_vect_win = .05; %50 ms gauss window - see what others do though
    npoints = resp_vect_win*fs;
    resp_vect_filt = gausswin(npoints); %alpha 2.5
    
    bin_win = round(0.05*fs);
    bin_step = round(0.01*fs);
    
    tape_lengths = zeros(4,1);
    for songIdx = 1:4
        tape_lengths(songIdx) = length(1:bin_step:(round(motif_length(songIdx)*fs)-bin_win)+1);
    end
    
    impt_motifs = [1 125 501 550;...
        126 250 551 600;...
        251 375 601 650;...
        376 500 651 700];

    resp_vector_s1 = zeros(length(good_unit_ids), 175, tape_lengths(1));
    resp_vector_s2 = zeros(length(good_unit_ids), 175, tape_lengths(2));
    resp_vector_s3 = zeros(length(good_unit_ids), 175, tape_lengths(3));
    resp_vector_s4 = zeros(length(good_unit_ids), 175, tape_lengths(4));
    
    for songIdx = 1:4
        this_motif_length = motif_length(songIdx);
        this_win_frame = [0 this_motif_length]; % seconds
        this_impt_line = impt_motifs(songIdx,:);
        these_motif_starts_1 = motifStarts(this_impt_line(1):this_impt_line(2));
        these_motif_starts_2 = motifStarts(this_impt_line(3):this_impt_line(4));
        these_motif_starts = [these_motif_starts_1'; these_motif_starts_2'];
        eval(['this_sType_str = ''s' num2str(songIdx) ''';'])
        
        these_good_units = sp.cids(find(sp.cgs==2)); %single units
        
        these_ids = good_unit_ids;
        
        eval(['this_resp_vector = resp_vector_' this_sType_str ';'])
        
        for unitIdx = 1:length(these_ids)
            
            
            this_unit_id = these_ids(unitIdx);
            disp(['song' num2str(songIdx) ', unit' num2str(unitIdx) ' of ' num2str(length(these_ids))])
            theseSpikes = sp.st(sp.clu==these_good_units(this_unit_id));
            
            for trialIdx = 1:175
                this_start = these_motif_starts(trialIdx);
                this_response_window = this_start + this_win_frame;
                
                relevant_spikes = find(theseSpikes<this_response_window(2) & theseSpikes>this_response_window(1)); % indices
                goodSpikes = round(theseSpikes(relevant_spikes)*fs) - round(this_start*fs); % samples
                
                this_train = zeros(round(motif_length(songIdx)*fs),1);
                this_train(goodSpikes+1) = 1;
                
                smooth_train = conv(this_train, resp_vect_filt, 'same');
                
                
                bin_starts = 1:bin_step:(length(smooth_train)-bin_win+1);
                bin_ends = bin_starts+bin_win-1;
                this_tape = zeros(length(bin_starts),1);
                
                for bin_idx = 1:length(bin_starts)
                    this_window = [bin_starts(bin_idx) bin_ends(bin_idx)];
                    these_resps = sum(smooth_train(this_window(1):this_window(2)));
                    this_tape(bin_idx) = these_resps;
                end
                
                this_resp_vector(unitIdx,trialIdx,:) = this_tape; % doesn't actually change size - i was clever with evals
                
                eval(['resp_vector_' this_sType_str  '= this_resp_vector;'])
            end
        end
    end
    
    big_rv1 = cat(1,big_rv1, resp_vector_s1); 
    big_rv2 = cat(1,big_rv2, resp_vector_s2); 
    big_rv3 = cat(1,big_rv3, resp_vector_s3);
    big_rv4 = cat(1,big_rv4, resp_vector_s4); 
    disp('bird done')
end

clearvars -except pipeline_output_save_target total_motif_resp...
    total_motif_resp_mu total_wfs motifInfo num_goods num_goods_mu ...
    big_rv1 big_rv2 big_rv3 big_rv4

save(pipeline_output_save_target)
disp('ssa across birds complete!')
    