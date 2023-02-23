% run from C:\Users\Fee/ Lab\Documents\MATLAB\NCM/ Project\useful/
% scripts\custom\

addpath(genpath('neuropix_pipeline'))

close all;
clear;
bird_id = '9049';
rec_date = '210724';
stim_date = '210610';
dataset_name = '9049_210724_LH_NCM_g0';

bin_name = [dataset_name '_t0.nidq.bin'];
path = ['D:\\' dataset_name];

DemoReadSGLXData(bin_name, path, rec_date)

%% we don't have imec if we didn't cluster...
myKsDir = [path '\' dataset_name '_imec0'];

% [spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
% disp('done with dir')

thisBinName = [dataset_name '_t0.nidq.bin'];
thisPath = ['D:\\' dataset_name];

[sp, b, Audio_eventTimes, IMEC_eventTimes] = correct_IMEC_Spiketimes(thisPath, thisBinName);

save([thisPath '\syncData_' rec_date '.mat'], 'sp', 'b', 'Audio_eventTimes', 'IMEC_eventTimes')

disp('done with sync!')
save('sync_checkpoint.mat')

%%
spwf = sp;
cids = sp.cids(find(sp.cgs==2));

best_wfs = zeros(length(cids),82);

myKsDir = [path '\' dataset_name '_imec0'];
for neuron = 1:length(cids)
    gwfparams.dataDir = myKsDir;    % KiloSort/Phy output folder
    apD = dir(fullfile(myKsDir, '*ap*.bin')); % AP band file from spikeGLX specifically
    gwfparams.fileName = apD(1).name;         % .dat file containing the raw
    gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
    gwfparams.nCh = 385;                      % Number of channels that were streamed to disk in .dat file
    gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf = 100;                    % Number of waveforms per unit to pull out
    gwfparams.spikeTimes = ceil(spwf.st(spwf.clu==cids(neuron))*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = spwf.clu(spwf.clu==cids(neuron));
    
    wf = getWaveForms(gwfparams);
    
    theseWFs = squeeze(wf.waveFormsMean);
    wfRange = range(theseWFs,2);
    [~, maxIdx] = max(wfRange);
    
    best_wf = theseWFs(maxIdx,:);
    
    best_wfs(neuron,:) = best_wf;
end

save([thisPath '\wfData_' rec_date '.mat'], 'best_wfs')

disp('done with waveforms')

%% timeline info
% motif times start here
[micData, fs] = audioread([path '\micData_' rec_date '.wav']);

load(['workspace_' stim_date '.mat'], 'stim_timeline', 'timeline')


save_info = [path '/timingData_' rec_date '.mat'];

%theoretically, this varies... (and it got us into trouble)...
minThresh = 0.015;

timeOffset = find(micData>minThresh,1);

motifStarts = stim_timeline.motifs{2,1};
motifStarts = motifStarts + (timeOffset/fs);

x_MotifStarts = motifStarts;

motifStarts_samples = round(motifStarts*fs);
motif_lengths = diff(motifStarts_samples);

mStarts_2 = zeros(1,length(motifStarts));
mStarts_2(1) = motifStarts_samples(1);

precess = [];

for mIdx = 2:length(motifStarts)
    prev_m_start = mStarts_2(mIdx-1); 
    prev_m_end = prev_m_start + motif_lengths(mIdx-1);
    windowTime1 = 0.02;
    windowTime2 = 0.02;
    
    end_check = prev_m_end - round(windowTime1*fs); % try 10ms step back - equal to shortest ISI (in case overshot)
    fr_samp = [end_check prev_m_end+round(windowTime2*fs)]; % also go forward just in case undershot (can be arbitrarily forward bc only first will be found)
    next_m_start = find(micData(fr_samp(1):fr_samp(2)) > minThresh,1);
%     figure(3); plot(micData(fr_samp(1):fr_samp(2))); pause;
    mStarts_2(mIdx) = end_check + next_m_start;
    
    precess = [precess; prev_m_end-mStarts_2(mIdx)];
end
    
mStarts_samples_2 = mStarts_2;
mStarts_2 = mStarts_samples_2 / fs;

% good for checking!
for n = 1:length(mStarts_2)
    figure(1); 
    if n < length(mStarts_2)
        plot(micData((mStarts_samples_2(n)):mStarts_samples_2(n+1)))
    else
        plot(micData((mStarts_samples_2(n)):mStarts_samples_2(n)+fs))
    end
    title(num2str(n))
    pause(0.05);%(0.05);
end

% ALRIGHT, mStarts_samples_2 and mStarts_2 look good! - 
% maybe now we could do the other analyses

timingData = cell(length(mStarts_2), 2);
timingData(:,1) = stim_timeline.motifs{1};
timingData(:,2) = mat2cell(mStarts_2', ones(1,length(mStarts_2)),1);

motifTimingData = timingData; % this is what needs to be saved!

save(save_info, 'motifTimingData');

disp('ALL DONE WITH TIMING DATA!')

%% load checkpoint
load([path '\timingData_' rec_date '.mat']);
load([path '\syncData_' rec_date '.mat']);
load(['workspace_' stim_date '.mat'], 'timeline', 'exp_info');

fs = 40000; % update fs to that of microphone



[micData, ~] = audioread([path '\micData_' rec_date '.wav']);


% song 6-10-21
if strcmp(stim_date, '210610')
    syllable_lengths = [2397; 5777; 3735; 4824; 5491] / 40000;
elseif strcmp(stim_date, '210613')
    syllable_lengths = [3954; 3482; 6598; 6327; 4992] / 40000;
end

syllable_offsets = [0; syllable_lengths(1)+.02; sum(syllable_lengths(1:2))+.04; ...
    sum(syllable_lengths(1:3))+.06; sum(syllable_lengths(1:4))+.08];

motif_length = sum(syllable_lengths)+(4*0.02);


motifInfo = cell(length(motifTimingData(:,1)),3);
syllableInfo = cell(5*length(motifTimingData(:,1)),4);
for mIdx = 1:length(motifTimingData(:,1))
    motifInfo(mIdx,1) = motifTimingData(mIdx,2);
    motifInfo(mIdx,2) = motifTimingData(mIdx,1);
    motifInfo(mIdx,3) = {cell2mat(motifInfo(mIdx,1))+motif_length};
    
    mType = motifTimingData{mIdx,1};
    
    if strcmp(stim_date, '210610')
        sTypes = {'std_a', 'std_b', 'std_c', 'std_d', 'std_f'};
    elseif strcmp(stim_date, '210613')
        sTypes = {'std_a', 'std_b', 'std_c', 'std_d', 'std_h'};
    end
    
    if ~strcmp(mType, 'std_A')
        pSyll = mType(2);
        switch pSyll
            case 'b'
                sTypes{2} = mType;
            case 'c'
                sTypes{3} = mType;
            case 'd'
                sTypes{4} = mType;
            case 'f'
                sTypes{5} = mType;
            case 'h'
                sTypes{5} = mType;
        end
    end
    
    for sub_idx = 1:5
        this_syll_idx = (mIdx*5)-5+sub_idx;
        
        syllableInfo(this_syll_idx,1) = ...
            {cell2mat(motifTimingData(mIdx,2)) + syllable_offsets(sub_idx)};
        syllableInfo{this_syll_idx,2} = mType;
        syllableInfo{this_syll_idx,3} = sTypes{sub_idx};
        syllableInfo(this_syll_idx,4) = ...
            {cell2mat(syllableInfo(this_syll_idx,1))...
            + syllable_lengths(sub_idx)};
    end
end

save('motif_checkpoint.mat')
disp('Ready for action!')

%% check unit quality etc before plotting!
%% unit quality screen - sometimes ksDriftmap runs out of memory??? no idea what to do in that case...
% i could probably go all the way back to kilosort and disable some
% channels? ksdriftmap memory requirement depends on number of spikes
valid_contacts = [1 240];

[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);
disp('drift info here and ready')

valid_spikes = find(spikeSites<=240);

spikeSites = spikeSites(valid_spikes);
spikeDepths = spikeDepths(valid_spikes);
spikeAmps = spikeAmps(valid_spikes);
spikeTimes = spikeTimes(valid_spikes);

sp.st = sp.st(valid_spikes);
sp.spikeTemplates = sp.spikeTemplates(valid_spikes);
sp.clu = sp.clu(valid_spikes);

%% now, create response matrices
goodUnits = sp.cids(find(sp.cgs==2));
good_unit_ids = 1:length(goodUnits);

resp_win = [0 0.2]; % resp win extends 200ms past syllable end. NOT NECESSARILY THE RIGHT THING TO DO
iti_win = [-0.5 0]; % should always be at least .75 between motifs, but we want to allow up to 250ms for motif response i think...

motif_types = motifInfo(:,2);
inter_motif_intervals = cell2mat(motifInfo(2:end,1)) - cell2mat(motifInfo(1:end-1,3));
bout_start_motifs = [1; find(inter_motif_intervals > 0.5)+1]; % first one always new 'bout'
motif_starts = cell2mat(motifInfo(:,1));
syllable_types = syllableInfo(:,3);
unique_sylls = unique(syllable_types);

% should always be the same for stim 6/10
if strcmp(stim_date, '210610')
    pert_transition = 2376;
    pert_end = 12475;
    post_transition = 14476;
elseif strcmp(stim_date, '210613')
    pert_transition = 2376;
    pert_end = 12435;
    post_transition = 14436;
end

% initialize matrix
all_trial_resp = zeros(length(sTypes), size(motifInfo,1), length(good_unit_ids));
iti_resp = zeros(length(bout_start_motifs), length(good_unit_ids));
motif_resp = zeros(size(motifInfo,1), length(good_unit_ids));
std_motif_sylls = find(strcmp(syllableInfo(:,2), 'std_A'));

pert_types = unique_sylls(1:18);

avg_resp_by_type_expmt = zeros(length(unique_sylls), length(good_unit_ids));

all_resp_by_type_prepert = zeros(length(pert_types), 25, length(good_unit_ids)); 
all_resp_by_type_pert = zeros(length(pert_types), 25, length(good_unit_ids));
all_resp_by_type_postpert = zeros(length(pert_types), 25, length(good_unit_ids));

arbt_prepert_iti = zeros(length(pert_types),25,length(good_unit_ids));
arbt_pert_iti = arbt_prepert_iti;
arbt_postpert_iti = arbt_pert_iti;


for unitIdx = 1:length(good_unit_ids)
    disp(['resp matrices for unit ' num2str(unitIdx) ' of ' num2str(length(good_unit_ids)) '.'])
    this_id = good_unit_ids(unitIdx);
    theseSpikes = sp.st(sp.clu==goodUnits(this_id));
    
    % create massive response matrix! 
    %(syllable archetype x motif num x neuron)
    for syllIdx = 1:size(syllableInfo,1)
        this_motif_num = ceil(syllIdx/length(sTypes));
        this_syll_num = mod(syllIdx, length(sTypes));
        if (this_syll_num == 0)
            this_syll_num = 5;
        end
        this_win = cell2mat(syllableInfo(syllIdx,[1 4])) + resp_win;
        syll_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(syll_spikes) / (this_win(2) - this_win(1));
        
        all_trial_resp(this_syll_num, this_motif_num, unitIdx) = resp_fr;
    end
    
    % for each inter-bout period (when motifs are separated by at least
    % 750ms), find the firing rate during that period (within iti_win)
    for itiIdx = 1:length(bout_start_motifs)
        this_motif = bout_start_motifs(itiIdx);
        this_win = motif_starts(this_motif) + iti_win;
        iti_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(iti_spikes) / (this_win(2) - this_win(1));
        
        iti_resp(itiIdx, unitIdx) = resp_fr;
    end
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for motifIdx = 1:size(motifInfo,1)
        this_motif_win = cell2mat(motifInfo(motifIdx, [1 3])) + resp_win;
        motif_spikes = find((theseSpikes > this_motif_win(1)) & (theseSpikes < this_motif_win(2)));
        resp_fr = length(motif_spikes) / (this_motif_win(2) - this_motif_win(1));
        
        motif_resp(motifIdx, unitIdx) = resp_fr;
    end
    
    %for each syllable type, find average response to that syllable simply
    %by finding response to that syllable each time it occurs and averaging
    %the firing rates
    for typeIdx = 1:length(unique_sylls)
        this_type = unique_sylls(typeIdx);
        these_sylls = find(strcmp(syllable_types, this_type));
        
        temp_resp = [];
        
        for syllIdx = 1:length(these_sylls)
            this_win = cell2mat(syllableInfo(these_sylls(syllIdx),[1 4])) + resp_win;
            syll_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
            resp_fr = length(syll_spikes) / (this_win(2) - this_win(1));
        
            temp_resp = [temp_resp; resp_fr];
        end
        avg_resp_by_type_expmt(typeIdx, unitIdx) = mean(temp_resp);
    end
    
    
    % for each p type, for each occurrence of that ptype, find the previous
    % standard trial and normalize syllable response by that. average 
    % across these, separated by phase (prepert, pert, postpert)
    % also do the same with iti normalization
    for ptypeIdx = 1:length(pert_types)
        this_p_type = pert_types{ptypeIdx};
        % to id ITI, find motif number of pert, and then find iti
        % "bout_start_motifs" entry that has preceding iti
        
        
        these_stds = intersect(std_motif_sylls, find(strcmp(syllableInfo(:,3),['std_' this_p_type(2)])));
        these_perts = find(strcmp(syllableInfo(:,3),this_p_type));
        these_motif_nums = ceil(these_perts/length(sTypes));
        key_itis = zeros(length(these_motif_nums),1);
        for this_motif_idx = 1:length(these_motif_nums)
            this_this = these_motif_nums(this_motif_idx);
            this_list = bout_start_motifs - this_this;
            this_list(this_list>0) = -999999;
            [~,key_iti] = max(this_list);
            key_itis(this_motif_idx) = key_iti;
        end
        these_itis = iti_resp(key_itis);
        
        temp_resp = [];
        temp_resp_iti = [];
        for pIdx = 1:length(these_perts)
            this_p = these_perts(pIdx);
            this_std_cand = find(these_stds < this_p);
            this_std = this_std_cand(end);
            
            this_p_win = cell2mat(syllableInfo(this_p,[1 4])) + resp_win;
            syll_spikes = find((theseSpikes > this_p_win(1)) & (theseSpikes < this_p_win(2)));
            resp_fr_p = length(syll_spikes) / (this_p_win(2) - this_p_win(1));
            
            this_bl_win = cell2mat(syllableInfo(this_std,[1 4])) + resp_win;
            syll_spikes = find((theseSpikes > this_bl_win(1)) & (theseSpikes < this_bl_win(2)));
            resp_fr_bl = length(syll_spikes) / (this_bl_win(2) - this_bl_win(1));
            
            resp_fr_iti = these_itis(pIdx);
            
            temp_resp = [temp_resp; resp_fr_p - resp_fr_bl];
            temp_resp_iti = [temp_resp_iti; resp_fr_p - resp_fr_iti];
        end
        
        this_prepert = find(these_perts < pert_transition);
        this_pert = find(these_perts >= pert_transition & these_perts < pert_end);
        this_postpert = find(these_perts >= post_transition);
        
        all_resp_by_type_prepert(ptypeIdx, :, unitIdx) = temp_resp(this_prepert);
        all_resp_by_type_pert(ptypeIdx, :, unitIdx) = temp_resp(this_pert);
        all_resp_by_type_postpert(ptypeIdx, :, unitIdx) = temp_resp(this_postpert);
        
        arbt_prepert_iti(ptypeIdx,:,unitIdx) = temp_resp_iti(this_prepert);
        arbt_pert_iti(ptypeIdx,:,unitIdx) = temp_resp_iti(this_pert);
        arbt_postpert_iti(ptypeIdx,:,unitIdx) = temp_resp_iti(this_postpert);
    end        
end


%% eliminate low fr
fr_thresh = 1/60; % 1 spike per minute

% exclude low fr - dangerous because neuron may ONLY spike during pert
low_fr = [];

goodUnits = sp.cids(find(sp.cgs==2)); % 1 for MUs, 2 for isolated
for unitIdx = 1:length(goodUnits)
    theseSpikes = sp.st(sp.clu==goodUnits(unitIdx));
    this_fr = length(theseSpikes) / (length(micData)/fs);
    if this_fr < fr_thresh
        low_fr = [low_fr; unitIdx];
    end
end

disp(['found ' num2str(length(low_fr)) ' low fr units.'])

%% eliminate units without auditory responses
% really what we need is average response for each syllable type, then keep
% doing what we're doing (except instead of 2*sd, use 2.87*sd for
% bonferroni correction - pvalue of responsiveness is then .05/23)

% technically, to claim that a unit is definitely auditory, we'd want 2.87
% sds because of bonferroni. But really we just want to eliminate units
% that are clearly non-auditory, so we'll use a single std. dev

% we just can't make claims about nonexistence of certain units, but we
% wouldn't be able to do that anyway

mean_iti = mean(iti_resp,1);
sd_iti = std(iti_resp,1); 

iti_floor = mean_iti - sd_iti; % bonferroni corrected bounds
iti_ceil = mean_iti + sd_iti;

% just looking for syllable type that drives greatest response
max_resp = max(avg_resp_by_type_expmt, [], 1);
min_resp = min(avg_resp_by_type_expmt, [], 1);

non_aud = find(iti_floor < min_resp & iti_ceil > max_resp);

disp(['found  ' num2str(length(non_aud)) ' non-aud'])

%% eliminate units with significant, long-lasting fr changes
drift_thresh= 10;

high_fr_drift = [];

win_size = 100;
step_size = 5;
moving_avg = [];

for this_start = 1:step_size:size(iti_resp,1)-(win_size-1)
    this_end = this_start + (win_size-1);
    
    this_row = mean(iti_resp(this_start:this_end,:),1);
    moving_avg = [moving_avg; this_row];
end

mv_avg_denom = moving_avg;
mv_avg_denom(mv_avg_denom == 0) = 0.1;

for this_col = 1:size(moving_avg,2)
    if ((max(moving_avg(:,this_col))/min(moving_avg(:,this_col)))/mean(moving_avg(:,this_col))) > drift_thresh
        high_fr_drift = [high_fr_drift; this_col];
    end
end

disp(['found ' num2str(length(high_fr_drift)) ' with fr drift'])
    
% 1. use motif_resp
% 2. take moving average of that (10 trials, step of 1 trial)
% 3. fr_drift = diff of that moving average
% 4. eliminate neurons where fr_drift exceeds 5 (maybe 10) hz (experiment)
% 5. update all_trial_resp, iti_resp and good_unit_ids accordingly

%% eliminate units that drift significantly in space
sd_win_size = 10*fs; % 10 second window
sd_step_size = 2*fs; % 1 second 
spatial_drift_thresh = 30;

high_spatial_drift = [];

% deffo_bad = load('def_bad_210119b.mat');
tape_length = round(max(round(spikeTimes*fs))+(5*fs), 5, 'significant');
win_starts = 1:sd_step_size:(tape_length-sd_win_size)+1;
win_ends = win_starts + sd_win_size - 1;
depth_traces = zeros(length(goodUnits),length(win_starts));
density_traces = depth_traces;
ctr_samples = (sd_win_size/2):sd_step_size:(tape_length-(sd_win_size/2));

depth_means = zeros(length(goodUnits),1);
depth_sds = depth_means;

for unitIdx = 1:length(goodUnits)
    theseDepths = spikeDepths(sp.clu==goodUnits(unitIdx));
    theseSpikes = spikeTimes(sp.clu==goodUnits(unitIdx));
    
    depth_means(unitIdx) = mean(theseDepths);
    depth_sds(unitIdx) = std(theseDepths);
    
    for winIdx = 1:length(win_starts)
        this_win_start = win_starts(winIdx)/fs;
        this_win_end  = win_ends(winIdx)/fs;
        these_spike_ids = find(theseSpikes>this_win_start & theseSpikes<this_win_end);
        these_sds = theseDepths(these_spike_ids);
        depth_traces(unitIdx,winIdx) = mean(these_sds);
        density_traces(unitIdx,winIdx) = length(these_spike_ids)/(this_win_end-this_win_start);
    end
end

total_drift = range(depth_traces,2);
high_spatial_drift = find(total_drift>=spatial_drift_thresh);

disp(['found ' num2str(length(high_spatial_drift)) ' unstable units'])

% 1. signal_loc is average depth of unit across 10 sec window with 2 sec
% steps
% 2. eliminate neurons where range of signal_loc exceeds 30u (experiment)
% 3. update all_trial_resp, iti_resp and good_unit_ids accordingly

%% CHECK IN - HOW MANY NEURONS REMAIN? ARE THEY REASONABLE?

bad_units = union(union(low_fr, non_aud), union(high_fr_drift, high_spatial_drift));
good_unit_ids = setdiff(1:length(goodUnits), bad_units);

disp(['there are ' num2str(length(good_unit_ids)) ' good single units'])


%% 
save(['D:\\pipeline_outputs\pipeline_' dataset_name '.mat'])
disp('pipeline complete!')