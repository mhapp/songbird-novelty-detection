% run from C:\Users\Fee/ Lab\Documents\MATLAB\NCM/ Project\useful/
% scripts\custom\

% for speed, it is important to have raw data local when processing

addpath(genpath('neuropix_pipeline'))

close all;
clear;
bird_id = '9074';
rec_date = '210811';
stim_date = '210810';
dataset_name = '9074_210811_LH_NCM_g0';

bin_name = [dataset_name '_t0.nidq.bin'];
path = ['D:\\' dataset_name];

DemoReadSGLXData(bin_name, path, rec_date)

myKsDir = [path '\' dataset_name '_imec0'];

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

%%
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

load(['workspace_' stim_date '.mat'], 'stim_timeline', 'timeline', 'syll_phase')


save_info = [path '/timingData_' rec_date '.mat'];

%theoretically, this varies... (and it got us into trouble)...
minThresh = 0.0149; % should be fine, but change if there are issues!
if strcmp(bin_name, '9074_210811_LH_NCM_g0_t0.nidq.bin')
    purge_data = [1.401e8 1.5e8+1.823e6];
    micData(purge_data(1):purge_data(2)) = 0;
elseif strcmp(bin_name, '9063_210812_LH_NCM_g0_t0.nidq.bin')
    purge_data = [1 787800];
    pd2 = [14.1e7 15e7];
    micData(purge_data(1):purge_data(2)) = 0;
    micData(pd2(1):pd2(2)) = 0;
end

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

% at end of motif period, there's a syllable period. this always starts at
% the same time, and we just need to check when.
    
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

% nb we also want syllable timing data from end of experiment

end_of_motifs = mStarts_samples_2(end) + 10000000;
this_syll_start_win = [end_of_motifs end_of_motifs+5000000];
syll_phase_starts = zeros(length(syll_phase),1);
for sIdx = 1:length(syll_phase)
    if sIdx > 1
        prev_s_start = syll_phase_starts(sIdx-1); 
        prev_s_end = prev_s_start + (.5*fs); %just before earliest conceivable start of next syll
        
        end_check = prev_s_end + (2*fs); % if we hit here, we've gone too far
    
        fr_samp = [prev_s_end end_check]; % also go forward just in case undershot (can be arbitrarily forward bc only first will be found)
    else
        fr_samp = this_syll_start_win;
    end
        
    next_s_start = find(micData(fr_samp(1):fr_samp(2)) > minThresh,1);
    
    syll_phase_starts(sIdx) = fr_samp(1) + next_s_start;
end
    
    

sStarts_samples = syll_phase_starts;
sStarts = sStarts_samples / fs;

% good for checking!
for n = 1:length(sStarts)
    figure(1); 
    if n < length(sStarts)
        plot(micData((sStarts_samples(n)):sStarts_samples(n+1)))
    else
        plot(micData((sStarts_samples(n)):sStarts_samples(n)+(.5*fs)))
    end
    title(num2str(n))
    pause(0.05);%(0.05);
end

% ALRIGHT, mStarts_samples_2 and mStarts_2 look good! - 
% maybe now we could do the other analyses

syll_phase_timingData = cell(length(sStarts), 2);
syll_phase_timingData(:,1) = syll_phase;
syll_phase_timingData(:,2) = mat2cell(sStarts, ones(1,length(sStarts)),1);
    

save(save_info, 'motifTimingData', 'syll_phase_timingData');

disp('ALL DONE WITH TIMING DATA!')

%% load checkpoint
load([path '\timingData_' rec_date '.mat']);
load([path '\syncData_' rec_date '.mat']);
load(['workspace_' stim_date '.mat'], 'timeline', 'exp_info');

fs = 40000; % update fs to that of microphone

% oh no - we'll need syllable lengths etc for all syllables. let's get this
% from the stimuli script

[micData, ~] = audioread([path '\micData_' rec_date '.wav']);

% song 6-10-21
if strcmp(stim_date, '210810')
    syllable_lengths = [4359
        1516
        4248
        3222
        5505
        3839
        3742
        3759
        9613
        8268
        7274
        2261
        2459
        4645
        5909
        6976
        4113
        11403
        10760
        3947
        5504
        7099] / 44100;
end

% matrix of offsets for each motif
syll_offsets = zeros(8, 5);

syll_offsets(1, 1:5) = [0, syllable_lengths(1)+.02, ...
    sum(syllable_lengths([1 6]))+.04, sum(syllable_lengths([1 6 11]))+.06,...
    sum(syllable_lengths([1 6 11 16]))+.08];

syll_offsets(2, 1:5) = [0, syllable_lengths(2)+.02, ...
    sum(syllable_lengths([2 7]))+.04, sum(syllable_lengths([2 7 12]))+.06,...
    0];

syll_offsets(3, 1:5) = [0, syllable_lengths(4)+.02, ...
    sum(syllable_lengths([3 8]))+.04, sum(syllable_lengths([3 8 13]))+.06,...
    sum(syllable_lengths([3 8 13 18]))+.08];

syll_offsets(4, 1:5) = [0, syllable_lengths(4)+.02, ...
    sum(syllable_lengths([4 9]))+.04, sum(syllable_lengths([4 9 14]))+.06,...
    0];

syll_offsets(5, 1:5) = [0, syllable_lengths(1)+.02, ...
    sum(syllable_lengths([1 9]))+.04, sum(syllable_lengths([1 9 11]))+.06,...
    0];

syll_offsets(6, 1:5) = [0, syllable_lengths(5)+.02, ...
    sum(syllable_lengths([5 6]))+.04, sum(syllable_lengths([5 6 11]))+.06,...
    sum(syllable_lengths([5 6 11 16]))+.08];

syll_offsets(7, 1:5) = [0, syllable_lengths(4)+.02, ...
    sum(syllable_lengths([4 10]))+.04, sum(syllable_lengths([4 10 15]))+.06,...
    0];

syll_offsets(8, 1:5) = [0, syllable_lengths(17)+.02, ...
    sum(syllable_lengths([17 12]))+.04, sum(syllable_lengths([17 12 7]))+.06,...
    0];


motif_length = syll_offsets(:,4) + syllable_lengths([21; 17; 22; 19; 16; 21; 20; 2]);

motifInfo = cell(length(motifTimingData(:,1)),3);

% how many sylls? it'll always be the same - should be 5825...
syllableInfo = cell(5825,4); 
this_syll_idx = 1;
for mIdx = 1:length(motifTimingData(:,1))
    motifInfo(mIdx,1) = motifTimingData(mIdx,2); %start time
    motifInfo(mIdx,2) = motifTimingData(mIdx,1); %start 

    mType = motifTimingData{mIdx,1};
    
    
    if strcmp(stim_date, '210810')
        sTypes = {...
            'a1','a2','a3','a4','a5',...
            'b1','b2','b3','b4','b5',...
            'c1','c2','c3','c4','c5',...
            'd1','d2','d3','d4','d5',...
            'e1','e3'};
    end
    
    
    % before this, we need to figure out which motif it is
    switch mType
        case 'std_A' % idk what the labels are
            mType_num = 1;
            this_syll_ids = [1 6 11 16 21];
            max_sub = 5;
        case 'std_B'
            mType_num = 2;
            this_syll_ids = [2 7 12 17];
            max_sub = 4;
        case 'std_C'
            mType_num = 3;
            this_syll_ids = [3 8 13 18 22];
            max_sub = 5;
        case 'std_D'
            mType_num = 4;
            this_syll_ids = [4 9 14 19];
            max_sub = 4;
        case 'std_E'
            mType_num = 5;
            this_syll_ids = [1 9 11 16];
            max_sub = 4;
        case 'std_F'
            mType_num = 6;
            this_syll_ids = [5 6 11 16 21];
            max_sub = 5;
        case 'std_G'
            mType_num = 7;
            this_syll_ids = [4 10 15 20];
            max_sub = 4;
        case 'std_H'
            mType_num = 8;
            this_syll_ids = [17 12 7 2];
            max_sub = 4;
    end
    motifInfo(mIdx,3) = {cell2mat(motifInfo(mIdx,1))+motif_length(mType_num)};
    these_sylls = sTypes(this_syll_ids);

    for sub_idx = 1:max_sub
        syllableInfo(this_syll_idx,1) = ...
            {cell2mat(motifTimingData(mIdx,2)) + syll_offsets(mType_num,sub_idx)};
        syllableInfo{this_syll_idx,2} = mType;
        syllableInfo{this_syll_idx,3} = sTypes{this_syll_ids(sub_idx)}; % this is wrong!
        syllableInfo(this_syll_idx,4) = ...
            {cell2mat(syllableInfo(this_syll_idx,1))...
            + syllable_lengths(this_syll_ids(sub_idx))};
        
        this_syll_idx = this_syll_idx + 1;
    end
end

for sp_idx = 1:length(syll_phase)
    real_id = this_syll_idx - 1 + sp_idx;
    syllableInfo(real_id,1) = {cell2mat(syll_phase_timingData(sp_idx,2))};
    syllableInfo(real_id,2) = {'NA'};
    syllableInfo(real_id,3) = syll_phase_timingData(sp_idx,1);
    
    this_syll_id = find(strcmp(syll_phase_timingData(sp_idx,1), sTypes));
    this_len = syllable_lengths(this_syll_id);
    
    syllableInfo(real_id,4) = {cell2mat(syllableInfo(real_id,1)) + this_len};
end

save('motif_checkpoint.mat')
disp('Ready for action!')

%% unit quality screen - sometimes ksDriftmap runs out of memory??? no idea what to do in that case...
% i could probably go all the way back to kilosort and disable some
% channels? ksdriftmap memory requirement depends on number of spikes

% what do i do when this happens? reload from motif_checkpoint and keep
% trying... if issues aren't solved, try different dataset

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
% this changes because we may now have a different response mat for each of
% the 8 motif types. i think that's the idea, and then we can use eval
% statements to make comparisons across these tensors

% although do we even want it? or do we want motif response mats (standard
% motifs, 4xtrialsxNeurons and alternate motfs mats, 4xtrialsxneurons).
% that sorts us for the basic stuff, and then we can have one for syllables
% x trials x neurons, along with a matrix that tells where transition
% points are for each syllable.

goodUnits = sp.cids(find(sp.cgs==2));
good_unit_ids = 1:length(goodUnits);

resp_win = [0 0.2]; % resp win extends 200ms past syllable end. NOT NECESSARILY THE RIGHT THING TO DO

% we could do better here, but at the moment we're using same matrix for
% iti during syllable phase
iti_win_syll = [-0.250 0]; % should always be at least .5 between motifs, but we want to allow up to 250ms for motif response i think...

iti_win_bouts = [-7 -2]; % at least 8 seconds between bouts - check 5 sec window in middle

% std motifs from motifs 1-700 (with first repeats 1-500, 125 each)
% alt motifs from motifs 701-1200 (125 each, no repeats)

iti_starts_syll = [syllableInfo(end-549:end,1)];
bout_starts = [motifInfo(1:5:end,1)];

% should always be the same for stim 6/10
if strcmp(stim_date, '210810')
    %idk what matters yet
end

% initialize matrix
all_motif_resp = zeros(length(motifInfo(:,1)), length(good_unit_ids));
all_syll_resp = zeros(length(syll_phase(:,1)), length(good_unit_ids));
all_bout_resp = zeros(length(bout_starts), length(good_unit_ids));

iti_resp = zeros(length(bout_starts), length(good_unit_ids));
iti_resp_syll = zeros(length(iti_starts_syll), length(good_unit_ids));


for unitIdx = 1:length(good_unit_ids)
    disp(['resp matrices for unit ' num2str(unitIdx) ' of ' num2str(length(good_unit_ids)) '.'])
    this_id = good_unit_ids(unitIdx);
    theseSpikes = sp.st(sp.clu==goodUnits(this_id));
    
    % create massive response matrix! 
    %(syllable archetype x motif num x neuron)
    for motifIdx = 1:size(motifInfo,1)
        this_win = cell2mat(motifInfo(motifIdx, [1 3])) + resp_win;
        
        motif_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(motif_spikes) / (this_win(2) - this_win(1));
        
        all_motif_resp(motifIdx, unitIdx) = resp_fr;
    end
    
    % for each inter-bout period (when motifs are separated by at least
    % 750ms), find the firing rate during that period (within iti_win)
    for syllIdx = 1:length(syllableInfo(:,1))
        this_win = cell2mat(syllableInfo(syllIdx, [1 4])) + resp_win;
        
        syll_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(syll_spikes) / (this_win(2) - this_win(1));
        
        all_syll_resp(syllIdx, unitIdx) = resp_fr;
    end
    
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for boutIdx = 1:length(bout_starts)
        this_win = cell2mat(bout_starts(boutIdx)) + resp_win;
        bout_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(bout_spikes) / (this_win(2) - this_win(1));
        
        all_bout_resp(boutIdx, unitIdx) = resp_fr;
    end
    

    for itiIdx = 1:length(bout_starts)
        this_win = cell2mat(bout_starts(itiIdx)) + iti_win_bouts;
        iti_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(iti_spikes) / (this_win(2) - this_win(1));
        
        iti_resp(itiIdx, unitIdx) = resp_fr;
    end
    
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for itiIdx = 1:length(iti_starts_syll)
        this_win = cell2mat(iti_starts_syll(itiIdx)) + iti_win_syll;
        iti_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(iti_spikes) / (this_win(2) - this_win(1));
        
        iti_resp_syll(itiIdx, unitIdx) = resp_fr;
    end
end



%% repeat resp mat calculations for multi-units

goodUnits_mu = sp.cids(find(sp.cgs==1));
good_unit_ids_mu = 1:length(goodUnits_mu);


% initialize matrix
all_motif_resp_mu = zeros(length(motifInfo(:,1)), length(good_unit_ids_mu));
all_syll_resp_mu = zeros(length(syll_phase(:,1)), length(good_unit_ids_mu));
all_bout_resp_mu = zeros(length(bout_starts), length(good_unit_ids_mu));

iti_resp_mu = zeros(length(bout_starts), length(good_unit_ids_mu));
iti_resp_syll_mu = zeros(length(iti_starts_syll), length(good_unit_ids_mu));


for unitIdx = 1:length(good_unit_ids_mu)
    disp(['resp matrices for unit ' num2str(unitIdx) ' of ' num2str(length(good_unit_ids_mu)) '.'])
    this_id = good_unit_ids_mu(unitIdx);
    theseSpikes = sp.st(sp.clu==goodUnits_mu(this_id));
    
    % create massive response matrix! 
    %(syllable archetype x motif num x neuron)
    for motifIdx = 1:size(motifInfo,1)
        this_win = cell2mat(motifInfo(motifIdx, [1 3])) + resp_win;
        
        motif_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(motif_spikes) / (this_win(2) - this_win(1));
        
        all_motif_resp_mu(motifIdx, unitIdx) = resp_fr;
    end
    
    % for each inter-bout period (when motifs are separated by at least
    % 750ms), find the firing rate during that period (within iti_win)
    for syllIdx = 1:length(syllableInfo(:,1))
        this_win = cell2mat(syllableInfo(syllIdx, [1 4])) + resp_win;
        
        syll_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(syll_spikes) / (this_win(2) - this_win(1));
        
        all_syll_resp_mu(syllIdx, unitIdx) = resp_fr;
    end
    
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for boutIdx = 1:length(bout_starts)
        this_win = cell2mat(bout_starts(boutIdx)) + resp_win;
        bout_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(bout_spikes) / (this_win(2) - this_win(1));
        
        all_bout_resp_mu(boutIdx, unitIdx) = resp_fr;
    end
    
    
    for itiIdx = 1:length(bout_starts)
        this_win = cell2mat(bout_starts(itiIdx)) + iti_win_bouts;
        iti_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(iti_spikes) / (this_win(2) - this_win(1));
        
        iti_resp_mu(itiIdx, unitIdx) = resp_fr;
    end
    
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for itiIdx = 1:length(iti_starts_syll)
        this_win = cell2mat(iti_starts_syll(itiIdx)) + iti_win_syll;
        iti_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(iti_spikes) / (this_win(2) - this_win(1));
        
        iti_resp_syll_mu(itiIdx, unitIdx) = resp_fr;
    end
end

%% eliminate low fr (do this for MU too)
fr_thresh = 1/60; % 1 spike per minute

% exclude low fr - dangerous because neuron may ONLY spike during pert
low_fr = [];
low_fr_mu = [];

goodUnits = sp.cids(find(sp.cgs==2)); % 1 for MUs, 2 for isolated
for unitIdx = 1:length(goodUnits)
    theseSpikes = sp.st(sp.clu==goodUnits(unitIdx));
    this_fr = length(theseSpikes) / (length(micData)/fs);
    if this_fr < fr_thresh
        low_fr = [low_fr; unitIdx];
    end
end

disp(['found ' num2str(length(low_fr)) ' low fr single units.'])


for unitIdx = 1:length(goodUnits_mu)
    theseSpikes = sp.st(sp.clu==goodUnits_mu(unitIdx));
    this_fr = length(theseSpikes) / (length(micData)/fs);
    if this_fr < fr_thresh
        low_fr_mu = [low_fr_mu; unitIdx];
    end
end

disp(['found ' num2str(length(low_fr_mu)) ' low fr multiunits.'])

%% eliminate units without auditory responses, also for multiunits
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

mean_iti_mu = mean(iti_resp_mu,1);
sd_iti_mu = std(iti_resp_mu,1); 

iti_floor = mean_iti - sd_iti; % bonferroni corrected bounds
iti_ceil = mean_iti + sd_iti;

iti_floor_mu = mean_iti_mu - sd_iti_mu; % bonferroni corrected bounds
iti_ceil_mu = mean_iti_mu + sd_iti_mu;

% just looking for syllable type that drives greatest response
max_resp = max(all_motif_resp, [], 1);
min_resp = min(all_motif_resp, [], 1);

max_resp_mu = max(all_motif_resp_mu, [], 1);
min_resp_mu = min(all_motif_resp_mu, [], 1);

non_aud = find(iti_floor < min_resp & iti_ceil > max_resp);

non_aud_mu = find(iti_floor_mu < min_resp_mu & iti_ceil_mu > max_resp_mu);

disp(['found  ' num2str(length(non_aud)) ' single non-aud'])

disp(['found  ' num2str(length(non_aud_mu)) ' multi non-aud'])

%% eliminate units with significant, long-lasting fr changes, SINGLE AND MULTI
drift_thresh= 10;

high_fr_drift = [];
high_fr_drift_mu = [];

win_size = 100;
step_size = 5;
moving_avg = [];
moving_avg_mu = [];

for this_start = 1:step_size:size(iti_resp,1)-(win_size-1)
    this_end = this_start + (win_size-1);
    
    this_row = mean(iti_resp(this_start:this_end,:),1);
    moving_avg = [moving_avg; this_row];
end


for this_start = 1:step_size:size(iti_resp_mu,1)-(win_size-1)
    this_end = this_start + (win_size-1);
    
    this_row = mean(iti_resp_mu(this_start:this_end,:),1);
    moving_avg_mu = [moving_avg_mu; this_row];
end

mv_avg_denom = moving_avg;
mv_avg_denom(mv_avg_denom == 0) = 0.1;

mv_avg_denom_mu = moving_avg_mu;
mv_avg_denom_mu(mv_avg_denom_mu == 0) = 0.1;

for this_col = 1:size(moving_avg,2)
    if ((max(moving_avg(:,this_col))/min(moving_avg(:,this_col)))/mean(moving_avg(:,this_col))) > drift_thresh
        high_fr_drift = [high_fr_drift; this_col];
    end
end

for this_col = 1:size(moving_avg_mu,2)
    if ((max(moving_avg_mu(:,this_col))/min(moving_avg_mu(:,this_col)))/mean(moving_avg_mu(:,this_col))) > drift_thresh
        high_fr_drift_mu = [high_fr_drift_mu; this_col];
    end
end

disp(['found ' num2str(length(high_fr_drift)) ' su with fr drift'])
disp(['found ' num2str(length(high_fr_drift_mu)) ' mu with fr drift'])
    
% 1. use motif_resp
% 2. take moving average of that (10 trials, step of 1 trial)
% 3. fr_drift = diff of that moving average
% 4. eliminate neurons where fr_drift exceeds 5 (maybe 10) hz (experiment)
% 5. update all_trial_resp, iti_resp and good_unit_ids accordingly

%% eliminate units that drift significantly in space, SINGLE ONLY
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




depth_traces = zeros(length(goodUnits_mu),length(win_starts));
density_traces = depth_traces;
ctr_samples = (sd_win_size/2):sd_step_size:(tape_length-(sd_win_size/2));

depth_means = zeros(length(goodUnits_mu),1);
depth_sds = depth_means;

for unitIdx = 1:length(goodUnits_mu)
    theseDepths = spikeDepths(sp.clu==goodUnits_mu(unitIdx));
    theseSpikes = spikeTimes(sp.clu==goodUnits_mu(unitIdx));
    
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
high_spatial_drift_mu = find(total_drift>=spatial_drift_thresh);

disp(['found ' num2str(length(high_spatial_drift)) ' unstable single units'])
disp(['found ' num2str(length(high_spatial_drift_mu)) ' unstable multiunits'])

% 1. signal_loc is average depth of unit across 10 sec window with 2 sec
% steps
% 2. eliminate neurons where range of signal_loc exceeds 30u (experiment)
% 3. update all_trial_resp, iti_resp and good_unit_ids accordingly

%% CHECK IN - HOW MANY NEURONS REMAIN? ARE THEY REASONABLE?

bad_units = union(union(low_fr, non_aud), union(high_fr_drift, high_spatial_drift));
good_unit_ids = setdiff(1:length(goodUnits), bad_units);

bad_units_mu = union(union(low_fr_mu, non_aud_mu), union(high_fr_drift_mu, high_spatial_drift_mu));
good_unit_ids_mu = setdiff(1:length(goodUnits_mu), bad_units_mu);

disp(['there are ' num2str(length(good_unit_ids)) ' good single units'])
disp(['there are ' num2str(length(good_unit_ids_mu)) ' good multiunits'])

