% run from C:\Users\Fee/ Lab\Documents\MATLAB\NCM/ Project\useful/
% scripts\custom\

addpath(genpath('neuropix_pipeline'))

close all;
clear;
bird_id = '9138';
rec_date = '220228';
stim_date = '220224';
dataset_name = '9138_220228_RH_NCM_g0';

bin_name = [dataset_name '_t0.nidq.bin'];
path = ['D:\\' dataset_name];

DemoReadSGLXData(bin_name, path, rec_date)

%% we don't have imec if we didn't cluster...
% not much to be done if we don't even hve memory here...

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

% for some reason, we can get waveforms here in 9062, 8/12 - just skipping
% individual problematic neurons. not too clear what consequences will be.
% could be the case that some neurons just have rates that are too high?
% would be great if it were easy to check this... let's figure it out
% doesn't seem to be the case... it doesn't necessarily stop on neurons
% that have high firing rates
%
% in fact, it's extremely unclear what is happening. we don't get issues on
% the same unit every time, and sometimes it just works... it's not some
% rolling memory issue, since sometimes, upon restarting, it stops much
% sooner than the original time


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

% this needs to be altered because IBI is included in motif timeline now


% motif times start here
[micData, fs] = audioread([path '\micData_' rec_date '.wav']);

load([path '/workspace_' stim_date '.mat'], 'stim_timeline', 'timeline')

save_info = [path '/timingData_' rec_date '.mat'];

purge_until = 0;

%theoretically, this varies... (and it got us into trouble)...
if strcmp(dataset_name, '9206_222201_LH_NCM_g0')
    minThresh = 0.0149;
elseif strcmp(dataset_name, '9131_220227_RH_NCM_g0')
    purge_until = 1.36e6;
    minThresh = 0.0125; % should be fine
elseif strcmp(dataset_name, '9138_220228_RH_NCM_g0')
    minThresh = 0.0125;
    purge_until = 2e6;
else
    minThresh = 0.0125;
end

micData(1:purge_until) = 0;

timeOffset = find(micData>minThresh,1);

motifStarts = stim_timeline.motifs{2,1};
motifStarts = motifStarts + (timeOffset/fs);

x_MotifStarts = motifStarts;

motifStarts_samples = round(motifStarts*fs);
motif_lengths = diff(motifStarts_samples);

motif_ids = stim_timeline.motifs{1,:};

ibi_ids = find(contains(motif_ids,'IBI'));
true_ids = setdiff(1:length(motif_ids), ibi_ids);

real_motif_lengths = motif_lengths(true_ids);

mStarts_2 = zeros(1,length(motifStarts)-length(ibi_ids));
mStarts_2(1) = motifStarts_samples(1);

precess = [];

for mIdx = 2:length(mStarts_2)
    prev_m_start = mStarts_2(mIdx-1); 
    prev_m_end = prev_m_start + real_motif_lengths(mIdx-1);
    
%     if mIdx >= 10
%         disp('hold on')
%     end
    
    
    this_motif_str = motif_ids{true_ids(mIdx)-1};
    if contains(this_motif_str, 'IBI')
        this_ibi_len = (str2double(this_motif_str(end-3:end))/1000)*fs;
        prev_m_end = prev_m_end + this_ibi_len;
    end
    
    % these are genuinely useful bc onset detection is tricky
    windowTime1 = 0.02;
    windowTime2 = 0.05;
    
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

%good for checking!
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
timingData(:,1) = motif_ids(true_ids);
timingData(:,2) = mat2cell(mStarts_2', ones(1,length(mStarts_2)),1);

motifTimingData = timingData; % this is what needs to be saved!


save(save_info, 'motifTimingData');

disp('ALL DONE WITH TIMING DATA!')


%% load checkpoint
load([path '\timingData_' rec_date '.mat']);
load([path '\syncData_' rec_date '.mat']);
load([path '\workspace_' stim_date '.mat'], 'timeline', 'exp_info');

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
elseif strcmp(stim_date, '222201')
     syllable_lengths = [4359
        3839
        7274
        6976
        5504
        2643
        6369
        4118
        5319
        6054
        4248
        3759
        2459
        11403
        7099
        ] /44100;
elseif strcmp(stim_date, '220224')
     syllable_lengths = [4359 % song 1
        3839
        7274
        6976
        5504
        2643 % song 2a
        6369 % 2b
        4118 %2c
        5319 %2d
        6054 %2e
        3475 % song 3a, 11
        2927 %3b
        2581 %3c
        8852 %3d
        9600 %3e
        7380 %3f
        8614 %3g
        3815 %3h
        11170 %3i, 19
        2302 % song 4a
        5487] /44100;
end

% matrix of offsets for each motif
syll_offsets = zeros(4, 9);

syll_offsets(1, 1:5) = [0, syllable_lengths(1)+.02, ...
    sum(syllable_lengths([1 2]))+.04, sum(syllable_lengths([1 2 3]))+.06,...
    sum(syllable_lengths([1 2 3 4]))+.08];

syll_offsets(2, 1:5) = [0, syllable_lengths(2)+.02, ...
    sum(syllable_lengths([6 7]))+.04, sum(syllable_lengths([6 7 8]))+.06,...
    sum(syllable_lengths([6 7 8 9]))+.08];

syll_offsets(3, 1:9) = [0, syllable_lengths(11)+.03, ...
    sum(syllable_lengths([11 12]))+.06, sum(syllable_lengths([11 12 13]))+.09,...
    sum(syllable_lengths([11 12 13 14]))+.12, sum(syllable_lengths([11:15]))+.15,...
    sum(syllable_lengths([11:16]))+.18, sum(syllable_lengths([11:17]))+.21, ...
    sum(syllable_lengths([11:18]))+.24,];

syll_offsets(4, 1:2) = [0, syllable_lengths(20)+.02];

motif_length(1) = syll_offsets(1,5) + syllable_lengths(5);
motif_length(2) = syll_offsets(2,5) + syllable_lengths(10);
motif_length(3) = syll_offsets(3,9) + syllable_lengths(19);
motif_length(4) = syll_offsets(4,2) + syllable_lengths(21);

motifInfo = cell(length(motifTimingData(:,1)),3);

% how many sylls? it'll always be the same - should be 5825...
syllableInfo = cell(14925,4); 
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
    elseif strcmp(stim_date, '222201')
        sTypes = {...
            'a1','b1','c1','d1','e1',...
            'a2','b2','c2','d2','e2',...
            'a3','b3','c3','d3','e3'};
    elseif strcmp(stim_date, '220224')
        sTypes = {...
            'a1','b1','c1','d1','e1',...
            'a2','b2','c2','d2','e2',...
            'a3','b3','c3','d3','e3','f3','g3','h3','i3',...
            'a4','b4'};
    end
    
    
    % before this, we need to figure out which motif it is
    switch mType
        case 'std_A' % idk what the labels are
            mType_num = 1;
            this_syll_ids = [1 2 3 4 5];
            max_sub = 5; %idk what max sub is - max syll id?
        case 'std_B'
            mType_num = 2;
            this_syll_ids = [6 7 8 9 10];
            max_sub = 5;
        case 'std_C'
            mType_num = 3;
            this_syll_ids = [11 12 13 14 15 16 17 18 19];
            max_sub = 9;
        case 'std_D'
            mType_num = 4;
            this_syll_ids = [20 21];
            max_sub = 2;
    end
    motifInfo(mIdx,3) = {cell2mat(motifInfo(mIdx,1))+motif_length(mType_num)};
    these_sylls = sTypes(this_syll_ids);

    for sub_idx = 1:max_sub
        syllableInfo(this_syll_idx,1) = ...
            {cell2mat(motifTimingData(mIdx,2)) + syll_offsets(mType_num,sub_idx)};
        syllableInfo{this_syll_idx,2} = mType;
        syllableInfo{this_syll_idx,3} = sTypes{this_syll_ids(sub_idx)}; % this is wrong! % ??
        syllableInfo(this_syll_idx,4) = ...
            {cell2mat(syllableInfo(this_syll_idx,1))...
            + syllable_lengths(this_syll_ids(sub_idx))};
        
        this_syll_idx = this_syll_idx + 1;
    end
end

bout_starts = cell2mat(motifInfo(1:5:end,1));

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

disp('we survived memory issues!')

%% now, create response matrices
% this changes because we may now have a different response mat for each of
% the 8 motif types. i think that's the idea, and then we can use eval
% statements to make comparisons across these tensors

% although do we even want it? or do we want motif response mats (standard
% motifs, 4xtrialsxNeurons and alternate motfs mats, 4xtrialsxneurons).
% that sorts us for the basic stuff, and then we can have one for syllables
% x trials x neurons, along with a matrix that tells where transition
% points are for each syllable.

% NB: can't really do ITI windows during first epoch

goodUnits = sp.cids(find(sp.cgs==2));
good_unit_ids = 1:length(goodUnits);

resp_win = [0 0.2]; % resp win extends 200ms past syllable end. NOT NECESSARILY THE RIGHT THING TO DO

bout_starts = [motifInfo(1:5:end,1)];
bout_lengths = [cell2mat(motifInfo(5,3)) - cell2mat(motifInfo(1,1));
    cell2mat(motifInfo(1255,3)) - cell2mat(motifInfo(1251,1));
    cell2mat(motifInfo(2505,3)) - cell2mat(motifInfo(2501,1));
    cell2mat(motifInfo(2630,3)) - cell2mat(motifInfo(2626,1))];


% prep for getting non-aud responses
pre_frame_size = mean(motif_length) * fs;
pre_expmt = [1 timeOffset];
pre_win_starts = (pre_expmt(1):round(pre_frame_size/2):pre_expmt(2)-pre_frame_size) ./ fs;
pre_win_ends = (pre_win_starts + (pre_frame_size/fs));

% initialize matrix
all_motif_resp = zeros(length(motifInfo(:,1)), length(good_unit_ids));
all_bout_resp = zeros(length(bout_starts), length(good_unit_ids));
all_pre_resp = zeros(length(pre_win_starts), length(good_unit_ids));


for unitIdx = 1:length(good_unit_ids)
    disp(['resp matrices for unit ' num2str(unitIdx) ' of ' num2str(length(good_unit_ids)) '.'])
    this_id = good_unit_ids(unitIdx);
    theseSpikes = sp.st(sp.clu==goodUnits(this_id));
    
    for pre_win_idx = 1:length(pre_win_starts)
        this_win = [pre_win_starts(pre_win_idx) pre_win_ends(pre_win_idx)];
        pre_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        pre_fr = length(pre_spikes) / (this_win(2) - this_win(1));
        all_pre_resp(pre_win_idx,unitIdx) = pre_fr;
    end
    
    % create massive response matrix! 
    %(syllable archetype x motif num x neuron)
    for motifIdx = 1:size(motifInfo,1)
        this_win = cell2mat(motifInfo(motifIdx, [1 3])) + resp_win;
        
        motif_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(motif_spikes) / (this_win(2) - this_win(1));
        
        all_motif_resp(motifIdx, unitIdx) = resp_fr;
    end
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for boutIdx = 1:length(bout_starts)
        if boutIdx < 251  % i know but it's late and i can't be bothered
            typeIdx = 1;
        elseif boutIdx < 501
            typeIdx = 2;
        elseif boutIdx < 526
            typeIdx = 3;
        elseif boutIdx < 551
            typeIdx = 4;
        elseif boutIdx < 561
            typeIdx = 1;
        elseif boutIdx < 571
            typeIdx = 2;
        elseif boutIdx < 581
            typeIdx = 3;
        else
            typeIdx = 4;
        end
        
        
        this_win = [cell2mat(bout_starts(boutIdx)) ...
            cell2mat(bout_starts(boutIdx))+bout_lengths(typeIdx)]  + resp_win;
        bout_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(bout_spikes) / (this_win(2) - this_win(1));
        
        all_bout_resp(boutIdx, unitIdx) = resp_fr;
    end
end

%% repeat resp mat calculations for multi-units

goodUnits_mu = sp.cids(find(sp.cgs==1));
good_unit_ids_mu = 1:length(goodUnits_mu);


% initialize matrix
all_motif_resp_mu = zeros(length(motifInfo(:,1)), length(good_unit_ids_mu));
all_bout_resp_mu = zeros(length(bout_starts), length(good_unit_ids_mu));
all_pre_resp_mu = zeros(length(pre_win_starts), length(good_unit_ids_mu));


for unitIdx = 1:length(good_unit_ids_mu)
    disp(['resp matrices for unit ' num2str(unitIdx) ' of ' num2str(length(good_unit_ids_mu)) '.'])
    this_id = good_unit_ids_mu(unitIdx);
    theseSpikes = sp.st(sp.clu==goodUnits_mu(this_id));
    
    for pre_win_idx = 1:length(pre_win_starts)
        this_win = [pre_win_starts(pre_win_idx) pre_win_ends(pre_win_idx)];
        pre_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        pre_fr = length(pre_spikes) / (this_win(2) - this_win(1));
        all_pre_resp_mu(pre_win_idx,unitIdx) = pre_fr;
    end
    
    % create massive response matrix! 
    %(syllable archetype x motif num x neuron)
    for motifIdx = 1:size(motifInfo,1)
        this_win = cell2mat(motifInfo(motifIdx, [1 3])) + resp_win;
        
        motif_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(motif_spikes) / (this_win(2) - this_win(1));
        
        all_motif_resp_mu(motifIdx, unitIdx) = resp_fr;
    end
    
    
    
    % for each motif, find the response to that motif (again allowing a
    % resp_win response window)
    for boutIdx = 1:length(bout_starts)
        if boutIdx < 251  % i know but it's late and i can't be bothered
            typeIdx = 1;
        elseif boutIdx < 501
            typeIdx = 2;
        elseif boutIdx < 526
            typeIdx = 3;
        elseif boutIdx < 551
            typeIdx = 4;
        elseif boutIdx < 561
            typeIdx = 1;
        elseif boutIdx < 571
            typeIdx = 2;
        elseif boutIdx < 581
            typeIdx = 3;
        else
            typeIdx = 4;
        end
        
        this_win = [cell2mat(bout_starts(boutIdx)) ...
            cell2mat(bout_starts(boutIdx))+bout_lengths(typeIdx)]  + resp_win;
        bout_spikes = find((theseSpikes > this_win(1)) & (theseSpikes < this_win(2)));
        resp_fr = length(bout_spikes) / (this_win(2) - this_win(1));
        
        all_bout_resp_mu(boutIdx, unitIdx) = resp_fr;
    end
    
end

%% eliminate low fr (do this for MU too)
fr_thresh = 1/10; % 1 spike per 10 seconds, overall
phase_thresh = 20*10; % 200 total spikes in phase 2

% exclude low fr - dangerous because neuron may ONLY spike during pert
low_fr = [];
low_fr_mu = [];

p2_times = [motifInfo{2751,1} motifInfo{2950,3}];

goodUnits = sp.cids(find(sp.cgs==2)); % 1 for MUs, 2 for isolated
for unitIdx = 1:length(goodUnits)
    theseSpikes = sp.st(sp.clu==goodUnits(unitIdx));
    
    p2_spikes = find(theseSpikes > p2_times(1) & theseSpikes < p2_times(2));
    
    this_fr = length(theseSpikes) / (length(micData)/fs);
    if this_fr < fr_thresh
        low_fr = [low_fr; unitIdx];
    elseif p2_spikes < phase_thresh
        low_fr = [low_fr; unitIdx];
    end
end

disp(['found ' num2str(length(low_fr)) ' low fr single units.'])


for unitIdx = 1:length(goodUnits_mu)
    theseSpikes = sp.st(sp.clu==goodUnits_mu(unitIdx));
    p2_spikes = find(theseSpikes > p2_times(1) & theseSpikes < p2_times(2));
    this_fr = length(theseSpikes) / (length(micData)/fs);
    if this_fr < fr_thresh
        low_fr_mu = [low_fr_mu; unitIdx];
    elseif p2_spikes < phase_thresh
        low_fr_mu = [low_fr_mu; unitIdx];
    end
end

disp(['found ' num2str(length(low_fr_mu)) ' low fr multiunits.'])

%% eliminate units without auditory responses, also for multiunits
% do this with pre-experiment time

mean_iti = mean(all_pre_resp,1);
sd_iti = std(all_pre_resp,1); 

mean_iti_mu = mean(all_pre_resp_mu,1);
sd_iti_mu = std(all_pre_resp_mu,1); 

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

% for each song type, check phase 1 and 2 and confirm similar
% probably way too gross - not eliminating anything

all_bout_means = zeros(2,4,size(all_bout_resp,2));
all_bout_sems = zeros(2,4,size(all_bout_resp,2));

all_bout_means_mu = zeros(2,4,size(all_bout_resp_mu,2));
all_bout_sems_mu = zeros(2,4,size(all_bout_resp_mu,2));

first_bouts = [1 251 501 526;
    551 561 571 581];
last_bouts = [250 500 525 550;
    560 570 580 590];

for uType = 1:2
    if uType == 1
        this_summary = all_bout_resp_mu;
        this_means = all_bout_means_mu;
        this_sems = all_bout_sems_mu;
    else
        this_summary = all_bout_resp;
        this_means = all_bout_means;
        this_sems = all_bout_sems;
    end
    for pIdx = 1:2
        for songIdx = 1:4
            this_first = first_bouts(pIdx,songIdx);
            this_last = last_bouts(pIdx,songIdx);
            
            this_slice = this_summary(this_first:this_last,:);
            this_means(pIdx,songIdx,:) = mean(this_summary, 1);
            this_sems(pIdx,songIdx,:) = std(this_summary, [], 1) ./ sqrt(this_last-this_first+1);
        end
    end
    
    if uType == 1
        all_bout_means_mu = this_means;
        all_bout_sems_mu = this_sems;
    else
        all_bout_means = this_means;
        all_bout_sems = this_sems;
    end
end

all_floors = all_bout_means - all_bout_sems;
all_ceils = all_bout_means + all_bout_sems;

all_floors_mu = all_bout_means_mu - all_bout_sems_mu;
all_ceils_mu = all_bout_means_mu + all_bout_sems_mu;

fr_increase_su = squeeze(all_floors(1,:,:)) > squeeze(all_ceils(2,:,:));
fr_decrease_su = squeeze(all_floors(2,:,:)) > squeeze(all_ceils(1,:,:));
fr_change_su = (fr_increase_su | fr_decrease_su);
frc_su_allsongs = sum(fr_change_su,1) > 0;

fr_increase_mu = squeeze(all_floors_mu(1,:,:)) > squeeze(all_ceils_mu(2,:,:));
fr_decrease_mu = squeeze(all_floors_mu(2,:,:)) > squeeze(all_ceils_mu(1,:,:));
fr_change_mu = (fr_increase_mu | fr_decrease_mu);
frc_mu_allsongs = sum(fr_change_mu,1) > 0;

high_fr_drift = find(frc_su_allsongs);
high_fr_drift_mu = find(frc_mu_allsongs);

disp(['found  ' num2str(length(high_fr_drift)) ' su with fr drift'])

disp(['found  ' num2str(length(high_fr_drift_mu)) ' mu with fr drift'])

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

%% CHECK IN - HOW MANY NEURONS REMAIN? ARE THEY REASONABLE?

bad_units = union(union(low_fr, non_aud), union(high_fr_drift, high_spatial_drift));
good_unit_ids = setdiff(1:length(goodUnits), bad_units);
% good_unit_ids = goodUnits(good_unit_ids);

bad_units_mu = union(union(low_fr_mu, non_aud_mu), union(high_fr_drift_mu, high_spatial_drift_mu));
good_unit_ids_mu = setdiff(1:length(goodUnits_mu), bad_units_mu);
% good_unit_ids_mu = goodUnits_mu(good_unit_ids_mu); % try this???

disp(['there are ' num2str(length(good_unit_ids)) ' good single units'])
disp(['there are ' num2str(length(good_unit_ids_mu)) ' good multiunits'])


save(['D:\ibi_expmts\pipelines\' dataset_name '.mat'])
disp('early pipeline saved!')


%% save pipeline!
close all;
save(['D:\\ibi_expmts\pipelines\pipeline_' dataset_name '.mat'])
disp('pipeline complete!')