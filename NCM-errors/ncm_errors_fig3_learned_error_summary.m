% 1. load error detection pipeline output
% 2. plot learned error summary imagesc
% 3. across all birds, generate learned error response index histogram
% 4. for given bird,unit_type,generate example raster (probably need to
%   re-load individual bird pipeline output for this)


%% prep and set parameters
clear
close all
error_detection_pipeline = 'D:\error_expmts\multibird_summary_ncm_errors.mat';

example_info = [1, 390, 15; 3, 458, 11]; % bird, unit, perturbation. as many as you want

load(error_detection_pipeline)

rng(23) % just affects the sorting of summary fig
%% plot learned error imagesc
% extremely straightforward - all the heavy lifting is done in the error
% detection script
h=figure(4);

imagesc(flipud(big_fli_mat(randperm(length(big_fli_mat(:,1))),:)'))
title('Familiar Pert Summary')
ylabel('Trial')
xlabel('Neuron x Pert Type Combo')

print('figure_pieces/fig3_learnedSummary', '-dsvg', '-r300')
saveas(h, 'figure_pieces/fig3_learnedSummary.fig')

%% plot learned error response indices - color red for significant ids

a1_total = mean(big_resp_mat(:, 1:25),2);
odd_total = mean(big_resp_mat(:,26:50),2);
a2_total = mean(big_resp_mat(:,51:75),2);

err_idx = a2_total - ((a1_total +  odd_total)/2);
[cts, edges] = histcounts(err_idx);


h2=figure(20);
histogram(err_idx, 'BinEdges', edges)
hold on;
histogram(err_idx(big_fli_ids),'BinEdges',edges)
hold off;
title('Learned Error Indices - All Interactions')

print('figure_pieces/fig3_learnedIdxSummary', '-dsvg', '-r300')
saveas(h2, 'figure_pieces/fig3_learnedIdxSummary.fig')


% now we want 'max per neuron'
% find max pert response per neuron, given that typically neuron responds
% to just 1
num_neurons = length(err_idx)/18;
ei2 = zeros(num_neurons,1);

sig_nrn_ids = [];
for nIdx = 1:num_neurons
    start_idx = ((nIdx-1)*18)+1;
    end_idx = start_idx + 17;
    
    sig_nrn_bool = ~isempty(intersect(find(big_fli_ids >= start_idx), find(big_fli_ids <= end_idx)));
    if sig_nrn_bool
        sig_nrn_ids = [sig_nrn_ids; nIdx];
    end
    
    ei2(nIdx) = max(err_idx(start_idx:end_idx));
end

[cts2, edges2] = histcounts(ei2);
h3=figure(15);
histogram(ei2, 'BinEdges', edges2)
hold on;
histogram(ei2(sig_nrn_ids),'BinEdges',edges2)
hold off;
title('Learned Error Indices - Max per Neuron')

print('figure_pieces/fig3_learnedIdxSummary_mpn', '-dsvg', '-r300')
saveas(h3, 'figure_pieces/fig3_learnedIdxSummary_mpn.fig')


%% plot raster examples
% nb: this script is more involved

for ex_idx = 1:length(example_info(:,1))
    this_bird = example_info(ex_idx,1);
    load(the_good_birds{this_bird});
    
    %     pClass = {'_W', '_I'};
    %     pType = {'S', 'H', 'N'};
    
    % window is slightly before motif to slightly after
    
    motif_1_length_secs = motif_length;
    timeWindow_1 = [-0.1 motif_1_length_secs+0.225]; % see much after because we perturb final syllable
    
    pids = exp_info.pertSyllIDs_1;
    
    pert_time = exp_info.sub_inserts_1; % % same for both!
    bStart = (cell2mat(syllableInfo(pids(1),1)) - cell2mat(syllableInfo(1,1)));
    dStart = (cell2mat(syllableInfo(pids(2),1)) - cell2mat(syllableInfo(1,1)));
    fStart = (cell2mat(syllableInfo(pids(3),1)) - cell2mat(syllableInfo(1,1)));
    pert_time_1 = pert_time + [bStart dStart fStart];
    pert_time_2 = 0 + [bStart dStart fStart];
    
    bLen = cell2mat(syllableInfo(pids(1),4)) - cell2mat(syllableInfo(pids(1),1));
    dLen = cell2mat(syllableInfo(pids(2),4)) - cell2mat(syllableInfo(pids(2),1));
    fLen = cell2mat(syllableInfo(pids(3),4)) - cell2mat(syllableInfo(pids(3),1));
    
    
    
    if strcmp(stim_date, '210610')
        bl_1 = 1;
        pert = 476;
        hab = 2496;
        bl_2 = 2896;
    elseif strcmp(stim_date, '210613')
        bl_1 = 1;
        pert = 476;
        hab = 2488;
        bl_2 = 2888;
    end
    
    target_pert = pert_types{example_info(ex_idx,3)};
    target_syll = target_pert(2);
    target_type = target_pert(4);
    target_class = target_pert(5);
    
    switch target_syll
        case 'b'
            pert_time_2 = pert_time_2(1);
            syll_lens_1 = bLen;
        case 'd' 
            pert_time_2 = pert_time_2(2);
            syll_lens_1 = dLen;
        case 'f'
            pert_time_2 = pert_time_2(3);
            syll_lens_1 = fLen;
    end
    
    
    
    
    relevant_syll = syllableInfo; % this is perturb phase
    motif_starts = cell2mat(motifInfo(:,1));
    syllTypes = syllableInfo(:,3);
    mTypes = motifInfo(:,2);
    unique_types = unique(mTypes);
    
    % we now need to build structure with start times of each relevant
    % trial
    
    % what are relevant trials?
    relevant_types = {'@@@@@'}; % this just needs to be something absent from strings so that first cell is empty
    for uIdx = 1:length(unique_types)
        thisU = unique_types{uIdx};
        if contains(thisU, target_pert) % so is okay if target type doesn't exist...
            relevant_types{end+1} = thisU;
        end
    end
    
    
    start_times = cell(length(relevant_types),1);
    
    % first prehab
    for mIdx = bl_1:pert-1
        prevStart = motif_starts(19*(floor((mIdx-1)/19)) + 1);
        
        thisStart = motif_starts(mIdx);
        thisType = mTypes{mIdx};
        typeNum = [];
        relIdx = 1;
        while isempty(typeNum) && (relIdx <= length(relevant_types))
            if contains(thisType, relevant_types{relIdx})
                typeNum = relIdx;
            end
            relIdx = relIdx+1;
        end
        if ~isempty(typeNum)
            start_times{typeNum} = [start_times{typeNum}; thisStart];
            if typeNum == 2 % only add std once
                start_times{1} = [start_times{1}; prevStart]; % add previous standard motif
            end
        end
    end
    
    
    % pert phase
    for mIdx = pert:hab-1
        if mIdx-pert > 1
            prevStart = motif_starts(mIdx-1);
        end
        thisStart = motif_starts(mIdx);
        thisType = mTypes{mIdx};
        typeNum = [];
        relIdx = 1;
        while isempty(typeNum) && (relIdx <= length(relevant_types))
            if contains(thisType, relevant_types{relIdx})
                typeNum = relIdx;
            end
            relIdx = relIdx+1;
        end
        if ~isempty(typeNum)
            start_times{typeNum} = [start_times{typeNum}; thisStart];
            start_times{1} = [start_times{1}; prevStart]; % add previous standard motif
        end
    end
    
    
    % second prehab
    for mIdx = bl_2:length(motifTimingData(:,1))
        prevStart = motif_starts(19*(floor((mIdx-1)/19)) + 1);
        
        thisStart = motif_starts(mIdx);
        thisType = mTypes{mIdx};
        typeNum = [];
        relIdx = 1;
        while isempty(typeNum) && (relIdx <= length(relevant_types))
            if contains(thisType, relevant_types{relIdx})
                typeNum = relIdx;
            end
            relIdx = relIdx+1;
        end
        if ~isempty(typeNum)
            start_times{typeNum} = [start_times{typeNum}; thisStart];
            if typeNum == 2 % only add std once
                start_times{1} = [start_times{1}; prevStart]; % add previous standard motif
            end
        end
    end
    
    
    
    this_id = (example_info(ex_idx,2)); % not clearly grabbing correct unit? or else wrong pert
    theseSpikes = sp.st(sp.clu==goodUnits(this_id));

    perturb_response_plot_psth_vF(theseSpikes, start_times, this_id, timeWindow_1, micData, target_type, pert_time_2, syll_lens_1, 1, ex_idx)
    
    h5 = gcf;
    print(['figure_pieces/fig3_learnedExample' '_' num2str(ex_idx)], '-dsvg', '-r300')
    saveas(h5, ['figure_pieces/fig3_learnedExample'  '_' num2str(ex_idx) '.fig'])
end


disp('i think im done')
