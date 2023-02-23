function perturb_response_plot_psth_vF(spikeTimes, alignmentTimes, unitData, timeWindow, micData, thisType, pert_time, syll_lens, songIdx, figNum)
%%%%
% this plots a figure showing unit response to different perturbations
%%%%
% spike times is a simple vector of times when spike occurs (sec)
% alignment Times is cell array where each entry contains a vector of start
%   times for all perturbations of a given type
% unitData is unit #
% window is time (sec) before and after 0 for each plot
% micData is microphone data, used to generate spectrogram
%%%%
% function generates single figure window:
% - A DIFFERENT COLUMN FOR EACH PERTURBATION TYPE! THIS COULD BE AS MANY AS
%   19 COLUMNS! (but only 13 for first try because we messed up
%   perturbations)
% - each column has spectrogram at top, raster in middle, psth at bottom
%%%%
fs = 40000;


% do stuff based on type
if contains(thisType, 'W')
    pert_time = pert_time;% - 0.05; %0.05 commented out
    pert_lens = syll_lens + 0.2;
elseif contains(thisType, 'I')
    pert_time = pert_time + 0.05;
    pert_lens = 0.05 + 0.2;
end



%% first generate spectrograms
% find relevant times in mic data and compute average
for pertIdx = 1:length(alignmentTimes)
    these_ATs = alignmentTimes{pertIdx};
    spec_ATs = these_ATs(1:25); % to avoid blending
    micRows = zeros(length(spec_ATs),round((timeWindow(2)-timeWindow(1))*fs)+1);
    for trialIdx = 1:length(spec_ATs)
        thisIdxRange = round((these_ATs(trialIdx)+timeWindow)*fs);
        micRows(trialIdx,:) = micData(thisIdxRange(1):thisIdxRange(2))';
    end
    
    avgMic = mean(micRows,1);
    winLen = round(.005*fs); % 1ms window on spec
    freqs = 1:10:8001; % step to 8kHz in 10Hz bins
    
    figure(figNum);
    subplot(4,length(alignmentTimes),pertIdx)
    spectrogram(avgMic, winLen, floor(winLen/2), freqs, fs, 'yaxis');
    % this stuff to be inserted with actual plotting
    lim = caxis;
    span = diff(caxis);
    newspan = span/2;
    newlim1 = lim(2) - newspan;
    caxis([newlim1 lim(2)])
    set(gca, 'Xticklabel', []);
    set(gca, 'Yticklabel', []);
    title(['Unit ' num2str(unitData)])
    xlabel([])
    ylabel([])
    colorbar('off')
    
    %% now do PSTHs and ssa idx
    [~, bins, ~, ~, ~, ba] = psthAndBA(spikeTimes, these_ATs, timeWindow, 0.001); %1ms bin size
    
    [tr,b] = find(ba);
    [rasterX,yy] = rasterize(bins(b));
    rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); % yy is of the form [0 1 NaN 0 1 NaN...] so just need to add trial number to everything
    
    % scale the raster ticks
    rasterY(2:3:end) = rasterY(2:3:end)+floor(numel(these_ATs)/100);
    
    % PSTH smoothing filter
%     gw = gausswin(round(*6),3); % 5 ms gaussian window
%     smWin = gw./sum(gw);
%     
%     % smooth ba
%     baSm = conv2(smWin,1,ba', 'same')'./.001;
    psthSm = mean(ba); % mean(baSm) when smoothing
    
    if pertIdx == 1
        std_bins = bins;
        std_psth = psthSm;
    end
    
    
    figure(figNum);
    subplot(4,length(alignmentTimes),...
        [pertIdx+length(alignmentTimes), pertIdx+2*length(alignmentTimes),...
        pertIdx+3*length(alignmentTimes)])
    plot(rasterX, rasterY, 'k', 'MarkerSize', 1);
    xlim(timeWindow);
    xlabel('Time (s)')
    ylim([0 length(these_ATs)+1]);
%     makepretty;
%      if pertIdx ~= 1
        pTime = pert_time;
        pLen = pert_lens;
        hold on;
        lim_of_y = ylim;
        if songIdx == 1
            yl = [lim_of_y(1) 75+1];
        else
            yl = [lim_of_y(1) 50+1];
        end
        v1 = [pTime yl(1); pTime yl(2); pTime+pLen yl(1); pTime+pLen yl(2)];
        f1 = [1 2 4 3];
        
        patch('Faces', f1, 'Vertices', v1, 'FaceColor', 'k', 'FaceAlpha', 0.25);
             
        h1 = refline(0,25.5);
        h1.Color = 'r';
        if songIdx == 1
            h3 = refline(0,50.5);
            h3.Color = 'r';
        end
        hold off;
%      end
    
    if pertIdx == 1
        ylabel(['event number'])
    end
    
    
    xlabel(['Time (mot. on, s)'])
%     if pertIdx == 1
% %         ylabel(['Avg. Firing Rate (Hz)'])
%     end
    xlim(timeWindow)
end