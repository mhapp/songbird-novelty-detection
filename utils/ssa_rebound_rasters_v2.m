function ssa_rebound_rasters_v2(spikeTimes, alignmentTimes, unitData, timeWindow, micData, recDate, uType)
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
figPath = 'D:\ssa_expmts\figures'; % UPDATE THIS
fs = 40000;


    these_ATs = alignmentTimes;
    this_window = timeWindow;
    
    micRows = zeros(length(these_ATs),round((this_window(2)-this_window(1))*fs)+1);
    for trialIdx = 1:length(these_ATs)
        thisIdxRange = round((these_ATs(trialIdx)+this_window)*fs);
        micRows(trialIdx,:) = micData(thisIdxRange(1):thisIdxRange(2))';
    end
    
    avgMic = micRows(1,:);
    winLen = round(.005*fs); % 1ms window on spec
    freqs = 1:10:8001; % step to 8kHz in 10Hz bins
    
    figure(18);
    subplot(4,1,1)
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
    these_ATs = alignmentTimes;
    
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
    
    
    figure(18);
    subplot(4,1,[2 3 4])
    plot(rasterX, rasterY, 'k', 'MarkerSize', 1);
    xlim(timeWindow);
    xlabel('Time (s)')
    ylim([1 length(these_ATs)+1]);
%     makepretty;
%     if pertIdx ~= 1
%         pTime = pert_time(pertIdx);
%         pLen = pert_lens(pertIdx);
%         hold on;
%         lim_of_y = ylim;
%         if songIdx == 1
%             yl = [lim_of_y(1) 75+1];
%         else
%             yl = [lim_of_y(1) 50+1];
%         end
%         v1 = [pTime yl(1); pTime yl(2); pTime+pLen yl(1); pTime+pLen yl(2)];
%         f1 = [1 2 4 3];
%         
%         patch('Faces', f1, 'Vertices', v1, 'FaceColor', 'k', 'FaceAlpha', 0.25);
        
%         if songIdx == 1
%             yl = [75 lim_of_y(2)];
%         else
%             yl = [50 lim_of_y(2)];
%         end
%         
%         v1 = [pTime yl(1); pTime yl(2); pTime+pLen yl(1); pTime+pLen yl(2)];
%         f1 = [1 2 4 3];
%         
%         patch('Faces', f1, 'Vertices', v1, 'FaceColor', 'k', 'FaceAlpha', 0.25);
%         
            
        set(gca,'YTick',[])
        h1 = refline(0,6);
        h1.Color = 'k';
        h2 = refline(0,16);
        h2.Color = 'k';
        h3 = refline(0,11);
        h3.Color = 'r';
        hold off;
%     end
    
        ylabel(['event number'])
    
    
    xlabel(['Time (mot. on, s)'])
%     if pertIdx == 1
% %         ylabel(['Avg. Firing Rate (Hz)'])
%     end
    xlim(timeWindow)