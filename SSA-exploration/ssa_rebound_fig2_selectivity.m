% figure 2 - memory and habituation, song selectivity

%% prep and set variables

clear
close all
pipeline_out = 'D:\ssa_expmts\multibird_summary_ssa.mat';

%%
load(pipeline_out)

%% scatters for all songs

% unoptimized plotting, does 1 at a time for some reason

impt_motifs = [1 125 501 550;
    126 250 551 600;
    251 375 601 650;
    376 500 651 700];

all_mems = [];

for songIdx = 1:4
    tkm = impt_motifs(songIdx,:);
    
    this_resp_p1 = total_motif_resp(tkm(1):tkm(2),:);
    this_resp_p2 = total_motif_resp(tkm(3):tkm(4),:);
    spike_mean = mean([this_resp_p1; this_resp_p2],1);
    
    low_fr = find(spike_mean<0.5);
    ok_fr = setdiff(1:length(spike_mean),low_fr);
        
    
    p1_early = mean(this_resp_p1(1:5,ok_fr));
    p1_late = mean(this_resp_p1(end-4:end,ok_fr));
    p2_early = mean(this_resp_p2(1:5,ok_fr));
    
    this_hab = (p1_early - p1_late) ./ (p1_early + p1_late);
    this_mem = ((p1_early - p2_early) ./ (p1_early + p2_early));
    these_res_hab_mems = [this_hab this_mem];
    
    
    figure(1214411)
    subplot(2,2,songIdx)
    sc1 = scatter(this_hab, this_mem,...
        'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');%,...
    %'MarkerSize',20);
    sc1.MarkerFaceAlpha = 0.2;
    sc1.MarkerEdgeAlpha = 0.2;
    
    
    xlabel('Habituation')
    ylabel('Memory')
    title(['Restricted, Song ' num2str(songIdx)])
    refline(0,0)
    line([0 0], [-1 1])
    axis([-1.05 1.05 -1.05 1.05])
    
end

h2 = gcf;
print('figure_pieces/fig2_ssaScatters', '-dsvg', '-r300')
saveas(h2, 'figure_pieces/fig2_ssaScatters.fig')



%% song selectivity plots
% red dots for ssa units

impt_motifs = [1 125 501 550;
    126 250 551 600;
    251 375 601 650;
    376 500 651 700];

max_resp = [];
max_hab = [];
max_mem = [];
    
for songIdx = 1:4
    tkm = impt_motifs(songIdx,:);
    this_resp_p1 = total_motif_resp(tkm(1):tkm(2),:);
    
    this_max = max(this_resp_p1,[],1);
    
    p1_early = mean(this_resp_p1(1:5,:),1);
    p1_late = mean(this_resp_p1(end-4:end,:),1);
    p2_early = mean(this_resp_p2(1:5,:));
    
    this_hab = (p1_early - p1_late) ./ (p1_early + p1_late);
    this_mem = (p1_early - p2_early) ./ (p1_early + p2_early);
    
    max_resp = [max_resp; this_max];
    max_hab = [max_hab; this_hab];
    max_mem = [max_mem; this_mem];
end

resp_xval = max(max_resp, [], 1);
resp_yval = min(max_resp, [], 1);

hab_xval = max(max_hab, [], 1);
hab_yval = min(max_hab, [], 1);

mem_xval = max(max_mem, [], 1);
mem_yval = min(max_mem, [], 1);

figure(1240)
sc1 = scatter(resp_xval, resp_yval,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
sc1.MarkerFaceAlpha = 0.2;
sc1.MarkerEdgeAlpha = 0.2;
hold on;
refline(1,0)
hold off;
title(['Song Response Selectivity'])
xlabel('Maximal peak song response (hz)')
ylabel('Minimal peak song response (hz)')
h = gcf;
print('figure_pieces/fig2_songRespSelectivity', '-dsvg', '-r300')
saveas(h, 'figure_pieces/fig2_songRespSelectivity.fig')


figure(1234)
sc1 = scatter(hab_xval, hab_yval,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
sc1.MarkerFaceAlpha = 0.2;
sc1.MarkerEdgeAlpha = 0.2;
hold on;
refline(1,0)
hold off;
title(['Song Habituation Selectivity'])
xlabel('Maximal song habituation')
ylabel('Minimal song habituation')
h = gcf;
print('figure_pieces/fig2_songHabSelectivity', '-dsvg', '-r300')
saveas(h, 'figure_pieces/fig2_songHabSelectivity.fig')

figure(125)
sc1 = scatter(mem_xval, mem_yval,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
sc1.MarkerFaceAlpha = 0.2;
sc1.MarkerEdgeAlpha = 0.2;
hold on;
refline(1,0)
hold off;
title('Song Memory Selectivity')
xlabel('Maximal song memory')
ylabel('Minimal song memory')
h = gcf;
print('figure_pieces/fig2_songMemSelectivity', '-dsvg', '-r300')
saveas(h, 'figure_pieces/fig2_songMemSelectivity.fig')

