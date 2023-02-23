%%% figure 1 - ssa basic summary

% basic message is, let's investigate SSA a little more closely. First,
% here's a summary, along with the finding that memory component is present
% in most habituating neurons, and they even correlate in a way that makes
% sense (and it's therefore probably fine to ignore this memory parameter
% when thinking about SSA).

%% prep and set variables

clear
close all
pipeline_out = 'D:\ssa_expmts\multibird_summary_ssa.mat';

psth_bird_id = 2;
psth_unit_id = 22;
psth_unit_type = 'su';


%% ssa index summary    

load(pipeline_out)

big_ssa_mat = [];

p1_resp = total_motif_resp(1:500,:);
init_resps = p1_resp([1:5 126:130 251:255 376:380],:);
max_resps = max(init_resps,[],1);

% eliminate units with no initial resp to any song, to avoid infs
% NB: this misses units that potentiate, but we're not looking at those
% for this basic confirmatory analysis
valid_ids = find(max_resps > 0);
p1_resp_valid = p1_resp(:,valid_ids);
max_resps_valid = max_resps(valid_ids);

% normalize responses to strongest initial response
p1_norm = p1_resp_valid ./ max_resps_valid;

% get responses for each unit for each of 4 std songs
s1p1_resp = p1_resp(1:125,valid_ids);
s2p1_resp = p1_resp(126:250,valid_ids);
s3p1_resp = p1_resp(251:375,valid_ids);
s4p1_resp = p1_resp(376:500,valid_ids);


for unitIdx = 1:length(valid_ids)
    this_ssa_vec = [];
    for songIdx = 1:4 % at least one will be nonzero because of valid_units
        eval(['this_resp = s' num2str(songIdx) 'p1_resp(:,unitIdx);'])
        this_init = mean(this_resp(1:5));
        this_end = mean(this_resp(121:125));
        test_ssa = (this_init-this_end)/(this_init + this_end);
        if ~isinf(test_ssa)
            this_ssa_vec = [this_ssa_vec; test_ssa];
        else
            this_ssa_vec = [this_ssa_vec; 0];
        end
    end
    big_ssa_mat = [big_ssa_mat; max(this_ssa_vec)]; %max(this_ssa_vec) or not
end


h=figure; 
histogram(max(big_ssa_mat,0),20)
xlabel('SSA Index')
ylabel('Number of Neurons')
title('Distribution of SSA Index')
xlim([-0.2 1.2])

print('figure_pieces/fig1_ssaIdx', '-dsvg', '-r300')
saveas(h, 'figure_pieces/fig1_ssaIdx.fig')


%% plot psth

% unoptimized and counter-intuitive, sorry. 

% transition times for the stim
s1s2 = 125;
s2s3 = 250;
s3s4 = 375;
s4s1b = 500;
s1bs2b = 550;
s2bs3b = 600;
s3bs4b = 650;
std_end = 700;

if strcmp(psth_unit_type, 'su')
    if psth_bird_id > 1
        base = sum(num_goods(1:psth_bird_id-1));
    else
        base = 0;
    end
    actual_id = base + psth_unit_id;
    this_row = total_motif_resp(1:std_end,actual_id);
else
    this_row = total_motif_resp_mu(1:std_end,actual_id);
end

h2=figure(1006);
subplot(4,2,[1 2 3 4])
plot(1:s1s2, this_row(1:s1s2), 'r-')
hold on;
plot(s1s2+1:s1s2+50, this_row(s4s1b+1:s1bs2b), 'r.')
newstart = s1s2+51;
next_piece = this_row(s1s2+1:s2s3);
newend = newstart+length(next_piece)-1;
plot(newstart:newend, this_row(s1s2+1:s2s3), 'b-')

newstart = newend+1;
next_piece = this_row(s1bs2b+1:s2bs3b);
newend = newstart+length(next_piece)-1;
plot(newstart:newend, next_piece, 'b.')

newstart = newend+1;
next_piece = this_row(s2s3+1:s3s4);
newend = newstart+length(next_piece)-1;
plot(newstart:newend, next_piece, 'g-')

newstart = newend+1;
next_piece = this_row(s2bs3b+1:s3bs4b);
newend = newstart+length(next_piece)-1;
plot(newstart:newend, next_piece, 'g.')

newstart = newend+1;
next_piece = this_row(s3s4+1:s4s1b);
newend = newstart+length(next_piece)-1;
plot(newstart:newend, next_piece, 'k-')

newstart = newend+1;
next_piece = this_row(s3bs4b+1:std_end);
newend = newstart+length(next_piece)-1;
plot(newstart:newend, next_piece, 'k.')

hold off;
ylabel('resp strength (hz)')
title(['Std Motif Responses'])
xlim([0 std_end])

subplot(4,2,5)
plot(this_row(1:25), 'r*-')
xlim([1 25])
title('First Responses to Std 1')
ylim([min(0, min(this_row)) max(this_row)+10])

subplot(4,2,6)
plot(this_row(s1s2+1:s1s2+25), 'b*-')
xlim([1 25])
title('First Responses to Std 2')
ylim([min(0, min(this_row)) max(this_row)+10])

subplot(4,2,7)
plot(this_row(s2s3+1:s2s3+25), 'g*-')
xlim([1 25])
title('First Responses to Std 3')
ylim([min(0, min(this_row)) max(this_row)+10])

subplot(4,2,8)
plot(this_row(s3s4+1:s3s4+25), 'k*-')
xlim([1 25])
title('First Responses to Std 4, p1')
ylim([min(0, min(this_row)) max(this_row)+10])

print('figure_pieces/fig1_ssaExample', '-dsvg', '-r300')
saveas(h2, 'figure_pieces/fig1_ssaExample.fig')


%% imagescs showing p1 habituation for all songs

% slightly different for these plots
impt_motifs = [1 25;
    126 150;
    251 275;
    376 400];

% come up with p1 hab sort

for songIdx = 1:4
    tkm = impt_motifs(songIdx,:);
    
    
    this_resp_p1 = total_motif_resp(tkm(1):tkm(2), :);
    this_resp_norm = this_resp_p1 ./ max(this_resp_p1, [], 1);
    
        
    p1_early = mean(this_resp_p1(1:5,:));
    p1_late = mean(this_resp_p1(end-4:end,:));
    
    this_valid = find(p1_early>0);
    
    hab_idx = (p1_early - p1_late) ./ (p1_early + p1_late);
    
    hab_idx = hab_idx(this_valid);
    resp_valid = this_resp_norm(:,this_valid);
    
    
    [~,this_hab_sort] = sort(hab_idx,'descend');
    
    
    
    
    figure(88)
    subplot(2,2,songIdx)
    imagesc(resp_valid(:,this_hab_sort))
    xlabel('neuron id.')
    ylabel('trial')
    title(['Habituation Summary, Song ' num2str(songIdx)])
end

print('figure_pieces/fig1_habSummary', '-dsvg', '-r300')
saveas(h2, 'figure_pieces/fig1_habSummary.fig')

disp('all done!')