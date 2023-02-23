%%% pca + mismatch network

% y comes into pca network, which runs one iter. (immediately clusters) 
% output of pca network is unary context signal
% y and context signal are then fed into mismatch network

% pca network starts with relatively low excitability, and excitability
% increases every time total mismatch is over some threshold


%%% Important notes:
%%%     -this specific implementation could fail if we're unlucky with the
%%%     intial directions of the weights (ie if no weight is close enough
%%%     to input vector, we'd get into trouble)
%%%     -this version has been edited to work with very high dim inputs!
%%%     -v2 tries to introduce normalizations to mm weights
%%%
%%%

%%% Things to check!
%%%     -Mismatch error should be decreasing over time, although that may
%%%     take a while
%%%     -Rank of W should initially be close to 100 (or whatever
%%%     initialized), but should end up close to number of clusters
%%%     -Weight magnitudes (for W and M) shouldn't be too crazy
%%%     -Cluster IDs should be consistent - after running network once,
%%%     re-running on same points should give same clusters (unless network
%%%     takes very long to converge)
%%%     -Mismatch error should be much lower if same data run on trained network
%%%     -Context signal should be relatively sparse!

clear 
close all

h1 = figure(1);
set(h1, 'Position', [2209         -31         512       384])
h2 = figure(5);
set(h2, 'Position', [1370        -319         822         672])
% figure(10)
% figure(25)
pause
%% Configs
seed = 12345; %10 was working for many of these demos %56; %23, 10 for 100D 10 clusters
    % seed choice affects initial weight direction, which affects peaks of
    % error plots
%seed = 1123;
trackVars = 1; %turning on slows things WAY DOWN


% Select Dataset
%load('genData1.mat')
% load('genData2D_4c_200ppc.mat')
load('genData2D_3c_500ppc.mat')
% load('genData100D_10.mat')
% load('genData_realSpec_try2.mat')
% load('genData_realSpec_5c.mat')


% Define how data will be presented
twoClust_dm_unshuff = allPts_dm_unshuff(1:100, :);
twoClust_unshuff = allPts_unshuff(1:100, :);

shuffIdx = 1:100;%randperm(500);

part1_dm = twoClust_dm_unshuff(shuffIdx,:);
part1 = twoClust_unshuff(shuffIdx, :);

% bigY is input fed into clustering network - note that it's de-meaned
bigY = [allPts_dm_unshuff']; %[part1_dm' allPts_dm'];%_unshuff'];% allPts_dm' allPts_dm'];%(1:500,:)'];% allPts_dm'];

% mmY is input fed into mismatch network - just a non-de-meaned version of
% bigY
mmY = [allPts_unshuff']; %[part1' allPts'];%_unshuff'];% allPts' allPts'];%(1:500,:)'];%  allPts'];

% PCA CONFIGS
pcaNet.changeThresh = 1e-4; % for output convergence
initWeightMag = 0.1; %0.8 for some demos % should be 0.1
pcaNet.capW = 10;%500; %max weight to each output cell -- CHANGING THIS DIDN'T MAKE A DIFFERENCE
pcaNet.inhibCap = 10; %NEED 10 FOR SYNTHETIC, 25 FOR REAL... -- CHANGING THIS DOES MATTER

pcaNet.etaW = 0.05; %originally each was 0.5 -- MAKING THIS TOO HIGH DOESN'T MATTER MUCH
pcaNet.etaM = 0.05; % originally 0.5 -- CHANGING THIS DOESN'T MATTER MUCH
pcaNet.maxM = 1; 
pcaNet.maxW = 1;

% MISMATCH CONFIGS
mmNet.eta = 0.025; % CHANGING THIS DRASTICALLY AFFECTS SHAPE OF ERROR PLOTS
mmNet.thresh = 0.05; % CHANGING THIS AFFECTS SHAPE OF ERROR PLOTS

% MISTMATCH STRUCTURE CONFIGS
mmNet.signed_synapses = 1; %force positive weight coeffs 
mmNet.c_plastic = 1; % 
mmNet.y_plastic = 0;

% MISMATCH ARCHITECTURE CONFIGS
yE_type = 'rand'; %'rand' or 'randnorm'
cE_type = 'rand';  
yI_type = 'rand';
cI_type = 'rand';
yR_type = 'direct'; %'rand' or 'direct' 
cR_type = 'direct'; %'rand' or 'direct'

% MISMATCH NETWORK CONFIGS
nCells_y = size(mmY,1); %one input per input dim
nCells_c = 10; % (500 works on synthetic) 10 clusters in 100D case (nb 100d is a lot)
nCells_Ny = 5; %not super sure how high or low this needs to be
nCells_Nc = 5; %but 100 and 100 should be enough to give us convex cone
iterations = size(bigY,2);


% Total Configs -- this doesn't matter, we're always learning
% (learningThresh = 0)
sigmaThresh = 2.5; %1 works for seed 10 %2.5;%2.5;%%0.25;%25; %1 worked well 
deltaLearn = 0.2; %(0.1) % no buffer -> incorrect clustering (not obvious why needed for shuffled data)
learningThresh = 0;%1; %make it deltaLearn to basically eliminate deltaLearn %0.5; %always learning -> incorrect clustering
pcaNet.sigmaThresh = sigmaThresh;


pcaNet.learning = 1; %this is set to 0 if mismatch is low, 1 if high
learningSig = 1;
rng(seed) 
%% PCA - Setup
W_init = initWeightMag*randn(nCells_c, nCells_y);
M_init = 0*initWeightMag*rand(nCells_c); %initially 0

for idx = 1:nCells_c
    M_init(idx, idx) = 0; %nrns don't drive selves
end

D_init = zeros(nCells_c,1);
pcaNet.bigC = zeros(nCells_c, iterations);

pcaNet.W = W_init;
pcaNet.M = M_init;
pcaNet.D = D_init;

c = updateC_v5(...
    pcaNet.W, pcaNet.M, bigY(:,1), pcaNet.changeThresh);

[maxx, thisCluster] = max(c);

if maxx > 0
    pcaNet.clusters = thisCluster;
    pcaNet.bigC(thisCluster, 1) = c(thisCluster);
else
    pcaNet.clusters = 0;
end

cT = zeros(size(c));
cT(thisCluster) = c(thisCluster); % c at timestep t
y = bigY(:,1);

rng(seed);
%% Mismatch Setup
switch yE_type
    case 'rand'
        mmNet.we_yn = (rand(nCells_Ny, nCells_y))./nCells_y;
    case 'randnorm'
        mmNet.we_yn = (randn(nCells_Ny, nCells_y))./nCells_y;
end

switch yI_type
    case 'rand'
        mmNet.wi_yn = (rand(nCells_Nc, nCells_y))./nCells_y;
    case 'randnorm'
        mmNet.wi_yn = (randn(nCells_Nc, nCells_y))./nCells_y;
end

switch yR_type
    case 'rand'
        mmNet.r_yn = rand(nCells_y);
    case 'direct'
        mmNet.r_yn = eye(nCells_y);
end
    
switch cE_type
    case 'rand'
        mmNet.we_cn = (rand(nCells_Nc, nCells_c));%./nCells_c;
        %
    case 'randnorm'
        mmNet.we_cn = (randn(nCells_Nc, nCells_c));%./nCells_c;
end

switch cI_type
    case 'rand'
        mmNet.wi_cn = (rand(nCells_Ny, nCells_c));%./nCells_c;
        %
    case 'randnorm'
        mmNet.wi_cn = (randn(nCells_Ny, nCells_c));%./nCells_c;
end

switch cR_type
    case 'rand'
        mmNet.r_cn = rand(nCells_c);
    case 'direct'
        mmNet.r_cn = eye(nCells_c);
end

timesteps = length(bigY(1,:));

if trackVars == 1
    mmNet.Vs_y = zeros(nCells_Ny, timesteps);
    mmNet.Vs_c = zeros(nCells_Nc, timesteps);
    mmNet.Fs_y = zeros(nCells_Ny, timesteps);
    mmNet.Fs_c = zeros(nCells_Nc, timesteps);


    mmNet.yWs_e = zeros(nCells_Ny, nCells_y, timesteps+1);
    mmNet.cWs_e = zeros(nCells_Nc, nCells_c, timesteps+1);
    mmNet.yWs_i = zeros(nCells_Nc, nCells_y, timesteps+1);
    mmNet.cWs_i = zeros(nCells_Ny, nCells_c, timesteps+1);

    mmNet.wyChanges_e = zeros(nCells_Ny, nCells_y, timesteps);
    mmNet.wcChanges_e = zeros(nCells_Nc, nCells_c, timesteps);
    mmNet.wyChanges_i = zeros(nCells_Nc, nCells_y, timesteps);
    mmNet.wcChanges_i = zeros(nCells_Ny, nCells_c, timesteps);

    mmNet.yWs_e(:,:,1) = mmNet.we_yn;
    mmNet.cWs_e(:,:,1) = mmNet.we_cn;
    mmNet.yWs_i(:,:,1) = mmNet.wi_yn;
    mmNet.cWs_i(:,:,1) = mmNet.wi_cn;
end

mmNet.errors_y = zeros(timesteps, 1);
mmNet.errors_c = zeros(timesteps, 1);
mmNet.allErrors = zeros(timesteps,1);

% we've already done 1 iter of pca, so here do 1 iter of mm
mmNet = mismatchIter_v2(cT, mmY, 1, mmNet, trackVars);  


disp('setup complete')
%% Run through algorithm (make functions for PCAIter and MNIter)
for ts_idx = 2:iterations
    pcaNet = pcaIter_v6(bigY, ts_idx, pcaNet);

    
    pcaC = pcaNet.bigC(:, ts_idx);

    pcaC = pcaC>eps;
    
    mmNet = mismatchIter_v2(pcaC, mmY, ts_idx, mmNet, trackVars);
    
    sigmaNode = mmNet.allErrors(ts_idx);
    if sigmaNode > sigmaThresh
        learningSig = learningSig + deltaLearn;
        if learningSig > 1
            learningSig = 1;
        end
    else
        learningSig = learningSig - deltaLearn;
        if learningSig < 0
            learningSig = 0;
        end
    end
    
    if learningSig >= learningThresh
        pcaNet.learning = 1;
    else 
        pcaNet.learning = 0;
    end
    
%     agnoPlotting_2D(bigY, mmNet, pcaNet, ts_idx, iterations)
%     agnoPlotting_100D(bigY, mmNet, pcaNet, ts_idx, iterations)
end

agnoPlotting_2D(bigY, mmNet, pcaNet, ts_idx, iterations)

%% Plotting
% 
% % Cluster ID Distribution
% figure(2)
% subplot(2,2,[1 2])
% histogram(pcaNet.clusters(1:iterations), [-0.5:1:nCells_c+0.5])
% title('Inputs per Cluster - Total Dataset')
% xlabel('Cluster ID')
% ylabel('Number of Inputs')
% 
% subplot(2,2,3)
% histogram(pcaNet.clusters(1:iterations/2), [-0.5:1:nCells_c+0.5])
% title('Inputs per Cluster - First Half')
% xlabel('Cluster ID')
% ylabel('Number of Inputs')
% 
% subplot(2,2,4)
% histogram(pcaNet.clusters((iterations/2)+1:end), [-0.5:1:nCells_c+0.5])
% title('Inputs per Cluster - Second Half')
% xlabel('Cluster ID')
% ylabel('Number of Inputs')