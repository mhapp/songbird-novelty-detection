%% v14 wants two chains for odorants A1 B1, A2 B2
% very close...
% i need to 1) create separate variables for 1 and 2, to make the temporal
% jitter easier and 2) make the temporal jitter more extreme and varied

% now change the learning rule!

clear
close all
figure(5)
set(gca, 'Position', [0.1156    0.1100    0.6889    0.8150])
figure(10)
set(gca, 'Position', [0.1300    0.1148    0.7750    0.1312])
% pause

%%% configs

seed =2509;% 23; %2509; %23;

winSize = 40;

cdNet.nCD = 120;
cdNet.nMN = 20;
cdNet.nPD = 4;
cdNet.gamma = 0.0;  %0 before %  0.045;%045; %0.02 %recurrent inhibition
cdNet.beta = 0.025;%0.025;15;%0.025;%15; %0.115;  %feed forward inhibition
cdNet.epsilon = 0.2; %relative strength of heterosynaptic LTD
cdNet.tau = 30; % 30 %time constant of adaptation
cdNet.alpha = 20; % 50 % strength of adaptation
cdNet.max_wCD = 1; %max weight between cd neurons
cdNet.m = 25;% 25; %10-16 ideal %desired internal synapses per cd nrn
cdNet.Wmax = cdNet.max_wCD * cdNet.m;
cdNet.nSeed = 10; % 5 % number of seed neurons per pd output
cdNet.pin = 0; %.01; 

seedThresh = 0;
seedDrive = 1;

cdNet.fCoef = 0; 
cdNet.vCoef = 1;% - cdNet.fCoef;


cdNet.eta_cd = 0.25;%15; %0.001; %25;%25;%0.5; %initial learning rate
deltaM = -0.01;%.025;% -0.05;% -0.5; %change in desired # synapses per mismatch

% m = 25, deltaM = -0.05 is original setting

mmNet.mmThresh = 0; %mismatch below this is treated as 0
mmNet.sigmaThresh = 3; %total MM above this increases eta (not really)
mmNet.eta_mm = 0.1; %0.1%0.5; % learning rate for mismatch network


rng(seed);
%%% initialize weight mats
w0 = 0.0375 .* rand(cdNet.nCD); %0.25*(2*cdNet.Wmax/cdNet.nCD)

cdNet.w_CD_CD = w0;%0.1*rand(cdNet.nCD, cdNet.nCD);
cdNet.w_PD_CD = zeros(cdNet.nCD, cdNet.nPD);
cdNet.w_PD_CD(1:cdNet.nSeed, 1) = 1;
cdNet.w_PD_CD(cdNet.nSeed+1:2*cdNet.nSeed, 2) = 1;
cdNet.w_PD_CD((2*cdNet.nSeed)+1:3*cdNet.nSeed, 3) = 1;
cdNet.w_PD_CD((3*cdNet.nSeed)+1:4*cdNet.nSeed, 4) = 1;

mmNet.w_CD_MN = zeros(cdNet.nMN, cdNet.nCD); %start as zeros
mmNet.w_PD_MN = randn(cdNet.nMN, cdNet.nPD);


%%% set input activity
patternIn = seedThresh*ones(4, 10);
patternIn(1, 1) = seedDrive;
patternIn(2, 7) = seedDrive;
corePatt1 = patternIn;

patternIn = seedThresh*ones(4, 10);
patternIn(3, 1) = seedDrive;
patternIn(4, 7) = seedDrive;
corePatt2 = patternIn;

cdPatt1x1 = [seedThresh*ones(4,1) corePatt1 seedThresh*ones(4,8) corePatt2 seedThresh*ones(4,11)];
cdPatt1x2 = [seedThresh*ones(4,4) corePatt1 seedThresh*ones(4,7) corePatt2 seedThresh*ones(4,9)];
cdPatt1x3 = [seedThresh*ones(4,2) corePatt1 seedThresh*ones(4,11) corePatt2 seedThresh*ones(4,7)];
cdPatt1x4 = [seedThresh*ones(4,3) corePatt1 seedThresh*ones(4,9) corePatt2 seedThresh*ones(4,8)];
cdPatt1x5 = [seedThresh*ones(4,5) corePatt1 seedThresh*ones(4,7) corePatt2 seedThresh*ones(4,8)];

altCore = seedThresh*ones(4,40);
altCore(1,3) = seedDrive;
altCore(4,9) = seedDrive;
altCore(3,22) = seedDrive;


patternIn = zeros(4, 10);
patternIn(1,1:4) = 1;
patternIn(2, 7:10) = 1;
mmCore1 = patternIn;

patternIn = zeros(4,10);
patternIn(3,1:4) = 1;
patternIn(4, 7:10) = 1;
mmCore2 = patternIn;

mmPatt1x1 = [zeros(4,1) mmCore1 zeros(4,8) mmCore2 zeros(4,11)];
mmPatt1x2 = [zeros(4,4) mmCore1 zeros(4,7) mmCore2 zeros(4,9)];
mmPatt1x3 = [zeros(4,2) mmCore1 zeros(4,11) mmCore2 zeros(4,7)];
mmPatt1x4 = [zeros(4,3) mmCore1 zeros(4,9) mmCore2 zeros(4,8)];
mmPatt1x5 = [zeros(4,5) mmCore1 zeros(4,7) mmCore2 zeros(4,8)];

mmAlt = zeros(4,40);
mmAlt(1,3:6) = 1;
mmAlt(4, 9:12) = 1;
mmAlt(3, 22:25) = 1;

corePattCD = [cdPatt1x1 cdPatt1x2 cdPatt1x3 cdPatt1x4 cdPatt1x5];
corePattMM = [mmPatt1x1 mmPatt1x2 mmPatt1x3 mmPatt1x4 mmPatt1x5];

cdNet.input = [repmat(corePattCD, [1 39]) repmat(altCore, [1 5])];
mmNet.input = [repmat(corePattMM, [1 39]) repmat(mmAlt, [1 5])];

iterations = size(cdNet.input, 2);
%%% final prep for iterations

cdNet.old_cd = zeros(cdNet.nCD, 1);
cdNet.totalVMat = zeros(cdNet.nCD, iterations);
cdNet.old_v = zeros(cdNet.nCD, 1);
cdNet.old_y = zeros(cdNet.nCD, 1);
cdNet.totalActMat = zeros(cdNet.nCD, iterations);
mmNet.totalMMMat = zeros(cdNet.nMN, iterations);
mmNet.totalSigMat = zeros(iterations, 1);

mmHist = zeros(1, winSize);



inDiff = diff(mmNet.input, 1, 2);

mmNet.onInd1 = find(inDiff(1,:) == 1);
mmNet.onInd2 = find(inDiff(3,:) == 1);


% start the fun!
for ts_idx = 1:iterations
    
    cdNet = cdIter_v3(cdNet, ts_idx);
    
    mmNet = mmIter_v1(mmNet, cdNet, ts_idx);
    
    mmHist = [mmHist(1:winSize-1) mmNet.totalSigMat(ts_idx)];
    
    thisMM = sum(mmHist(:));
    if thisMM >= mmNet.sigmaThresh
        cdNet.m = cdNet.m + deltaM; %desired internal synapses per cd nrn
        cdNet.Wmax = cdNet.max_wCD * cdNet.m;
        if cdNet.m < 1
            cdNet.m = 1;
            cdNet.Wmax = cdNet.max_wCD * cdNet.m;
        end
    end
    
    [cdNet, mmNet] = updateWs_v2(cdNet, mmNet, ts_idx);
    
    

    % plotting?
    if ts_idx == 1
        plotInd = flipud([1;2;3;4;5;30;34;44;57;89;27;79;87;88;120;26;69;83;96;32;47;80;21;25;6;7;8;9;10;49;73;85;107;35;98;100;39;55;70;106;28;60;78;45;101;118;11;12;13;14;15;24;56;90;42;63;119;68;113;86;108;58;112;16;17;18;19;20;36;50;71;95;48;102;116;33;53;111;23;104;67;82;22;29;31;37;38;40;41;43;46;51;52;54;59;61;62;64;65;66;72;74;75;76;77;81;84;91;92;93;94;97;99;103;105;109;110;114;115;117]);
    end
    
    if ts_idx > 0
        if (ts_idx < winSize*5) || (mod(ts_idx, winSize) == 0)
            cdWeight_plotter_v2(cdNet, ts_idx, plotInd, 5)
            pause(eps)
            
            if mod(ts_idx, winSize) == 0%if  ts_idx >= 15
                seqCD_plotter_v6(cdNet, mmNet, ts_idx, winSize, plotInd, 10); %fig10 is mismatch and activity
%                 seqCD_plotter_forShow(cdNet, mmNet, ts_idx, winSize, plotInd, 10); %fig10 is mismatch and activity
                disp(cdNet.m)
                pause(eps)
            end
        end
    end

end