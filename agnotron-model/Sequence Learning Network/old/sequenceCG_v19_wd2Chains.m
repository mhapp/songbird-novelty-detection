%% v14 wants two chains for odorants A1 B1, A2 B2
% very close...
% i need to 1) create separate variables for 1 and 2, to make the temporal
% jitter easier and 2) make the temporal jitter more extreme and varied

% now change the learning rule!

clear
close all
figure(5)
figure(10)
pause

%%% configs

seed =2509;% 23; %2509; %23;

winSize = 40;

cdNet.nCD = 120;
cdNet.nMN = 20;
cdNet.nPD = 4;
cdNet.gamma = 0;%0.045;%045; %0.02 %recurrent inhibition
cdNet.beta = 0.025;%0.025;15;%0.025;%15; %0.115;  %feed forward inhibition
cdNet.epsilon = 0.2; %relative strength of heterosynaptic LTD
cdNet.tau = 30; %time constant of adaptation
cdNet.alpha = 50;
cdNet.max_wCD = 1.5; %max weight between cd neurons
cdNet.m = 25; %10-16 ideal %desired internal synapses per cd nrn
cdNet.Wmax = cdNet.max_wCD * cdNet.m;
cdNet.nSeed = 5; % number of seed neurons per pd output
cdNet.pin = 0; %.01; 

seedThresh = 0;
seedDrive = 1;

cdNet.fCoef = 0; 
cdNet.vCoef = 1;% - cdNet.fCoef;


cdNet.eta_cd = 0.25;%15; %0.001; %25;%25;%0.5; %initial learning rate
deltaM = -0.05;%.025;% -0.05;% -0.5; %change in desired # synapses per mismatch


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
%         plotInd = [30;53;57;66;71;72;86;87;100;112;22;24;35;50;108;62;64;29;77;28;81;96;97;33;51;105;120;25;42;47;60;109;16;17;18;19;20;45;67;79;92;52;104;107;110;117;76;88;99;119;23;55;84;93;106;27;31;54;68;118;11;12;13;14;15;85;94;114;115;21;34;44;26;32;43;59;101;75;91;113;48;65;70;80;83;103;111;6;7;8;9;10;36;37;39;56;90;38;40;78;102;74;82;89;116;41;46;58;69;95;98;49;61;63;73;1;2;3;4;5];
        plotInd = flipud([1;2;3;4;5;30;34;44;57;89;27;79;87;88;120;26;84;96;32;47;91;33;40;63;6;7;8;9;10;49;73;85;107;35;46;72;98;100;39;55;70;106;28;74;78;22;66;82;41;103;108;101;11;12;13;14;15;24;56;90;31;42;119;25;68;67;86;54;102;16;17;18;19;20;36;50;71;95;48;76;116;53;77;23;109;29;75;21;37;38;43;45;51;52;58;59;60;61;62;64;65;69;80;81;83;92;93;94;97;99;104;105;110;111;112;113;114;115;117;118]);
    end
    
    if ts_idx > 0
        if (ts_idx < winSize*5) || (mod(ts_idx, winSize) == 0)
            cdWeight_plotter_v2(cdNet, ts_idx, plotInd, 5)
            pause(0.005)
            
            if mod(ts_idx, winSize) == 0%if  ts_idx >= 15
                seqCD_plotter_v6(cdNet, mmNet, ts_idx, winSize, plotInd, 10); %fig10 is mismatch and activity
                disp(cdNet.m)
                pause(0.001)
            end
        end
    end

end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           