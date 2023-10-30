close all
clc
clear
%
%%
if ~exist('results','dir')
    mkdir('results');
end
%% Flags for plotting
% PlotFigs_2_3_4           
% PlotFigs_5          
% PlotFigs_6           
% ComputeTableFlag     
%%
SimDurVecMin = [];
%%
KFlags4Plot = 4;
%%
for CurFlagInd=1:KFlags4Plot
    CurFlagInd
    FlagsArray = zeros(1,KFlags4Plot);
    FlagsArray(CurFlagInd) = 1;
    
    PlotFigs_2_3_4      = FlagsArray(1); %
    PlotFigs_5          = FlagsArray(2);%
    PlotFigs_6          = FlagsArray(3);%
    ComputeTableFlag    = FlagsArray(4); %
    
    %%
    if ComputeTableFlag
        KInterferencePos    = 20;
        KTrainPoints        = 100;%100
        KArtPoints          = 200;%200
        KTestPoints         = 200;%200
        BetaVec             = 0.2;%
        SirMax_dB           = 0;
        SirMin_dB           = -20;
        IsInterferenceActive = 1;
        IsInterferenceFixed = 1;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 2;
    end
    if PlotFigs_2_3_4 %
        KInterferencePos    = 1;
        KTrainPoints        = 100;%100
        KArtPoints          = 200;%200
        KTestPoints         = 300;%300
        BetaVec             = [0.2 0.3 0.4 0.5 0.6];
        SirMax_dB           = -20;
        SirMin_dB           = -20;
        IsInterferenceActive = 0;
        IsInterferenceFixed = 1;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 1;
    end
    if PlotFigs_5 % 
        KInterferencePos    = 20;
        KTrainPoints        = 100;%100
        KArtPoints          = 200;%200
        KTestPoints         = 200;%200
        BetaVec             = 0.2;
        SirMax_dB           = -20;
        SirMin_dB           = -20;
        IsInterferenceActive = 1;
        IsInterferenceFixed = 1;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 2;
    end
    if PlotFigs_6 %
        KInterferencePos    = 1;
        KTrainPoints        = 100;%100
        KArtPoints          = 200;%200
        KTestPoints         = 300;%300
        BetaVec             = 0.2;
        SirMax_dB           = 0;
        SirMin_dB           = 0;
        IsInterferenceActive = 1;
        IsInterferenceFixed = 0;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 3;
    end
    %%
    MainDAforDoA_Acoustic_MC;
    clearvars -except FlagsArray KFlags4Plot CurFlagInd SimDurVecMin
end
