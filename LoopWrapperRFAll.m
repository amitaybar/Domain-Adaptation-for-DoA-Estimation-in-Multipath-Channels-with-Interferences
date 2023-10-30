close all
clc
clear
%%
if ~exist('results','dir')
    mkdir('results');
end
%% Flags for plotting
% PlotFigs_7_8_9         % 
% PlotFigs_10            % 
% PlotFigs_11_12    	   %
% ComputeTableSirFlag  
% ComputeTableSnrFlag  

%%
SimDurVecMin = [];
%%
KFlags4Plot = 4;
%%
for CurFlagInd=1:KFlags4Plot
    CurFlagInd
    FlagsArray = zeros(1,KFlags4Plot);
    FlagsArray(CurFlagInd) = 1;
    
    PlotFigs_7_8_9          = FlagsArray(1);
    PlotFigs_11_12          = FlagsArray(2);
    ComputeTableSirFlag     = FlagsArray(3);
	ComputeTableSnrFlag     = FlagsArray(4);
    
    if PlotFigs_7_8_9 %
        KInterferencePos    = 20;%20
        KTrainPoints        = 100;%100
        KArtPoints          = 200;%200
        KTestPoints         = 300;%300
        SirMax_dB           = -20;
        SirMin_dB           = -20;
        SnrMax_dB           = 20;
        SnrMin_dB           = 20;
        IsInterferenceFixed = 1;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 1;
    end
    if PlotFigs_11_12 %
        KInterferencePos    = 1;
        KTrainPoints        = 100;%100
        KArtPoints          = 200;%200
        KTestPoints         = 300;%300
        SirMax_dB           = -10;
        SirMin_dB           = -20;
        SnrMax_dB           = 20;
        SnrMin_dB           = 20;
        IsInterferenceFixed = 0;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 2;
    end
    if ComputeTableSirFlag
        KInterferencePos    = 20; %20
        KTrainPoints        = 100;%100;
        KArtPoints          = 200;%200
        KTestPoints         = 300;%300;
        SirMax_dB           = 0;
        SirMin_dB           = -20;
        SnrMax_dB           = 20;
        SnrMin_dB           = 20;
        IsInterferenceFixed = 1;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 2;
    end
    if ComputeTableSnrFlag
        KInterferencePos    = 20; %20
        KTrainPoints        = 100;%100;
        KArtPoints          = 200;%200
        KTestPoints         = 300;%300;
        SirMax_dB           = 0;%
        SirMin_dB           = 0;%
        SnrMax_dB           = 20;
        SnrMin_dB           = 0;
        IsInterferenceFixed = 1;
        SigDim4Music_pt     = 1;
        SigDim4Music        = 2;
    end
    %%    
    MainDAforDoA_RF_MC;
    clearvars -except FlagsArray KFlags4Plot CurFlagInd SimDurVecMin
end
