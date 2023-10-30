% close all
% clear all
clc
%%
tic
%%
AngBetweenVecs          = @(u,v) ( acos(real(u/norm(u)*v'/norm(v)))*180/pi );
%%
rng(10,'twister');
%%
% addpath([pwd '\STFT']);
% addpath([pwd '\Win4Stft']);
addpath(fullfile(pwd,'STFT'));
addpath(fullfile(pwd,'Win4Stft'));
%%
linewd                  = 0.8;
hcfontsize              = 12;
MarkerSize              = 9;
%%
BlueColDef              = [0 0.4470 0.7410];
OrangeColDef            = [0.9290    0.6940    0.1250];
RedColDef               = [0.8500    0.3250    0.0980];
PurpleColDef            = [0.4940, 0.1840, 0.5560];
GreenColDef             = [0.4660, 0.6740, 0.1880];
%%
KMics               = 9;
DistBetweenMics     = 0.12;
KInterPos           = 1;
Lx                  = 4;
Ly                  = 2;
Lx_Inter            = 1;
Ly_Inter            = 2;
DistBetweenTrain    = 0.2;
TolDeg              = 0.1; %[deg] 
Epsilon             = 1e-7;
%%
c                       = 340; %[m/s]               % Sound velocity in meters / second
RoomDim                 = [5.2 6.2 3.5];            % Room dimensions in meters
fs                      = 12000;%                                         
SigLength               = 2.5*  fs;                                 
WinLength               = 2*1024;
KSampleRir              = 2048;% 
OrderTr                 = -1;
%%
TrainPosXStart          = 0.5;
TrainPosYStart          = 3.5;
%%
BinInd                  = round(1/(2*DistBetweenMics/c)/fs*WinLength+0.5);% 
FreqBinInter        	= fs/WinLength*(BinInd-0.5); %  
Lambda                  = c/FreqBinInter;
%%
PerFiltOut  = 0.03;
SamplesOut  = KInterferencePos*KTestPoints*PerFiltOut;
%% MC params
SnrMax_dB           = 20;
SnrMin_dB           = 20;
SnrVec_dB           = SnrMin_dB:10:SnrMax_dB;
SirVec_dB           = SirMin_dB:10:SirMax_dB;
%%
TrainIsTestFlag     = 1;
IsOnlyInterfActiveDuringTrainFlag = 0;
IsTwoInterferenceSourcesActive = 1 && IsInterferenceActive && ~IsInterferenceFixed; %1: two interference sources are active simulatnously ONLY in the test phase
%%
Array               = 0:KMics-1;
%%
if ~TrainIsTestFlag
    KTrainPoints = (Lx/DistBetweenTrain+1)*(Ly/DistBetweenTrain+1);
end
FirstMicPos = [ 2 0.5 ] + 0*[1.5 0 ];
DirVec = [1 0];
for m=0:KMics-1
    MicsPosMat(m+1,:) = FirstMicPos + m*DistBetweenMics*DirVec;
end
%  
MicsPosMat = [MicsPosMat 1.5*ones(length(MicsPosMat),1)];
ArrayVec = MicsPosMat(2,:) - MicsPosMat(1,:);
MidArrayPos = mean(MicsPosMat,1);

 
%%
ThetaEstTsErr_SNR           = zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));
ThetaEstTs_ptErr_SNR      	= zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));
ThetaEstMvdrTsErr_SNR       = zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));
ThetaEstMvdrTs_ptErr_SNR    = zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));
NaiveEstErr_SNR             = zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));
SirDsTsTot_SNR              = zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));
SirDsTsTott_pt_SNR          = zeros(KInterferencePos*KTestPoints,length(SnrVec_dB),length(SirVec_dB));

for int_pos = 1: KInterferencePos
    
    if TrainIsTestFlag
        TrainPosMat = [Lx*rand(KTrainPoints,1) Ly*rand(KTrainPoints,1)] + [TrainPosXStart TrainPosYStart];
    else
        TrainPosX       = (TrainPosXStart:DistBetweenTrain:TrainPosXStart+Lx);
        TrainPosY       = (TrainPosYStart:DistBetweenTrain:TrainPosYStart+Ly);
        [TrainPosXMesh, TrainPosYMesh] = meshgrid(TrainPosX,TrainPosY);
        TrainPosMat     = [TrainPosXMesh(:) TrainPosYMesh(:)];
    end
    
        ArtPosMat       = [Lx*rand(KArtPoints,1) Ly*rand(KArtPoints,1)] + [TrainPosXStart TrainPosYStart];
        ArtPosMat       = [ArtPosMat 1.5*ones(KArtPoints,1)];

    
    if IsInterferenceFixed
        InterfPosTrain = [Lx*rand(KInterPos,1) Ly*rand(KInterPos,1)] + [TrainPosXStart TrainPosYStart];
        InterfPosTest   = InterfPosTrain;
    else
        InterStartPoint = [0 1];
        InterfPosTrain  = [Lx_Inter*rand(KTrainPoints,1) Ly_Inter*rand(KTrainPoints,1)] + InterStartPoint ;
        InterfPosTest   = [Lx_Inter*rand(KTestPoints,1)  Ly_Inter*rand(KTestPoints,1)]   + InterStartPoint ;
    end
    
    %      
    TestPosMat      = [Lx*rand(KTestPoints,1) Ly*rand(KTestPoints,1)] + [TrainPosXStart TrainPosYStart];
    %  
    TrainPosMat     = [TrainPosMat 1.5*ones(length(TrainPosMat),1)];
    InterfPosTrain  = [InterfPosTrain 1.5*ones(size(InterfPosTrain,1),1)];
    InterfPosTest   = [InterfPosTest 1.5*ones(size(InterfPosTest,1),1)];
    TestPosMat      = [TestPosMat 1.5*ones(length(TestPosMat),1)];
    
    InterCenterAreatDeg = AngBetweenVecs(mean(InterfPosTrain,1)        - MidArrayPos,ArrayVec);  % for null steering
    InterCenterAreatRad = InterCenterAreatDeg/180*pi;
    
    if  (PlotFigs_2_3_4 || PlotFigs_6) && int_pos==1     
        %% 3D plot of the scenario
        %
        figure; hold on
        p1 = plot3(TrainPosMat(:,1),TrainPosMat(:,2),1.5*ones(length(TrainPosMat)),'o','Color',BlueColDef,'MarkerSize',6,'MarkerFaceColor',BlueColDef);%
        p2 = plot3(TestPosMat(:,1),TestPosMat(:,2),1.5*ones(size(TestPosMat,1)),'p','Color',OrangeColDef,'MarkerFaceColor',OrangeColDef);
        if IsInterferenceActive
            p3 = plot3(InterfPosTest(:,1),InterfPosTest(:,2),InterfPosTest(:,3),'>','Color',RedColDef,'MarkerSize',7,'MarkerFaceColor',RedColDef);
        end
        plotcube(RoomDim,[ 0  0  0],0.01,[1 0 0]);
        p4 = plot3(MicsPosMat(:,1),MicsPosMat(:,2),1.5*ones(length(MicsPosMat)),'x','MarkerSize',12,'Color',PurpleColDef);
        if IsInterferenceActive
            legend([p1(1) p2(1) p3(1) p4(1)],{'Adaptation','Operational','Interference','Array'},'Location','Best');
        else
            legend([p1(1) p2(1) p4(1)],{'Adaptation','Operational','Array'},'Location','Best');
        end
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
%         set(gca, 'FontSize', hcfontsize);
        box on;
        axis equal;
        if PlotFigs_2_3_4
            saveas(gcf,fullfile(pwd,'results', 'Fig2_AcousticRoom3D.jpg'))
            saveas(gcf,fullfile(pwd,'results', 'Fig2_AcousticRoom3D.fig'))
        else % Figure 6
            saveas(gcf,fullfile(pwd,'results', 'Fig6_AcousticRoom3D.jpg'))
            saveas(gcf,fullfile(pwd,'results', 'Fig6_AcousticRoom3D.fig'))
        end
        %%
    end
    
    for sirInd = 1:length(SirVec_dB)
        sirInd;
        Sir_LinTrain        = 1/10^(SirVec_dB(sirInd)/20);
        Sir_LinTest         = Sir_LinTrain;
        for snrInd = 1:length(SnrVec_dB)
            snrInd;
            SnrTrain    = SnrVec_dB(snrInd);
            SnrTest     = SnrTrain;
            %%
            SigTrain        =  randn(SigLength,KTrainPoints);
            SigInterTrain   =  randn(SigLength,KTrainPoints);
            SigTest         =  randn(SigLength,KTestPoints);
            SigInterTest    =  randn(SigLength,KTestPoints);
            
            for b=1:length(BetaVec)
                b
                Beta = BetaVec(b);
                for SourceInd = 1:KTrainPoints
                    [ztr(:,:,SourceInd),CorrTrVec(:,:,SourceInd)]   = SigCorrAtMicsAcousticArrayFunc(TrainPosMat(SourceInd,:),MicsPosMat,KMics,Beta,fs,BinInd,RoomDim,OrderTr,KSampleRir,WinLength,SigTrain(:,SourceInd),SnrTrain);
                    if IsInterferenceActive
                        [InterSigAtArray,~] = SigCorrAtMicsAcousticArrayFunc(InterfPosTrain(randi(size(InterfPosTrain,1)),:),MicsPosMat,KMics,Beta,fs,BinInd,RoomDim,OrderTr,KSampleRir,WinLength,SigInterTrain(:,SourceInd),SnrTrain);
                    else
                        InterSigAtArray = 0;
                    end
                    SigTmpTrain(:,:,SourceInd)                      = (IsOnlyInterfActiveDuringTrainFlag-1)*ztr(:,:,SourceInd) + Sir_LinTrain*InterSigAtArray;
                    CorrTrVec(:,:,SourceInd)                        = 1/size(SigTmpTrain(:,:,SourceInd),2)*conj(SigTmpTrain(:,:,SourceInd)'*SigTmpTrain(:,:,SourceInd));
                    CorrTrVecInv(:,:,SourceInd)                     = inv(CorrTrVec(:,:,SourceInd));
                end
                
                for SourceInd = 1:KTestPoints
                    [zts(:,:,SourceInd),CorrTsVec(:,:,SourceInd)] = SigCorrAtMicsAcousticArrayFunc(TestPosMat(SourceInd,:),MicsPosMat,KMics,Beta,fs,BinInd,RoomDim,OrderTr,KSampleRir,WinLength,SigTest(:,SourceInd),SnrTest);
                    if IsInterferenceActive
                        [InterSigAtArray,~] = SigCorrAtMicsAcousticArrayFunc(InterfPosTest(randi(size(InterfPosTest,1)),:),MicsPosMat,KMics,Beta,fs,BinInd,RoomDim,OrderTr,KSampleRir,WinLength,SigInterTest(:,SourceInd),SnrTest);
                        if IsTwoInterferenceSourcesActive
                            [InterSigAtArray2,~] = SigCorrAtMicsAcousticArrayFunc(InterfPosTest(randi(size(InterfPosTest,1)),:),MicsPosMat,KMics,Beta,fs,BinInd,RoomDim,OrderTr,KSampleRir,WinLength,SigInterTest(:,SourceInd),SnrTest);
                            InterSigAtArray = InterSigAtArray + InterSigAtArray2;
                        end
                    else
                        InterSigAtArray = 0;
                    end
                    SigTmpTest = zts(:,:,SourceInd) + Sir_LinTest*InterSigAtArray;
                    CorrTsVec(:,:,SourceInd)                    = 1/size(SigTmpTest,2)*conj(SigTmpTest'*SigTmpTest);
                end
                
                for SourceInd = 1:KArtPoints
                    for m=1:KMics
                        DistSourceMic(m)    = norm(ArtPosMat(SourceInd,:) - MicsPosMat(m,:));
                        FsAtten(m)          = 1/DistSourceMic(m);
                    end
                    SteerVecPerSource(:,SourceInd)  = FsAtten.*exp(-2j*pi*1/(Lambda)*DistSourceMic);
                    CorArtVec(:,:,SourceInd)        = SteerVecPerSource(:,SourceInd)*SteerVecPerSource(:,SourceInd)'+Epsilon*eye(KMics);
                    CorArtVecInv(:,:,SourceInd)     = inv(CorArtVec(:,:,SourceInd));
                end
                
                %%
                CorrTrMean      = mean(CorrTrVec,3);
                CorArtMean      = mean(CorArtVec,3);
                CorrTrMeanInv   = mean(CorrTrVecInv,3);
                CorArtMeanInv   = mean(CorArtVecInv,3);
                E_pt            = CorArtMean^0.5*CorrTrMean^-0.5;
                Sigma_TrHalf    = CorrTrMean^-0.5;
                Sigma_ArtHalf   = CorArtMean^0.5;%               %
                E_ptInv         = CorArtMeanInv^0.5*CorrTrMeanInv^-0.5;
                
                %%
                ThetaVec     	= linspace(0,180,1e3+1);
                for SourceInd = 1:KTestPoints
                    
                    if SourceInd<=KTrainPoints
                        Gamma1Tr = CorrTrVec(:,:,SourceInd);
                    end
                    Gamma1Ts = CorrTsVec(:,:,SourceInd);
                    for i = 1:length(ThetaVec)
                        StiringVecSpectrum                  = exp(2j*pi*Array'/(Lambda)*DistBetweenMics*cos(ThetaVec(i)/180*pi));
                        MLSpectrumOfGamma1Ts(i)             = StiringVecSpectrum'*Gamma1Ts*StiringVecSpectrum;
                        MLSpectrumOfGamma1Ts_pt(i)          = StiringVecSpectrum'*E_pt*Gamma1Ts*E_pt'*StiringVecSpectrum; % should be active
                        
                        MLSpectrumOfE_pt(i)                = abs(StiringVecSpectrum'*(E_pt)*StiringVecSpectrum)^2;
                        MLSpectrumOfE_Sigma_TrHalf(i)       = abs(StiringVecSpectrum'*(Sigma_TrHalf)*StiringVecSpectrum)^2;
                        MLSpectrumOfE_Sigma_ArtHalf(i)       = abs(StiringVecSpectrum'*(Sigma_ArtHalf)*StiringVecSpectrum)^2;
 
                        %%
                        [EigvecMat, ~]                  = SortedEVD(Gamma1Ts);
                        [EigvecMat_pt, ~]               = SortedEVD(E_pt*Gamma1Ts*E_pt');
                        UNoise                          = EigvecMat(:,SigDim4Music+1:end);
                        UNoise_pt                       = EigvecMat_pt(:,SigDim4Music_pt+1:end);
                        MusicSpectrumOfGamma1Ts(i)       = 1/((StiringVecSpectrum'*UNoise)*(StiringVecSpectrum'*UNoise)');
                        MusicSpectrumOfGamma1Ts_pt(i)    = 1/((StiringVecSpectrum'*UNoise_pt)*(StiringVecSpectrum'*UNoise_pt)');
                        MvdrSpectrumOfGamma1Ts(i)       = 1/( StiringVecSpectrum'*(inv(Gamma1Ts))*StiringVecSpectrum);
                        MvdrSpectrumOfGamma1Ts_pt(i)    = 1/(StiringVecSpectrum'*(E_ptInv*inv(Gamma1Ts)*E_ptInv')*StiringVecSpectrum);
                    end
                    ThetaTestDeg(SourceInd)                 = AngBetweenVecs(TestPosMat(SourceInd,:)        - MidArrayPos,ArrayVec);
                    if size(InterfPosTest,1) > 1
                        ThetaInterferenceDeg(SourceInd)                    = AngBetweenVecs(InterfPosTest(SourceInd,:)       - MidArrayPos,ArrayVec);
                    else
                        ThetaInterferenceDeg(SourceInd)                    = AngBetweenVecs(InterfPosTest       - MidArrayPos,ArrayVec);
                    end
                    
                    ThetaEstVecTs(SourceInd)                = DoAFromSpectrumFunc(ThetaVec,MLSpectrumOfGamma1Ts);
                    ThetaEstVecTs_pt(SourceInd)             = DoAFromSpectrumFunc(ThetaVec,MLSpectrumOfGamma1Ts_pt);
                    
                    ThetaEstMvdrVecTs(SourceInd)            = DoAFromSpectrumFunc(ThetaVec,MvdrSpectrumOfGamma1Ts);
                    ThetaEstMvdrVecTs_pt(SourceInd)         = DoAFromSpectrumFunc(ThetaVec,MvdrSpectrumOfGamma1Ts_pt);
                    
                    ThetaEstMusicVecTs(SourceInd)            = DoAFromSpectrumFunc(ThetaVec,MusicSpectrumOfGamma1Ts);
                    ThetaEstMusicVecTs_pt(SourceInd)         = DoAFromSpectrumFunc(ThetaVec,MusicSpectrumOfGamma1Ts_pt);
                    
                    ThetaInterferenceDegCurr = ThetaInterferenceDeg(SourceInd);
                    
                    SirDsTs(SourceInd)                      = SirCalcFunc(ThetaVec,ThetaTestDeg(SourceInd),ThetaInterferenceDegCurr,TolDeg,MLSpectrumOfGamma1Ts);
                    SirDsTs_pt(SourceInd)                   = SirCalcFunc(ThetaVec,ThetaTestDeg(SourceInd),ThetaInterferenceDegCurr,TolDeg,MLSpectrumOfGamma1Ts_pt);
                    SirMvdrTs(SourceInd)                    = SirCalcFunc(ThetaVec,ThetaTestDeg(SourceInd),ThetaInterferenceDegCurr,TolDeg,MvdrSpectrumOfGamma1Ts);
                    SirMvdrTs_pt(SourceInd)                 = SirCalcFunc(ThetaVec,ThetaTestDeg(SourceInd),ThetaInterferenceDegCurr,TolDeg,MvdrSpectrumOfGamma1Ts_pt);
                    SirMusicTs(SourceInd)                   = SirCalcFunc(ThetaVec,ThetaTestDeg(SourceInd),ThetaInterferenceDegCurr,TolDeg,MusicSpectrumOfGamma1Ts);
                    SirMusicTs_pt(SourceInd)                = SirCalcFunc(ThetaVec,ThetaTestDeg(SourceInd),ThetaInterferenceDegCurr,TolDeg,MusicSpectrumOfGamma1Ts_pt);
                    
                    if PlotFigs_2_3_4 && Beta == 0.4 && SourceInd == 14
                        rLimScale = [-15 0];
                        figure;
                        polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfGamma1Ts_pt)/max(abs(MLSpectrumOfGamma1Ts_pt))),'LineWidth',3);hold on
                        polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfGamma1Ts)/max(abs(MLSpectrumOfGamma1Ts))),':','LineWidth',3);                        
                        thetalim([0 180])
                        rlim(rLimScale );
                        p = polarplot(ThetaTestDeg(SourceInd)/180*pi*ones(1,2),rLimScale,'LineWidth',2);
                        p.Color = 'black';
                        if IsInterferenceActive
                            p = polarplot(ThetaInterferenceDeg(SourceInd)/180*pi*ones(1,2),rLimScale,'--','LineWidth',2);
                            p.Color = 'black';
                            legend('Our DS','DS','Desired','Interference');
                        end
                        legend('Our DS','DS','Desired');
                        ax = gca;
                        ax.FontSize = 16;
                        ax.FontName = 'Times';
                        saveas(gcf,fullfile('results', 'Fig3Left_Spectrum.jpg')) 
                        saveas(gcf,fullfile('results', 'Fig3Left_Spectrum.fig')) 
                        figure;
                        polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfE_pt)/max(abs(MLSpectrumOfE_pt))),'LineWidth',3);hold on
                        polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfE_Sigma_TrHalf)/max(abs(MLSpectrumOfE_Sigma_TrHalf))),'--','LineWidth',3);
                        polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfE_Sigma_ArtHalf)/max(abs(MLSpectrumOfE_Sigma_ArtHalf))),'-.','LineWidth',3);
                        thetalim([0 180])
                        rlim([-30 0] );
                        legend('E','\Sigma_A^{-0.5}','\Sigma_S^{0.5}');
                        ax = gca;
                        ax.FontSize = 16;
                        ax.FontName = 'Times';
                        saveas(gcf,fullfile('results', 'Fig3Right_Spectrum.jpg')) 
                        saveas(gcf,fullfile('results', 'Fig3Right_Spectrum.fig')) 
                    end
                end
                %
                ThetaEstTsErr(:,int_pos)            = abs(ThetaEstVecTs      - ThetaTestDeg);
                ThetaEstTs_ptErr(:,int_pos)         = abs(ThetaEstVecTs_pt      - ThetaTestDeg);
                ThetaEstMvdrTsErr(:,int_pos)        = abs(ThetaEstMvdrVecTs      - ThetaTestDeg);
                ThetaEstMvdrTs_ptErr(:,int_pos)     = abs(ThetaEstMvdrVecTs_pt    - ThetaTestDeg);
                ThetaEstMusicTsErr(:,int_pos)       = abs(ThetaEstMusicVecTs      - ThetaTestDeg);
                ThetaEstMusicTs_ptErr(:,int_pos)    = abs(ThetaEstMusicVecTs_pt    - ThetaTestDeg);
                NaiveEstErr(:,int_pos)              = abs(mean([max(ThetaTestDeg(:)) min(ThetaTestDeg(:))])      - ThetaTestDeg);
                
                SirDsTsTot(:,int_pos)               = SirDsTs;
                SirDsTsTot_pt(:,int_pos)            = SirDsTs_pt;
                SirMvdrTsTot(:,int_pos)             = SirMvdrTs;
                SirMvdrTsTot_pt(:,int_pos)          = SirMvdrTs_pt;
                SirMusicTsTot(:,int_pos)            = SirMusicTs;
                SirMusicTsTot_pt(:,int_pos)         = SirMusicTs_pt;
                
                Inds = (int_pos-1)*KTestPoints + 1:int_pos*KTestPoints;
                ThetaEstTsErr_SNR(Inds,snrInd,sirInd,b)         = ThetaEstTsErr(:,int_pos);
                ThetaEstTs_ptErr_SNR(Inds,snrInd,sirInd,b)      = ThetaEstTs_ptErr(:,int_pos);
                ThetaEstMvdrTsErr_SNR(Inds,snrInd,sirInd,b)     = ThetaEstMvdrTsErr(:,int_pos);
                ThetaEstMvdrTs_ptErr_SNR(Inds,snrInd,sirInd,b)  = ThetaEstMvdrTs_ptErr(:,int_pos);
                ThetaEstMusicTsErr_SNR(Inds,snrInd,sirInd,b)     = ThetaEstMusicTsErr(:,int_pos);
                ThetaEstMusicTs_ptErr_SNR(Inds,snrInd,sirInd,b)  = ThetaEstMusicTs_ptErr(:,int_pos);
                NaiveEstErr_SNR(Inds,snrInd,sirInd,b)           = NaiveEstErr(:,int_pos);
                
                SirDsTsTot_SNR(Inds,snrInd,sirInd,b)            = SirDsTsTot(:,int_pos);
                SirDsTsTott_pt_SNR(Inds,snrInd,sirInd,b)        = SirDsTsTot_pt(:,int_pos);
                SirMvdrTsTot_SNR(Inds,snrInd,sirInd,b)          = SirMvdrTsTot(:,int_pos);
                SirMvdrTsTot_pt_SNR(Inds,snrInd,sirInd,b)       = SirMvdrTsTot_pt(:,int_pos);
                SirMusicTsTot_SNR(Inds,snrInd,sirInd)          = SirMusicTsTot(:,int_pos);
                SirMusicTsTot_pt_SNR(Inds,snrInd,sirInd)       = SirMusicTsTot_pt(:,int_pos);
                
            end %b
        end % snrInd
    end % sirInd    
end % int_pos
%%
%%
FileNameLabel = 'Acoustic';
if ComputeTableFlag % saves the table in Latex code format
    Print2FileCodeV2Acoust;
end
%%
if PlotFigs_5 || PlotFigs_6
    SnrInd4Plot = 1;
    SirInd4Plot = 1;
    BetaInd = 1;
    ThetaEst_ptErr4BoxPlot = ThetaEstTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot,BetaInd);
    ThetaEstErr4BoxPlot = ThetaEstTsErr_SNR(:,SnrInd4Plot,SirInd4Plot,BetaInd);
    ThetaEstMvdr_ptErr4BoxPlot = ThetaEstMvdrTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot,BetaInd);
    ThetaEstMvdrErr4BoxPlot = ThetaEstMvdrTsErr_SNR(:,SnrInd4Plot,SirInd4Plot,BetaInd);
    ThetaEstMusic_ptErr4BoxPlot = ThetaEstMusicTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot,BetaInd);
    ThetaEstMusicErr4BoxPlot = ThetaEstMusicTsErr_SNR(:,SnrInd4Plot,SirInd4Plot,BetaInd);
    
    Data4BoxPlot            = [ThetaEst_ptErr4BoxPlot(:) ThetaEstErr4BoxPlot(:) ThetaEstMvdr_ptErr4BoxPlot(:) ThetaEstMvdrErr4BoxPlot(:) ThetaEstMusic_ptErr4BoxPlot(:) ThetaEstMusicErr4BoxPlot(:) ];
    LabelsML1               = {'Our','Standard','Our', 'Standard','Our','Standard'};
    LabelsML2               = {'DS','DS','MVDR','MVDR','MUSIC','MUSIC'};
    Labels4BoxPLot          = {LabelsML1,LabelsML2};   
    save(['Data4BoxPlot_PlotFigs5Is_mvdrDA' num2str(PlotFigs_5) '_PlotFigs6Is' num2str(PlotFigs_6)],'Data4BoxPlot')
    save(['Labels4BoxPLot_PlotFigs5Is_mvdrDA' num2str(PlotFigs_5) '_PlotFigs6Is' num2str(PlotFigs_6)],'Labels4BoxPLot')
    MyBoxSubPlot4DA_Func(Data4BoxPlot, Labels4BoxPLot, 'DS and MVDR w and w/o PP ' ,3,[2,2,2],1);
    if PlotFigs_5    
        saveas(gcf,fullfile(pwd,'results', 'Fig5_DoAErrAcoust.jpg'))
        saveas(gcf,fullfile(pwd,'results', 'Fig5_DoAErrAcoust.fig'))
    else
        saveas(gcf,fullfile(pwd,'results', 'Fig6_DoAErrAcoust.jpg'))
        saveas(gcf,fullfile(pwd,'results', 'Fig6_DoAErrAcoust.fig'))
    end
end
%%
MarkerVec = ['p','o','<','s','d'];
AlphaVec = [1 0.8 0.6 0.4 0.2];
%% Beta graph for median
% mean - median
if PlotFigs_2_3_4
    SnrInd4Plot = 1;
    SirInd4Plot = 1;
    figure; hold on
    method = 1;    
    Vec2Plot = squeeze(median(ThetaEstTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot,:),1));
    s =  plot(BetaVec,Vec2Plot,'Color',BlueColDef,'LineWidth',3-0.4*method,'Marker',MarkerVec(method),'MarkerSize',6,'MarkerFaceColor',BlueColDef);
    s.Color(4) = AlphaVec(1);
    Vec2Plot = squeeze(median(ThetaEstTsErr_SNR(:,SnrInd4Plot,SirInd4Plot,:),1));
    s =  plot(BetaVec,Vec2Plot,'Color',RedColDef,'LineWidth',3-0.4*method,'Marker',MarkerVec(method),'MarkerSize',6,'MarkerFaceColor',RedColDef);
    s.Color(4) = AlphaVec(1);
    ylim([0 10]);
    legend('Our DS','DS','Location','best');
    xlabel('\beta [sec]');
    ylabel('Error[deg]');
%     set(gca, 'FontSize', hcfontsize/1.5);
    set(gca, 'LineWidth', linewd);
    box on ; grid on
    saveas(gcf,fullfile('results', 'Fig4Top_DoABetaDS.jpg')) 
    saveas(gcf,fullfile('results', 'Fig4Top_DoABetaDS.fig')) 
    figure; hold on
    method = 2;    
    Vec2Plot = squeeze(median(ThetaEstMvdrTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot,:),1));
    s =  plot(BetaVec,Vec2Plot,'Color',BlueColDef,'LineWidth',3-0.4*method,'Marker',MarkerVec(method),'MarkerSize',6,'MarkerFaceColor',BlueColDef);
    s.Color(4) = AlphaVec(1);
    Vec2Plot = squeeze(median(ThetaEstMvdrTsErr_SNR(:,SnrInd4Plot,SirInd4Plot,:),1));
    s =  plot(BetaVec,Vec2Plot,'Color',RedColDef,'LineWidth',3-0.4*method,'Marker',MarkerVec(method),'MarkerSize',6,'MarkerFaceColor',RedColDef);
    s.Color(4) = AlphaVec(1);
    ylim([0 10]);
    legend('Our MVDR','MVDR','Location','best');
    xlabel('\beta [sec]');
    ylabel('Error[deg]');
%     set(gca, 'FontSize', hcfontsize/1.5);
    set(gca, 'LineWidth', linewd);
    box on ; grid on
    method = 3;
    saveas(gcf,fullfile('results', 'Fig4Mid_DoABetaMvdr.jpg')) 
    saveas(gcf,fullfile('results', 'Fig4Mid_DoABetaMvdr.fig'))
    figure; hold on    
    Vec2Plot = squeeze(median(ThetaEstMusicTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot,:),1));
    s =  plot(BetaVec,Vec2Plot,'Color',BlueColDef,'LineWidth',3-0.4*method,'Marker',MarkerVec(method),'MarkerSize',6,'MarkerFaceColor',BlueColDef);
    s.Color(4) = AlphaVec(1);
    Vec2Plot = squeeze(median(ThetaEstMusicTsErr_SNR(:,SnrInd4Plot,SirInd4Plot,:),1));
    s =  plot(BetaVec,Vec2Plot,'Color',RedColDef,'LineWidth',3-0.4*method,'Marker',MarkerVec(method),'MarkerSize',6,'MarkerFaceColor',RedColDef);
    s.Color(4) = AlphaVec(1);
    ylim([0 10]);
    legend('Our MUSIC','MUSIC','Location','best');
    xlabel('\beta [sec]');
    ylabel('Error[deg]');
%     set(gca, 'FontSize', hcfontsize/1.5);
    set(gca, 'LineWidth', linewd);
    box on ; grid on
    saveas(gcf,fullfile('results', 'Fig4Bottom_DoABetaMusic.jpg')) 
    saveas(gcf,fullfile('results', 'Fig4Bottom_DoABetaMusic.fig'))

end

%%
SimDurMin = toc/60;
SimDurVecMin(end+1) = SimDurMin;
%%
function Sir = SirCalcFunc(ThetaVec,ThetaTarDeg,ThetaTarInter,TolDeg,Spectrum)

IndSpecTar = find(abs(ThetaVec-ThetaTarDeg)<TolDeg);
IndSpecInter = find(abs(ThetaVec-ThetaTarInter)<TolDeg);
SpectrumNorm = Spectrum/max(abs(Spectrum));
Sir = 10*log10(max(abs(SpectrumNorm(IndSpecTar)))) -  10*log10(max(abs(SpectrumNorm(IndSpecInter))));

end

function plotcube(varargin)
%
%

%
inArgs = { ...
    [10 56 100] , ... %
    [10 10  10] , ... %
    .7          , ... %
    [1 0 0]       ... %
    };
%
inArgs(1:nargin) = varargin;
%
[edges,origin,alpha,clr] = deal(inArgs{:});
XYZ = { ...
    [0 0 0 0]  [0 0 1 1]  [0 1 1 0] ; ...
    [1 1 1 1]  [0 0 1 1]  [0 1 1 0] ; ...
    [0 1 1 0]  [0 0 0 0]  [0 0 1 1] ; ...
    [0 1 1 0]  [1 1 1 1]  [0 0 1 1] ; ...
    [0 1 1 0]  [0 0 1 1]  [0 0 0 0] ; ...
    [0 1 1 0]  [0 0 1 1]  [1 1 1 1]   ...
    };
XYZ = mat2cell(...
    cellfun( @(x,y,z) x*y+z , ...
    XYZ , ...
    repmat(mat2cell(edges,1,[1 1 1]),6,1) , ...
    repmat(mat2cell(origin,1,[1 1 1]),6,1) , ...
    'UniformOutput',false), ...
    6,[1 1 1]);
cellfun(@patch,XYZ{1},XYZ{2},XYZ{3},...
    repmat({clr},6,1),...
    repmat({'FaceAlpha'},6,1),...
    repmat({alpha},6,1)...
    );
view(3);
end