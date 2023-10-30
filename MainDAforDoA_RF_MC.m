% close all
% clear all
% clc
%%
tic
%%
AngBetweenVecs          = @(u,v) ( acos(real(u/norm(u)*v'/norm(v)))*180/pi );
%%
rng(11,'twister');
%%
linewd                  = 0.8;
hcfontsize              = 12;
MarkerSize              = 9;
%%
BlueColDef              = [0 0.4470 0.7410];
OrangeColDef            = [0.9290    0.6940    0.1250];
RedColDef               = [0.8500    0.3250    0.0980];
PurpleColDef            = [0.4940, 0.1840, 0.5560];
%%
KMics               = 9;
c                   = 3e8; %[m/s]
Freq                = 2.4e9; %[Hz]
Lambda              = c/Freq;
DistBetweenMics     = Lambda/2;
SigLength           = 1e3; %[samples]
KInterPos           = 1;
Epsilon             = 1e-7; % 
Lx                  = 160;
Ly                  = 200;
Lx_Inter            = 2;
Ly_Inter            = 2;
DistBetweenTrain    = 20;
SirTrain_dB         = 0;
SirTest_dB          = 0;
Sir_LinTrain        = 1/10^(SirTrain_dB/20);
Sir_LinTest         = 1/10^(SirTest_dB/20);
TolDeg              = 0.1; %[deg] 

%% MC params
SnrVec_dB           = SnrMin_dB:10:SnrMax_dB;
SirVec_dB           = SirMin_dB:10:SirMax_dB;

%% Flags
TrainIsTestFlag     = 1;
IsOnlyInterfActiveDuringTrainFlag = 0;
IsTwoInterferenceSourcesActive = 0 && ~IsInterferenceFixed; %1: two interference sources are active simulatnously ONLY in the test phase
%% Generting sources 
RMin = 100; %[m]
RMax = 250; %[m]
if IsInterferenceFixed
    ThetaMin = 30; %[deg]
    ThetaMax = 150; %[deg]
    RMinInterf = 100; %[m]
    RMaxInterf = 250; %[m]
else
    ThetaMin = 70; %[deg]
    ThetaMax = 120; %[deg]
    RMinInterf = 100; %[m]
    RMaxInterf = 250; %[m]
    ThetaMinInterf = 130; %[deg]
    ThetaMaxInterf = 160; %[deg]
    ThetaMinInterf2 = 30;
end
%%
Array               = 0:KMics-1;
ErrPhaseVec         = 0*(360*rand(KMics,1)-180);
%%
if ~TrainIsTestFlag
    KTrainPoints = (Lx/DistBetweenTrain+1)*(Ly/DistBetweenTrain+1);
end

FirstMicPos = [ 0 0 ];
DirVec      = [rand  0.3*rand];
DirVec      = DirVec / norm(DirVec);
DirVec = [1 0];
for m=0:KMics-1
    MicsPosMat(m+1,:) = FirstMicPos + m*DistBetweenMics*DirVec;
end

ArrayVec = MicsPosMat(2,:)-MicsPosMat(1,:);
MidArrayPos = mean(MicsPosMat,1);

TrainPosXStart = FirstMicPos(1)-Lx/2;
TrainPosYStart = FirstMicPos(1)+ 120;


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
        
        RadiusTr     = (RMax-RMin)*rand(KTrainPoints,1) + RMin;
        ThetaTr      = (ThetaMax-ThetaMin)*rand(KTrainPoints,1)/180*pi + ThetaMin/180*pi;
        TrainPosMat  = [RadiusTr.*cos(ThetaTr) RadiusTr.*sin(ThetaTr)];
    else
        TrainPosX       = (TrainPosXStart:DistBetweenTrain:TrainPosXStart+Lx);
        TrainPosY       = (TrainPosYStart:DistBetweenTrain:TrainPosYStart+Ly);
        [TrainPosXMesh, TrainPosYMesh] = meshgrid(TrainPosX,TrainPosY);
        TrainPosMat     = [TrainPosXMesh(:) TrainPosYMesh(:)];
    end
    
    
    ArtPosMat       = [Lx*rand(KArtPoints,1) Ly*rand(KArtPoints,1)] + [TrainPosXStart TrainPosYStart];
    
    RadiusArt       = (RMax-RMin)*rand(KArtPoints,1) + RMin;
    ThetaArt        = (ThetaMax-ThetaMin)*rand(KArtPoints,1)/180*pi + ThetaMin/180*pi;
    ArtPosMat       = [RadiusArt.*cos(ThetaArt) RadiusArt.*sin(ThetaArt)];
    
    if IsInterferenceFixed
        InterfPosTrain = [Lx*rand(KInterPos,1) Ly*rand(KInterPos,1)] + [TrainPosXStart TrainPosYStart];
        RadiusInterf   	= (RMax-RMin)*rand + RMin;
        ThetaInterf   	= (ThetaMax-ThetaMin)*rand/180*pi + ThetaMin/180*pi;
        InterfPosTrain 	= [RadiusInterf.*cos(ThetaInterf) RadiusInterf.*sin(ThetaInterf)];
        InterfPosTest   = InterfPosTrain;
        InterfPosTest 	= repmat(InterfPosTest,KTestPoints,1);
    else
        pTrain=randi(2,KTrainPoints,1)-1;
        pTest =randi(2,KTestPoints,1)-1;
        RadiusInterf   	= (RMaxInterf-RMinInterf)*rand(KTrainPoints,1) + RMinInterf;
        ThetaInterf   	= (ThetaMaxInterf-ThetaMinInterf)*rand(KTrainPoints,1)/180*pi + pTrain*ThetaMinInterf/180*pi + (1-pTrain)*ThetaMinInterf2/180*pi;
        InterfPosTrain 	= [RadiusInterf.*cos(ThetaInterf) RadiusInterf.*sin(ThetaInterf)];
        RadiusInterf   	= (RMaxInterf-RMinInterf)*rand(KTestPoints,1) + RMinInterf;
        ThetaInterf   	= (ThetaMaxInterf-ThetaMinInterf)*rand(KTestPoints,1)/180*pi + pTest*ThetaMinInterf/180*pi  + (1-pTest)*ThetaMinInterf2/180*pi;
        InterfPosTest 	= [RadiusInterf.*cos(ThetaInterf) RadiusInterf.*sin(ThetaInterf)];
    end
    
    TestPosMat      = [Lx*rand(KTestPoints,1) Ly*rand(KTestPoints,1)] + [TrainPosXStart TrainPosYStart];
    RadiusTs     = (RMax-RMin)*rand(KTestPoints,1) + RMin;
    ThetaTs      = (ThetaMax-ThetaMin)*rand(KTestPoints,1)/180*pi + ThetaMin/180*pi;
    TestPosMat  = [RadiusTs.*cos(ThetaTs) RadiusTs.*sin(ThetaTs)];
    
    InterCenterAreatDeg = AngBetweenVecs(mean(InterfPosTrain,1)        - MidArrayPos,ArrayVec);  % for null steering
    InterCenterAreatRad = InterCenterAreatDeg/180*pi;
    
    if  (PlotFigs_7_8_9 || PlotFigs_11_12) && int_pos==1
        figure; hold on
        plot(TrainPosMat(:,1),TrainPosMat(:,2),'o','Color',BlueColDef,'MarkerSize',6,'MarkerFaceColor',BlueColDef);%
        % plot(ArtPosMat(:,1),ArtPosMat(:,2),'x');
        plot(TestPosMat(:,1),TestPosMat(:,2),'p','Color',OrangeColDef,'MarkerFaceColor',OrangeColDef);
        plot(InterfPosTest(:,1),InterfPosTest(:,2),'>','Color',RedColDef,'MarkerSize',10,'MarkerFaceColor',RedColDef);
        plot(MicsPosMat(:,1),MicsPosMat(:,2),'x','MarkerSize',12,'Color',PurpleColDef);
        h=legend('Adaptation','Operational','Interference','Array','Location','Best');
        h.FontSize = 12;
        % plot(InterferencePos(1),InterferencePos(2),'s');
        xlabel('[m]');
        ylabel('[m]');
        %         set(gca, 'FontSize', hcfontsize);
        % set(gca, 'LineWidth', linewd);
        box on; grid on
        axis equal
        if PlotFigs_7_8_9
            saveas(gcf,fullfile('results', 'Fig7_RFSetting.jpg'))
            saveas(gcf,fullfile('results', 'Fig7_RFSetting.fig'))
        end
        if PlotFigs_11_12
            saveas(gcf,fullfile('results', 'Fig11_RFSetting2Secs.jpg'))
            saveas(gcf,fullfile('results', 'Fig11_RFSetting2Secs.fig'))
        end
    end
    
    for sirInd = 1:length(SirVec_dB)
        SirVec_dB(sirInd)
        Sir_LinTrain        = 1/10^(SirVec_dB(sirInd)/20);
        Sir_LinTest         = Sir_LinTrain;
        for snrInd = 1:length(SnrVec_dB)
            SnrVec_dB(snrInd)
            SnrTrain    = SnrVec_dB(snrInd);
            SnrTest     = SnrTrain;
            %%
            SigTrain        =  randn(SigLength,KTrainPoints)+1j*randn(SigLength,KTrainPoints);
            SigInterTrain   =  randn(SigLength,KTrainPoints)+1j*randn(SigLength,KTrainPoints);
            %  
            SigTest         =  randn(SigLength,KTestPoints)+1j*randn(SigLength,KTestPoints);
            SigInterTest    =  randn(SigLength,KTestPoints)+1j*randn(SigLength,KTestPoints);
            
            [InterSigAtArray,~] = SigCorrAtMicsArrayFunc(InterfPosTrain(randi(size(InterfPosTrain,1)),:),MicsPosMat,KMics,Lambda,SigLength,SigTrain(:,1),SnrTrain);
            
            for SourceInd = 1:KTrainPoints
                
                [ztr(:,:,SourceInd),CorrTrVec(:,:,SourceInd)]   = SigCorrAtMicsArrayFunc(TrainPosMat(SourceInd,:),MicsPosMat,KMics,Lambda,SigLength,SigTrain(:,SourceInd),SnrTrain);
                
                [InterSigAtArray,~] = SigCorrAtMicsArrayFunc(InterfPosTrain(randi(size(InterfPosTrain,1)),:),MicsPosMat,KMics,Lambda,SigLength,SigInterTrain(:,SourceInd),SnrTrain);
                
                SigTmpTrain(:,:,SourceInd)                      = (IsOnlyInterfActiveDuringTrainFlag-1)*ztr(:,:,SourceInd) + Sir_LinTrain*InterSigAtArray;
                CorrTrVec(:,:,SourceInd)                        = 1/size(SigTmpTrain(:,:,SourceInd),2)*conj(SigTmpTrain(:,:,SourceInd)'*SigTmpTrain(:,:,SourceInd));
                CorrTrVecInv(:,:,SourceInd)                     = inv(CorrTrVec(:,:,SourceInd));
                
            end
            SigTmpTrainSum = mean(SigTmpTrain,3); %  
            CorrTrVecSingleSig = conj(SigTmpTrainSum'*SigTmpTrainSum);
            
            for SourceInd = 1:KTestPoints
                [zts(:,:,SourceInd),CorrTsVec(:,:,SourceInd)]   = SigCorrAtMicsArrayFunc(TestPosMat(SourceInd,:),MicsPosMat,KMics,Lambda,SigLength,SigTest(:,SourceInd),SnrTest);
                [InterSigAtArray,~]                             = SigCorrAtMicsArrayFunc(InterfPosTest(SourceInd,:),MicsPosMat,KMics,Lambda,SigLength,SigInterTest(:,SourceInd),SnrTrain);
                if IsTwoInterferenceSourcesActive
                    [InterSigAtArray2,~]                        = SigCorrAtMicsArrayFunc(InterfPosTest(randi(size(InterfPosTest,1)),:),MicsPosMat,KMics,Lambda,SigLength,SigInterTest(:,SourceInd),SnrTrain);
                    InterSigAtArray                             = InterSigAtArray + InterSigAtArray2;
                end
                SigTmpTest                                      = zts(:,:,SourceInd) + Sir_LinTest*InterSigAtArray;
                CorrTsVec(:,:,SourceInd)                        = 1/size(SigTmpTest,2)*conj(SigTmpTest'*SigTmpTest);
            end
            
            for SourceInd = 1:KArtPoints
                for m=1:KMics
                    DistSourceMic(m)            = norm(ArtPosMat(SourceInd,:) - MicsPosMat(m,:));
                    FsAtten(m)                  = 1/DistSourceMic(m);
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
            Sigma_ArtHalf   = CorArtMean^0.5;
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
                    MLSpectrumOfGamma1Tr(i)             = StiringVecSpectrum'*Gamma1Tr*StiringVecSpectrum;
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
                    MvdrSpectrumOfGamma1Ts_pt(i)    = 1/(StiringVecSpectrum'*inv(E_pt*(Gamma1Ts)*E_pt')*StiringVecSpectrum); %
                end
                
                ThetaTestDeg(SourceInd)                 = AngBetweenVecs(TestPosMat(SourceInd,:)        - MidArrayPos,ArrayVec);
                if size(InterfPosTest,1) > 1
                    ThetaInterferenceDeg(SourceInd)         = AngBetweenVecs(InterfPosTest(SourceInd,:)        - MidArrayPos,ArrayVec);
                else
                    ThetaInterferenceDeg(SourceInd)         = AngBetweenVecs(InterfPosTest        - MidArrayPos,ArrayVec);
                end
                
                ThetaEstVecTr(SourceInd)                = DoAFromSpectrumFunc(ThetaVec,MLSpectrumOfGamma1Tr);
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
                
                if PlotFigs_7_8_9 && SourceInd ==3 && int_pos ==1 % 
                    rLimScale = [-15 0];
                    figure;
                    polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfGamma1Ts_pt)/max(abs(MLSpectrumOfGamma1Ts_pt))),'LineWidth',3);hold on
                    polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfGamma1Ts)/max(abs(MLSpectrumOfGamma1Ts))),':','LineWidth',3);
                    thetalim([0 180])
                    rlim(rLimScale );
                    p = polarplot(ThetaTestDeg(SourceInd)/180*pi*ones(1,2),rLimScale,'LineWidth',2);
                    p.Color = 'black';
                    p = polarplot(ThetaInterferenceDeg(SourceInd)/180*pi*ones(1,2),rLimScale,'--','LineWidth',2);
                    p.Color = 'black';
                    legend('Our DS','DS','Desired','Interference');
                    ax = gca;
                    ax.FontSize = 18;
                    ax.FontName = 'Times';
                    saveas(gcf,fullfile('results', 'Fig8Left_Spectrum.jpg'))
                    saveas(gcf,fullfile('results', 'Fig8Left_Spectrum.fig'))
                    
                    figure;
                    polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfE_pt)/max(abs(MLSpectrumOfE_pt))),'LineWidth',3);hold on
                    polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfE_Sigma_TrHalf)/max(abs(MLSpectrumOfE_Sigma_TrHalf))),'--','LineWidth',3);
                    polarplot(ThetaVec/180*pi,10*log10(abs(MLSpectrumOfE_Sigma_ArtHalf)/max(abs(MLSpectrumOfE_Sigma_ArtHalf))),'-.','LineWidth',3);
                    p = polarplot(ThetaInterferenceDeg(SourceInd)/180*pi*ones(1,2),[-30 0],'--','LineWidth',2);
                    p.Color = 'black';
                    thetalim([0 180])
                    rlim([-40 0] );
                    legend('E','\Sigma_A^{-0.5}','\Sigma_S^{0.5}');
                    ax = gca;
                    ax.FontSize = 18;
                    ax.FontName = 'Times';
                    saveas(gcf,fullfile('results', 'Fig8Right_Spectrum.jpg'))
                    saveas(gcf,fullfile('results', 'Fig8Right_Spectrum.fig'))
                end
            end
            ThetaEstTsErr(:,int_pos)         = abs(ThetaEstVecTs      - ThetaTestDeg);
            ThetaEstTs_ptErr(:,int_pos)      = abs(ThetaEstVecTs_pt      - ThetaTestDeg);
            ThetaEstMvdrTsErr(:,int_pos)     = abs(ThetaEstMvdrVecTs      - ThetaTestDeg);
            ThetaEstMvdrTs_ptErr(:,int_pos)  = abs(ThetaEstMvdrVecTs_pt    - ThetaTestDeg);
            ThetaEstMusicTsErr(:,int_pos)     = abs(ThetaEstMusicVecTs      - ThetaTestDeg);
            ThetaEstMusicTs_ptErr(:,int_pos)  = abs(ThetaEstMusicVecTs_pt    - ThetaTestDeg);
            NaiveEstErr(:,int_pos)           = abs(mean([max(ThetaTestDeg(:)) min(ThetaTestDeg(:))])      - ThetaTestDeg);
            
            SirDsTsTot(:,int_pos)           = SirDsTs;
            SirDsTsTot_pt(:,int_pos)        = SirDsTs_pt;
            SirMvdrTsTot(:,int_pos)         = SirMvdrTs;
            SirMvdrTsTot_pt(:,int_pos)      = SirMvdrTs_pt;
            SirMusicTsTot(:,int_pos)         = SirMusicTs;
            SirMusicTsTot_pt(:,int_pos)      = SirMusicTs_pt;
            
            
            
            Inds = (int_pos-1)*KTestPoints + 1:int_pos*KTestPoints;
            ThetaEstTsErr_SNR(Inds,snrInd,sirInd)         = ThetaEstTsErr(:,int_pos);
            ThetaEstTs_ptErr_SNR(Inds,snrInd,sirInd)      = ThetaEstTs_ptErr(:,int_pos);
            ThetaEstMvdrTsErr_SNR(Inds,snrInd,sirInd)     = ThetaEstMvdrTsErr(:,int_pos);
            ThetaEstMvdrTs_ptErr_SNR(Inds,snrInd,sirInd)  = ThetaEstMvdrTs_ptErr(:,int_pos);
            ThetaEstMusicTsErr_SNR(Inds,snrInd,sirInd)     = ThetaEstMusicTsErr(:,int_pos);
            ThetaEstMusicTs_ptErr_SNR(Inds,snrInd,sirInd)  = ThetaEstMusicTs_ptErr(:,int_pos);
            NaiveEstErr_SNR(Inds,snrInd,sirInd)           = NaiveEstErr(:,int_pos);
            
            SirDsTsTot_SNR(Inds,snrInd,sirInd)            = SirDsTsTot(:,int_pos);
            SirDsTsTott_pt_SNR(Inds,snrInd,sirInd)        = SirDsTsTot_pt(:,int_pos);
            SirMvdrTsTot_SNR(Inds,snrInd,sirInd)          = SirMvdrTsTot(:,int_pos);
            SirMvdrTsTot_pt_SNR(Inds,snrInd,sirInd)       = SirMvdrTsTot_pt(:,int_pos);
            SirMusicTsTot_SNR(Inds,snrInd,sirInd)          = SirMusicTsTot(:,int_pos);
            SirMusicTsTot_pt_SNR(Inds,snrInd,sirInd)       = SirMusicTsTot_pt(:,int_pos);
            
        end % snrInd
    end % sirInd
    
end % int_pos
%%
FileNameLabel = 'RF';
if ComputeTableSirFlag
    FileNameLabel = 'RFSirTable';
    Print2FileCodeV2FixedSir_SIR;
end
if ComputeTableSnrFlag
    FileNameLabel = 'RFSnrTable';
    Print2FileCodeV2FixedSir_SNR;
end
%%
if PlotFigs_7_8_9
    SnrVec_dB;
    SirVec_dB;
    SnrInd4Plot = 1;
    SirInd4Plot = 1;
    ThetaEst_ptErr4BoxPlot = ThetaEstTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot);
    ThetaEstErr4BoxPlot = ThetaEstTsErr_SNR(:,SnrInd4Plot,SirInd4Plot);
    ThetaEstMvdr_ptErr4BoxPlot = ThetaEstMvdrTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot);
    ThetaEstMvdrErr4BoxPlot = ThetaEstMvdrTsErr_SNR(:,SnrInd4Plot,SirInd4Plot);
    ThetaEstMusic_ptErr4BoxPlot = ThetaEstMusicTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot);
    ThetaEstMusicErr4BoxPlot = ThetaEstMusicTsErr_SNR(:,SnrInd4Plot,SirInd4Plot);
    
    Data4BoxPlot            = [ThetaEst_ptErr4BoxPlot(:) ThetaEstErr4BoxPlot(:) ThetaEstMvdr_ptErr4BoxPlot(:) ThetaEstMvdrErr4BoxPlot(:) ThetaEstMusic_ptErr4BoxPlot(:) ThetaEstMusicErr4BoxPlot(:) ];
    
    LabelsML1               = {'Our','Standard','Our', 'Standard','Our','Standard'};
    LabelsML2               = {'DS','DS','MVDR', 'MVDR','MUSIC','MUSIC'};%LabelsML2               = {'','','','','',''};
    Labels4BoxPLot          = {LabelsML1,LabelsML2};
    
    save('Data4BoxPlot_PlotFigs_7_8_9Is1','Data4BoxPlot');
    save('Labels4BoxPLot_PlotFigs_7_8_9Is1','Labels4BoxPLot');
    MyBoxSubPlot4DA_Func(Data4BoxPlot, Labels4BoxPLot, 'DS and MVDR w and w/o PP ' ,3,[2,2,2],1);
    
    saveas(gcf,fullfile(pwd,'results', 'Fig9_DoAErrRF.jpg'))
    saveas(gcf,fullfile(pwd,'results', 'Fig9_DoAErrRF.fig'))
end
% ylim([0 10]);

if PlotFigs_11_12
    SnrInd4Plot = 1;
    for SirInd4Plot = 1:2
        ThetaEst_ptErr4BoxPlot = ThetaEstTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot);
        ThetaEstErr4BoxPlot = ThetaEstTsErr_SNR(:,SnrInd4Plot,SirInd4Plot);
        ThetaEstMvdr_ptErr4BoxPlot = ThetaEstMvdrTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot);
        ThetaEstMvdrErr4BoxPlot = ThetaEstMvdrTsErr_SNR(:,SnrInd4Plot,SirInd4Plot);
        ThetaEstMusic_ptErr4BoxPlot = ThetaEstMusicTs_ptErr_SNR(:,SnrInd4Plot,SirInd4Plot);
        ThetaEstMusicErr4BoxPlot = ThetaEstMusicTsErr_SNR(:,SnrInd4Plot,SirInd4Plot);
        
        Data4BoxPlot            = [ThetaEst_ptErr4BoxPlot(:) ThetaEstErr4BoxPlot(:) ThetaEstMvdr_ptErr4BoxPlot(:) ThetaEstMvdrErr4BoxPlot(:) ThetaEstMusic_ptErr4BoxPlot(:) ThetaEstMusicErr4BoxPlot(:) ];
        
        LabelsML1               = {'Our','Standard','Our', 'Standard','Our','Standard'};
        LabelsML2               = {'DS','DS','MVDR', 'MVDR','MUSIC','MUSIC'};
        Labels4BoxPLot          = {LabelsML1,LabelsML2};
        save(['Data4BoxPlot_PlotFigs_11_12Is1_Ind' num2str(SirInd4Plot)],'Data4BoxPlot');
        save(['Labels4BoxPLot_PlotFigs_11_12Is1_Ind'  num2str(SirInd4Plot)],'Labels4BoxPLot');
        MyBoxSubPlot4DA_Func(Data4BoxPlot, Labels4BoxPLot, 'DS and MVDR w and w/o PP ' ,3,[2,2,2],1);
        saveas(gcf,fullfile(pwd,'results', ['Fig12_DoAErrRFSir_' num2str(abs(SirVec_dB(SirInd4Plot))) '.jpg']))
        saveas(gcf,fullfile(pwd,'results', ['Fig12_DoAErrRFSir_' num2str(abs(SirVec_dB(SirInd4Plot))) '.fig']))
    end
end


%%
MarkerVec = ['p','o','<','s','d'];
AlphaVec = [1 0.8 0.6 0.4 0.2];
MCxAxis = SnrVec_dB;
Snr4Str = SirVec_dB';
xAxisLabel = 'SNR[dB]';
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

