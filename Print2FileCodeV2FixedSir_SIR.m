%% Functions
RoundFunc          = @(x) ( round(x*10)/10 ) ;
%%

% FileNameLabel = 'AcoustNewTable';
%%
FileId = fopen(fullfile('results',['Results_V_' FileNameLabel '.txt']),'w');
%-------------------SNR 20dB----------------------------
fprintf(FileId,['V1' num2str(1)]);
fprintf(FileId,'\n');
fprintf(FileId,'\n');

fprintf(FileId,' \\hline ');fprintf(FileId,'\n');

for sirInd = size(ThetaEstTs_ptErr_SNR,3):-1:1
    % DS
    fprintf(FileId,'  SIR=$%d\\text{dB}$  ',SirVec_dB(sirInd));
    fprintf(FileId,' & ');
    for snrInd=1:size(ThetaEstTs_ptErr_SNR,2)
        MedianData      = median(ThetaEstTsErr_SNR(:,snrInd,sirInd));
        MedianData_pt   = median(ThetaEstTs_ptErr_SNR(:,snrInd,sirInd));
        IqrData         = iqr(ThetaEstTsErr_SNR(:,snrInd,sirInd));
        IqrData_pt      = iqr(ThetaEstTs_ptErr_SNR(:,snrInd,sirInd));
        fprintf(FileId,'%.1f (%.1f)', RoundFunc(MedianData),RoundFunc(IqrData));
        fprintf(FileId,'     & ');
        fprintf(FileId,'\\textbf{%.1f} (%.1f)', RoundFunc(MedianData_pt),RoundFunc(IqrData_pt));
        fprintf(FileId,'     & ');
    end
    % MVDR
    for snrInd=1:size(ThetaEstTs_ptErr_SNR,2)
        MedianData      = median(ThetaEstMvdrTsErr_SNR(:,snrInd,sirInd));
        MedianData_pt   = median(ThetaEstMvdrTs_ptErr_SNR(:,snrInd,sirInd));
        IqrData         = iqr(ThetaEstMvdrTsErr_SNR(:,snrInd,sirInd));
        IqrData_pt      = iqr(ThetaEstMvdrTs_ptErr_SNR(:,snrInd,sirInd));
        fprintf(FileId,'%.1f (%.1f)', RoundFunc(MedianData),RoundFunc(IqrData));
        fprintf(FileId,'     & ');
        fprintf(FileId,'\\textbf{%.1f} (%.1f)', RoundFunc(MedianData_pt),RoundFunc(IqrData_pt));
        fprintf(FileId,'     & ');
    end
    % MUSIC
    for snrInd=1:size(ThetaEstTs_ptErr_SNR,2)
        MedianData      = median(ThetaEstMusicTsErr_SNR(:,snrInd,sirInd));
        MedianData_pt   = median(ThetaEstMusicTs_ptErr_SNR(:,snrInd,sirInd));
        IqrData         = iqr(ThetaEstMusicTsErr_SNR(:,snrInd,sirInd));
        IqrData_pt      = iqr(ThetaEstMusicTs_ptErr_SNR(:,snrInd,sirInd));
        fprintf(FileId,'%.1f (%.1f)', RoundFunc(MedianData),RoundFunc(IqrData));
        fprintf(FileId,'     & ');
        fprintf(FileId,'\\textbf{%.1f} (%.1f)', RoundFunc(MedianData_pt),RoundFunc(IqrData_pt));
        %         fprintf(FileId,'     & ');
        if snrInd > 1
            %             fprintf(FileId,'     & ');
        end
    end
    fprintf(FileId,'  \\\\');
    
    fprintf(FileId,'\n');
    fprintf(FileId,'\n');
    
    if snrInd < size(ThetaEstTs_ptErr_SNR,2)
        fprintf(FileId,' \\hline ');fprintf(FileId,'\n');
    end
end
 
 
%%
fprintf(FileId,'\n');
fprintf(FileId,'\n');
fprintf(FileId,'\n');

fprintf(FileId,'KInterferencePos = %.2f  ;  KTrainPoints = %.2f  ;  KTestPoints = %.2f ;  KTestPoints = %.2f \n',KInterferencePos,KTrainPoints,KTestPoints,KArtPoints);
fprintf(FileId,'IsInterferenceFixed = %.2f  ;  IsTwoInterferenceSourcesActive = %.2f \n',IsInterferenceFixed, IsTwoInterferenceSourcesActive);
if exist('Beta','var')
    fprintf(FileId,'Beta = %.2f',Beta);
end


%%
fclose(FileId);

