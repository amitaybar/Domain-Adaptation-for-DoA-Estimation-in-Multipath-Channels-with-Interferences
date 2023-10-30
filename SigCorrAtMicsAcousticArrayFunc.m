
function [Sig,Corr] = SigCorrAtMicsAcousticArrayFunc(SourcePos,MicsPosMat,KMics,Beta,fs,BinInd,RoomDim,Order,KSamplesRir,WinLength,Signal,Snr)

SigLength = length(Signal);
hr  =  RIR(SourcePos,fs,MicsPosMat,KSamplesRir,Beta,RoomDim,Order);

SigAtMic  =  zeros(SigLength,KMics);
for m = 1:KMics
    SigAtMic(:,m) = filter(hr(m,:),1,Signal);
end

nl                      = randn(SigLength,KMics);
varn                    = mean(var(nl));
vars                    = mean(var(SigAtMic));
Gi_Tr                   = sqrt(vars/varn*10^(-Snr/10));
SigNoise                     = SigAtMic+ Gi_Tr*nl;
%%

for m=1:KMics
    SigTrPerMicStft            = stft(SigNoise(:,m),WinLength);
    SigTrPerMicStftFreq(m,:) = SigTrPerMicStft(BinInd,:);
end
Sig = conj(SigTrPerMicStftFreq)';
 
Corr                    = 1/size(Sig,2)*conj(Sig'*Sig);

end



