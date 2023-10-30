
function [Sig,Corr] = SigCorrAtMicsArrayFunc(SourcePos,MicsPosMat,KMics,Lambda,SigLength,SigTrain,SnrTrain)


%%
RTmp                    = SourcePos - MicsPosMat;
R                       = sqrt(diag(RTmp*RTmp'));



RandPhase               = 2*pi*rand(1,1);


hLosVec                 = 1./R.*exp(-2j*pi*1/(Lambda).*R);


SigAtMic                = hLosVec*SigTrain' ;
SigAtMic                = conj(SigAtMic');
nl                      = randn(SigLength,KMics);

varn                    = mean(var(nl));
vars                    = mean(var(SigAtMic));
Gi_Tr                   = sqrt(vars/varn*10^(-SnrTrain/10));

Sig                     = SigAtMic+ Gi_Tr*nl;

Corr                    = 1/size(Sig,2)*conj(Sig'*Sig);

end