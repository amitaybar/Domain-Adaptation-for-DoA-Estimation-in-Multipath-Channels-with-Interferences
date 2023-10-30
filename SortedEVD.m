function [EigvecMat, EigValVec] = SortedEVD(Mat)

[Phi, LambdaMat]    = eig(Mat);
[~,ind]             = sort(diag(LambdaMat),'desc');
LambdaMat           = LambdaMat(ind,ind);
EigvecMat        	= Phi(:,ind);
EigValVec           = diag(LambdaMat);

end