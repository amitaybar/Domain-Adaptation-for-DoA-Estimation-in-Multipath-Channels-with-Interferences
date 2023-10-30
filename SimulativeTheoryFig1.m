close all
clear
clc
%%
M           = 9;
KDim        = M;
Array       = 0:M-1;
SirLin1     = 1;
I           = 0.01*eye(M);
Theta_I     = 30;
Theta_I2    = 100;
%% 
for ii=0:1
    clearvars -except M KDim Array SirLin1 I Theta_I Theta_I2 ii
    rng(11,'twister');
    IsTwoInterferences = ii;
    %%
    a_I         = exp(2j*pi*Array'/2*cos(Theta_I/180*pi));
    a_I2        = IsTwoInterferences*exp(2j*pi*Array'/2*cos(Theta_I2/180*pi));
    K = 100;
    for k=1:K
        a_j = exp(2j*pi*Array'/2*cos(180*rand/180*pi));
        a_l = exp(2j*pi*Array'/2*cos(180*rand/180*pi));
        Gamma_j(:,:,k)     = I + a_j*a_j';
        Gamma_l(:,:,k)     = I + a_l*a_l' + SirLin1^2*a_I*a_I' + a_I2*a_I2';
    end
    Gamma_T = mean(Gamma_j,3);
    Gamma_S = mean(Gamma_l,3);
    
    ThetaVec = linspace(0,180,180*4+1);
    for tt=1:length(ThetaVec)
        a_s = exp(2j*pi*Array'/2*cos(ThetaVec(tt)/180*pi));      
        Gamma       = I + (SirLin1^2*a_I*a_I' + a_I2*a_I2' + a_s*a_s');       
        
        E_pt        = (Gamma_T*Gamma_S^-1)^0.5;
        E_Coral     = Gamma_T^0.5*Gamma_S^-0.5;
       
        if sum(abs(a_I2))>0            
            Res_pt(tt)      = real(a_s'*E_pt*Gamma*E_pt'*a_s) / real(a_I'*E_pt*Gamma*E_pt'*a_I) ;
            Res_coral(tt)   = real(a_s'*E_Coral*Gamma*E_Coral'*a_s) / real(a_I'*E_Coral*Gamma*E_Coral'*a_I);
        else
            Res_pt(tt)      = ( real(a_s'*E_pt*Gamma*E_pt'*a_s) / real(a_I'*E_pt*Gamma*E_pt'*a_I));
            Res_coral(tt)   = ( real(a_s'*E_Coral*Gamma*E_Coral'*a_s) / real(a_I'*E_Coral*Gamma*E_Coral'*a_I) );
        end
    end
    
    figure; hold on
    plot(ThetaVec,10*log10(Res_coral),'LineWidth',3);
    plot(ThetaVec,10*log10(Res_pt),'--','LineWidth',3);
    xline(Theta_I,'LineWidth',3);
    legend('Proposed','PT','Interference');
    if sum(abs(a_I2))>0
        xline(Theta_I2,'LineWidth',3);
        legend('Proposed','PT','Interference');
    end    
    xlabel('[deg]');
    ylabel('[dB]');    
    xlim([0 180]);
    ylim([0 12]);
    box on; grid on
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Times';
end


