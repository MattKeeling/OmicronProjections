function [Detection, Susceptibility, gamma] = FIT_ALPHA_TAU_MEX( alpha, tau, gamma, IncP )

% alpha=1  All heterogeneity on Susceptibility, Detection=90%
% alpha=0  All heterogeneity on Detection

if nargin<3
    gamma=0.5;
end
if nargin<4
    IncP=0.2;
end


%load Starting_Data
load To_Fit_Alpha_n_Tau

%Col=jet(100); clf;
%plot([1:21],Q,'*'); hold on;
for i=1:100
    d=Q.^(1-alpha);  s=Q.^alpha;
    
    Symp=(UK_Ages);  Asymp=(UK_Ages./d)-Symp;
    
    NI = Q.* ( (Symp + tau*Asymp)*UK_from_to );
    nQ = UK_Ages ./ ( (Symp + tau*Asymp)*UK_from_to );
    
    nQ=nQ/max(nQ);
    
    Q=0.9*Q + 0.1*nQ;
    %plot([1:21],Q,'.-b','Color',Col(i,:));
end
%plot([1:21],Q,'o'); hold off

Detection=0.9 * Q.^(1-alpha); Susceptibility=(1/0.9) * Q.^(alpha);

a=IncP;  % 1/latent period
dgamma=gamma*0.5; flag=0;
%gamma=1/2.35; % to get DT=3.3
for i=1:12
    Tmp_Matrix1=UK_from_to*gamma;
    %    [T,S,E,D,nD,U , R01, Da]=Extended_ODEs(Tmp_Matrix1, a, gamma, Susceptibility, Detection, tau, UK_PP, 1);
    [T,~,E,~,~,~ , R01, ~, ~] = newerExtended_ODEs_mex(UK_from_toH*0 , gamma*(UK_from_toH + UK_from_toW + UK_from_toS + UK_from_toO), a, gamma, Susceptibility, Detection, tau, 0, UK_PP , 2, -ones(1,441));
    
    Susceptibility=Susceptibility*2.7/R01; %fine tune if necessary
    %    [T,S,E,D,nD, U, R02, Da]=Extended_ODEs(Tmp_Matrix1, a, gamma, Susceptibility, Detection, tau, UK_PP, 20);
    [T,~,E,~,~,~ , R02, ~, ~] = newerExtended_ODEs_mex(UK_from_toH*0 , gamma*(UK_from_toH + UK_from_toW + UK_from_toS + UK_from_toO), a, gamma, Susceptibility, Detection, tau, 0, UK_PP , 20, -ones(1,441));
    
    g=mean(E(end,:)./E(end-1,:)); DT=log(2)/log(g);
       
    if DT>3.3
        gamma=gamma+dgamma;
    else
        gamma=gamma-dgamma; flag=1;
    end
    if flag
        dgamma=dgamma/2;
    else
        dgamma=dgamma/2;
    end
end


