% VI=[0.63 0.80];
% VS=[0.63 0.87];
% VS_fiddle=[0.86 0.91];
% VH=[0.8 0.92];v
% VD=[0.8 0.96];

% Vaccine Efficacy Against ALPHA (dose 1,2 Pf, dose 1,2 AZ)
VIPfAZ=[63 80 63 78]/100;
VSPfAZ=[63 88 63 80]/100;
VHPfAZ=[80 93 80 90]/100;
VDPfAZ=[80 97 80 95]/100;

% Vaccine Efficacy Against DELTA  (dose 1,2 Pf, dose 1,2 AZ)
Dlta_VIPfAZ=[55 85 45 70]/100;
Dlta_VSPfAZ=[55 90 45 70]/100;
Dlta_VHPfAZ=[80 95 80 95]/100;
Dlta_VDPfAZ=[80 98 80 98]/100;


%Wildtype  VEff(Pf / AZ, Varient, dose)
VEffI(1,1,1:2)=VIPfAZ([1 2]); VEffI(2,1,1:2)=VIPfAZ([3 4]);
VEffS(1,1,1:2)=VSPfAZ([1 2]); VEffS(2,1,1:2)=VSPfAZ([3 4]);
VEffH(1,1,1:2)=VHPfAZ([1 2]); VEffH(2,1,1:2)=VHPfAZ([3 4]);
VEffD(1,1,1:2)=VDPfAZ([1 2]); VEffD(2,1,1:2)=VDPfAZ([3 4]);
%Alpha
VEffI(1,2,1:2)=VIPfAZ([1 2]); VEffI(2,2,1:2)=VIPfAZ([3 4]);
VEffS(1,2,1:2)=VSPfAZ([1 2]); VEffS(2,2,1:2)=VSPfAZ([3 4]);
VEffH(1,2,1:2)=VHPfAZ([1 2]); VEffH(2,2,1:2)=VHPfAZ([3 4]);
VEffD(1,2,1:2)=VDPfAZ([1 2]); VEffD(2,2,1:2)=VDPfAZ([3 4]);
%Delta
VEffI(1,3,1:2)=Dlta_VIPfAZ([1 2]); VEffI(2,3,1:2)=Dlta_VIPfAZ([3 4]);
VEffS(1,3,1:2)=Dlta_VSPfAZ([1 2]); VEffS(2,3,1:2)=Dlta_VSPfAZ([3 4]);
VEffH(1,3,1:2)=Dlta_VHPfAZ([1 2]); VEffH(2,3,1:2)=Dlta_VHPfAZ([3 4]);
VEffD(1,3,1:2)=Dlta_VDPfAZ([1 2]); VEffD(2,3,1:2)=Dlta_VDPfAZ([3 4]);

%Booster = Better than 2nd dose, by factor 1-F
F=0.8;
VEffI(1,:,3)=1-F*(1-VEffI(1,:,2)); VEffI(2,:,3)=1-F*(1-VEffI(1,:,2));
VEffS(1,:,3)=1-F*(1-VEffS(1,:,2)); VEffS(2,:,3)=1-F*(1-VEffS(1,:,2));
VEffH(1,:,3)=1-F*(1-VEffH(1,:,2)); VEffH(2,:,3)=1-F*(1-VEffH(1,:,2));
VEffD(1,:,3)=1-F*(1-VEffD(1,:,2)); VEffD(2,:,3)=1-F*(1-VEffD(1,:,2));

% And what happens for recovereds.
VEffI(:,:,4)=1-1e-4;
VEffS(:,:,4)=1-5e-5;
VEffH(:,:,4)=1;
VEffD(:,:,4)=1;

Transmission_VE(1:2,1:2,1:3)=0.45; Transmission_VE(1:2,3,1:3)=0.3; Transmission_VE(:,:,4)=0.9;
Transmission_Reduction=1-Transmission_VE;   % this is the transmission scaling higher is more transmission.

ALREADY_VACC=[6:21]; %25+
UT=0.95*ones(1,21); UT(15:21)=0.95; UT(5:10)=0.75; UT(4)=UT(5);  UT(3)=0.6*UT(5); % now with 12-15s.
UT(1:2)=0.0; 

BUT(15:21)=0.95;   % Boosters are 95% of uptake, in over 70
BUT(11:14)=0.90;   % Boosters are 90% of uptake, in over 50-69
BUT(1:10)=0.8;  % 80% in the under 50s.
BUT(1:4)=0.4*0.8;   % no boosting in under 18s.

UT(ALREADY_VACC)=0;


% ASSUME WANING TO ZERO - AS THIS GIVES A BETTER FIT.

% Waning is [VE_I(R) VE_I(V) OR_H(R) OR_H(V) OR_D(R) OR_D(V) Red_trans(R) Red_trans(V)]
% if ~exist('number')
%     WaningEfficacy=[0.4 0.4 0.3 0.3 0.3 0.3 0.8 0.8];
%     WaningSpeed=[0*0.005 1/20 0*0.005 1/250];
% else
%     if number<=10
%         WaningEfficacy=[0.5 0.5 0.3 0.3 0.3 0.3 0.8 0.8];
%         WaningSpeed=[0*0.005 1/60 0*0.005 1/120];
%     end
%     if number>10 & number<=20
%         WaningEfficacy=[0.3 0.3 0.3 0.3 0.3 0.3 0.8 0.8];
%         WaningSpeed=[0*0.005 1/60 0*0.005 1/250];
%     end
%     if number>20
        WaningEfficacy=[0.0 0.0 0.3 0.3 0.3 0.3 0.8 0.8];
        WaningSpeed=[0*0.005 1/60 0*0.005 1/400];
%     end
    WaningSpeed([1 3])=[1/60 1/250]/6; %% natural infection wanes over 5 years.
% end

WaningSpeed=WaningSpeed*1.25;  % Gives a better fit to the data.

for TYPE=2:3,  WaningEfficacy(TYPE,:)=WaningEfficacy(1,:); end


