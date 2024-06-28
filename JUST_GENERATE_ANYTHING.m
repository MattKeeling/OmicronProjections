function []=JUST_GENERATE_ANYTHING(FILE_STR, RUN_TYPE, CHANGE_DATE, Output, HOW_MANY_TO_RUN)
%%
%  JUST_GENERATE_ANYTHING(FILE_STR, RUN_TYPE, CHANGE_DATE, Output, HOW_MANY_TO_RUN)
%  
%
%  Output is [XL_output, Matlab_output] 0 or 1s
%%

OUTPUTS=[1:11];

if nargin==0
    RUN_TYPE=0;
end
RUN_TYPE=reshape(RUN_TYPE,1,[]);

if nargin<2
    CHANGE_DATE=floor(datenum(now));
end

W=whos('CHANGE_DATE');
if strncmp(W.class,'char',4) % input a date_str
    CHANGE_DATE=floor(datenum(CHANGE_DATE));
end

if CHANGE_DATE>datenum(2020,1,1)
    CHANGE_DATE=CHANGE_DATE+1-datenum(2020,1,1);
end


Remove_New_Variant=0;

Remove_BA2=1;

if nargin<3
    Output=[1 0];
end
XL_output=Output(1);
Matlab_Save=Output(2);

Rmax=8;
Rtot=[2:8];

if nargin<4
    HOW_MANY_TO_RUN=100;
end

RUN_UNTIL=datenum(2022,6,1);
%RUN_UNTIL=datenum(now)+10*7;

RESCALE=1;

MaxType=5;

%load r_Estimates.mat


nALL2PP=0; Transmission_Reduction=0;
VEffI=0; VEffS=0; VEffD=0; VEffH=0;
WaningSpeed=0; WaningEfficacy=0; UT=0; BUT=0; START_DATE=0;

KL=0;
for LOOP=1:length(RUN_TYPE)
    KL=KL+1;
    ZZ=1;
 
        load(['MCMC_' FILE_STR]);
        VERSION_TXT='8.1'; MODEL_NAME='FullMCMCODE_with_VoC';
        
%         nDefault_Vacc_Parameters; nAdd_Omicron3_VE;
%         [V1,V2,Phase, RatioPf]=Vaccination_over_Time_MTP_UK(UT);
%         [V3,Phase]=Boosters_over_Time_MTP_UK(BUT,V2);
        load Vaccination_Data.mat
        
        
        %%
        clear nALL nDEATHS SYMPTOMS HOSPITALISED nHOSP_AD nCASES PILLAR2
        
        load Starting_Data
        load Distributions2.mat
        Distribution_Symptoms_to_Hospital(1)=0;  % Don't beleive 1 day from Symp to Hosp.
        Distribution_Symptoms_to_Hospital=Distribution_Symptoms_to_Hospital/sum(Distribution_Symptoms_to_Hospital);
        load Probabilities_new.mat
        
        load Regional_PP.mat
        Region_PP(Region_PP==0)=100;
        
        Names={'England','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland','United Kingdom'};
        
        gamma=0.7; % just so it can be read in!
            
        
        %%
        % generate run_stops by region
        Comps = zeros(4,Num_Comp+1,11);
        
        tmpRS=RUN_STOP;  Z=(RUN_STOP(end-1)+7):7:(ceil(RUN_UNTIL+1-datenum(2020,1,1)));
        tmpRS(Num_Comp+[1:length(Z)]-1)=Z;
        
        RUN_STOPs = zeros(11,Num_Comp+1);
        for Region = 2:11
            RUN_STOPs(Region,1:length(tmpRS)) = tmpRS;
        end
        
        
        Keep_RUN_STOPs=RUN_STOPs(2,:);
        
        if size(V1,2)<max(RUN_STOPs,[],'all')
            mx=max(RUN_STOPs,[],'all')+100;
            for r=1:11  for a=1:21
                    V1(r,end:mx,a)=V1(r,end,a);
                    V2(r,end:mx,a)=V2(r,end,a);
                    RC(r,end:mx,a)=RC(r,end,a);
                    RCH(r,end:mx,a)=RCH(r,end,a);
                    RCD(r,end:mx,a)=RCD(r,end,a);
                end
            end
        end
         
        
        %%
        maxtime = max(RUN_STOPs(end,:))+30; % needs to be big enough to take account of any lags
        NQZ=105;
        nALL = zeros(11,maxtime,NQZ,length(nALPHA));
        nDEATHS = zeros(11,maxtime,NQZ,length(nALPHA));
        nHOSP_AD = zeros(11,maxtime,NQZ,length(nALPHA));
        nHOSP_OCC = zeros(11,maxtime,NQZ,length(nALPHA));
        nICU_AD = zeros(11,maxtime,NQZ,length(nALPHA));
        nICU_OCC = zeros(11,maxtime,NQZ,length(nALPHA));
        nCASES = zeros(11,maxtime,NQZ,length(nALPHA));
        PILLAR2 = zeros(11,maxtime,length(nALPHA));
        SWAB_PREV= zeros(11,maxtime,NQZ,length(nALPHA));
        REst = zeros(11,length(nALPHA));
        
        STEP=floor(length(nALPHA)/HOW_MANY_TO_RUN); STEP(STEP<1)=1;
        
        if length(nALPHA)<HOW_MANY_TO_RUN
            WHAT_TO_RUN=[length(nALPHA):-1:1];
        else
            WHAT_TO_RUN=sort(randperm(length(nALPHA),HOW_MANY_TO_RUN),'descend');
        end
       
        for QQ=WHAT_TO_RUN
            fprintf(1,'(%d %d) %d out of %d \n',ZZ,LOOP,QQ,length(nALPHA));
            
            H_STRETCH=nH_STRETCH(:,QQ);
            I_STRETCH=nI_STRETCH(:,QQ);
            
            FACTOR=nFACTOR(:,QQ);
            D_FACTOR=nD_FACTOR(:,QQ);
            H_FACTOR=nH_FACTOR(:,QQ);
            I_FACTOR=nI_FACTOR(:,QQ);
            TAU=nTAU(QQ);
            SCALING=nSCALING(:,QQ);
            COMPLIANCE=nCOMPLIANCE(:,:,QQ);
            try 
                COMPLIANCEO=nCOMPLIANCEO(:,:,QQ);
            catch me
                COMPLIANCEO=0*COMPLIANCE;
            end
            LOGLIKE=nALPHA(:,QQ);
            ALL2PP=nALL2PP(:,QQ);
            
            NV_SPEED=nNV_SPEED(:,QQ);
            NV_BETA=nNV_BETA(:,QQ);
            IMPORT_LEVEL=nIMPORT_LEVEL(:,QQ);
            IMPORT_DAY=nIMPORT_DAY(:,QQ);
            if Remove_BA2 IMPORT_DAY(34:end)=5000; end
            
            LAG=nLAG(:,QQ);
            INC_P=nINC_P(QQ);
            ALPHA=nALPHA(QQ);
            GAMMA=nGAMMA(QQ);
            
            DH_SCALING=nDH_SCALING(:,QQ);
            RE_INF=nRE_INF(:,QQ);
            number=nNumber(QQ);          
            
            Detection=nDetection(QQ,:);
            Susceptibility=nSusceptibility(QQ,:);
            gamma=nGAMMA(QQ);
            
            % generate Comp by region
            clear Comps
            
            RUN_STARTs=[82*ones(11,1) RUN_STOPs(:,1:(end-1))];
            
            Christmas=0;
            
            MAX_REGION=11;
                
            NC=Num_Comp;
            for Region = Rtot
                Comps(Region,1:NC)=COMPLIANCE(1:NC,Region);
                Comps(Region,NC:size(RUN_STOPs,2))=Comps(Region,NC);             
                CompsO(Region,1:NC)=COMPLIANCEO(1:NC,Region);
                CompsO(Region,NC:size(RUN_STOPs,2))=CompsO(Region,NC);                 
                MAX_REGION=Rmax;
            end
            
            if RUN_TYPE(LOOP)==0
                % Don't Change Comps, Output with 0
                TXT_STR=['Warwick_Omicron_MTP_' num2str(0)];
                EXTRA_STR='MTP';
            end
                       
            if RUN_TYPE(LOOP)>0
                % Run with fixed R values
                TXT_STR=['Warwick_Omicron_MTP' num2str(RUN_TYPE(LOOP),'%2.1f\n')]; TXT_STR(TXT_STR=='.')=[];
                EXTRA_STR=['MTP R' num2str(RUN_TYPE(LOOP),'%2.1f\n')];
                for Region= Rtot
                    Comps(Region,RUN_STOPs(Region,:)>=CHANGE_DATE+6)=-RUN_TYPE(LOOP);
                end
            end
            
            
            if RUN_TYPE(LOOP)==-1
                TXT_STR=['Warwick_Omicron_MTP_max'];
                EXTRA_STR='MTP';
                for Region= Rtot
                    %D=squeeze(mean(nCOMPLIANCE(NC+[-6:0],Region,:)-nCOMPLIANCE(NC+[-7:-1],Region,:),[1 3]));
                    %for p=(NC+1):size(Comps,2)
                    %    Comps(Region,p)=max(Comps(Region,p-1)+D,0.01);                     
                    %end
                    D=Comps(Region,RUN_STOPs(Region,:)>=CHANGE_DATE);
                    Comps(Region,RUN_STOPs(Region,:)>=CHANGE_DATE)=D(1);
                end
                Comps(Comps>0.99)=0.99;
            end
            
            
            Christmas=0;
            
            if Remove_New_Variant
                IMPORT_DAY(:)=max(max(RUN_STOPs))+100;
            end
            
            for Region=2:11
                NVEI(:,:,:,Region)=VEffI; NVES(:,:,:,Region)=VEffS;
                NVEH(:,:,:,Region)=VEffH; NVED(:,:,:,Region)=VEffD;         
                for j=1:MaxType
                    jj=j-1;
                    NVEI(:,j,4,Region)=1-(1-VEffI(:,j,3))*(1-RE_INF(Region+jj*11));
                    NVES(:,j,4,Region)=1-(1-VEffS(:,j,3))*(1-RE_INF(Region+jj*11));
                    NVEH(:,j,4,Region)=1-(1-VEffH(:,j,3))*(1-RE_INF(Region+jj*11));
                    NVED(:,j,4,Region)=1-(1-VEffD(:,j,3))*(1-RE_INF(Region+jj*11));
                end
            end
            
            %parfor Region=[2:8]
            for Region=2:8
                [nDC, nDCD, nDCH, nHospital, inHospital, nICU, inICU, nDeaths, ~, All, FS2]=SIMULATE_ONE_REGION(Region, TAU, ALPHA, INC_P, SCALING(Region), FACTOR(Region), H_FACTOR(Region:11:end), I_FACTOR(Region:11:end), D_FACTOR(Region:11:end),  ...
                H_STRETCH(Region:11:end), I_STRETCH(Region:11:end), LAG(Region:11:end), START_DATE(Region)+1, 1, Comps(Region,:), CompsO(Region,:), RUN_STOPs(Region,:), NV_BETA(Region:11:end)', NV_SPEED(Region:11:end)', IMPORT_DAY(Region:11:end),IMPORT_LEVEL(Region:11:end),...
                squeeze(V1(Region,:,:)), squeeze(V2(Region,:,:)), squeeze(V3(Region,:,:)), Transmission_Reduction, squeeze(NVEI(:,:,:,Region)), squeeze(NVES(:,:,:,Region)), squeeze(NVEH(:,:,:,Region)), squeeze(NVED(:,:,:,Region)), RatioPf,...
                Detection, Susceptibility, gamma, WaningSpeed, WaningEfficacy, [0.1 0], MaxType, [], 0);   
                
                nDeaths=nDeaths.* ((1+DH_SCALING(Region)*sum(inHospital,2)./sum(Region_PP(Region,:),2))*ones(1,size(nDeaths,2)));
                
                T=1:size(nDC,1);
                
                padding = maxtime-length(T);
                nALL(Region,:,:,QQ)=[abs(All); zeros(padding,NQZ)];
                nDEATHS(Region,:,:,QQ)=[abs(nDeaths); zeros(padding,NQZ)];
                nHOSP_OCC(Region,:,:,QQ)=[abs(inHospital); zeros(padding,NQZ)];
                nHOSP_AD(Region,:,:,QQ) = [abs(nHospital); zeros(padding,NQZ)];
                nCASES(Region,:,:,QQ) = [abs(nDC); zeros(padding,NQZ)];
                nICU_OCC(Region,:,:,QQ) = [abs(inICU); zeros(padding,NQZ)];
                nICU_AD(Region,:,:,QQ) = [abs(nICU); zeros(padding,NQZ)];
                              
                if RUN_TYPE(LOOP)>0
                    y=sum(abs(All),2);
                    r=(log(y(CHANGE_DATE+21))-log(y(CHANGE_DATE+7)))/14;
                    SS=NV_SPEED(Region:11:end);
                    R=max(((1+r./(3*INC_P*SS(end))).^3).*(1+r./(gamma*SS(end))),[],'omitnan');
                    REst(Region,QQ)=R;
                end
                
            end
            
            nALL(1,:,:,QQ)=sum(nALL(2:8,:,:,QQ),1);
            nALL(12,:,:,QQ)=sum(nALL(2:11,:,:,QQ),1);
            
        end
  
        %XL_output=1;
        
        %% GET RID OF ZEROS
        
        Q=sum(sum(sum(nALL,1),2),3);
        m=find(Q==0);
        nALL(:,:,:,m)=[];  nDEATHS(:,:,:,m)=[]; nHOSP_OCC(:,:,:,m)=[];
        nHOSP_AD(:,:,:,m)=[];  nCASES(:,:,:,m)=[]; 
        nICU_OCC(:,:,:,m)=[];  nICU_AD(:,:,:,m)=[];
        REst(:,m)=[];
        
        %% NEED TO COMBINE THE OLD & NEW VARIANT !
        KALL=nALL;
        
        A=1:size(nALL,3); a=1:21; b=21+a; c=42+a; d=63+a; e=84+a;
        oHOSP_AD=squeeze(sum(nHOSP_AD(:,:,d,:),3));
        nALL=nALL(:,:,a,:)+nALL(:,:,b,:)+nALL(:,:,c,:)+nALL(:,:,d,:)+nALL(:,:,e,:);
        nDEATHS=nDEATHS(:,:,a,:)+nDEATHS(:,:,b,:)+nDEATHS(:,:,c,:)+nDEATHS(:,:,d,:)+nDEATHS(:,:,e,:);
        nHOSP_OCC=nHOSP_OCC(:,:,a,:)+nHOSP_OCC(:,:,b,:)+nHOSP_OCC(:,:,c,:)+nHOSP_OCC(:,:,d,:)+nHOSP_OCC(:,:,e,:);
        nHOSP_AD=nHOSP_AD(:,:,a,:)+nHOSP_AD(:,:,b,:)+nHOSP_AD(:,:,c,:)+nHOSP_AD(:,:,d,:)+nHOSP_AD(:,:,e,:);
        nCASES=nCASES(:,:,a,:)+nCASES(:,:,b,:)+nCASES(:,:,c,:)+nCASES(:,:,d,:)+nCASES(:,:,e,:);
        nICU_OCC=nICU_OCC(:,:,a,:)+nICU_OCC(:,:,b,:)+nICU_OCC(:,:,c,:)+nICU_OCC(:,:,d,:)+nICU_OCC(:,:,e,:);
        nICU_AD=nICU_AD(:,:,a,:)+nICU_AD(:,:,b,:)+nICU_AD(:,:,c,:)+nICU_AD(:,:,d,:)+nICU_AD(:,:,e,:);  
        
        
        %% SAVE OUTPUTS TO SPREADSHEET
        
  
      
        
        if Matlab_Save
            
            %save(['For_Trystan_' DATE_STR '.mat'],'nALL','nDEATHS','nHOSP_OCC','nHOSP_AD','nCASES','SWAB_PREV','nICU_OCC','nICU_AD','RUN_STARTs');
           
            byAge_ALL=squeeze(mean(nALL(:,:,:,:),4));
            byAge_HospAd=squeeze(mean(nHOSP_AD(:,:,:,:),4));
            byAge_Death=squeeze(mean(nDEATHS(:,:,:,:),4));
 
            nALL=squeeze(sum(nALL(:,:,:,:),3));
            nDEATHS=squeeze(sum(nDEATHS(:,:,:,:),3));
            nHOSP_OCC=squeeze(sum(nHOSP_OCC(:,:,:,:),3));
            nHOSP_AD=squeeze(sum(nHOSP_AD(:,:,:,:),3));
            nCASES=squeeze(sum(nCASES(:,:,:,:),3));
            nICU_OCC=squeeze(sum(nICU_OCC(:,:,:,:),3));
            nICU_AD=squeeze(sum(nICU_AD(:,:,:,:),3));
            Ratio=squeeze(sum(KALL(2:8,:,21+[1:21],:),[1 3])+sum(KALL(2:8,:,63+[1:21],:),[1 3]))./squeeze(sum(KALL(2:8,:,:,:),[1 3]));
            Ratio2=squeeze(sum(KALL(:,:,21+[1:21],:),[3])+sum(KALL(:,:,63+[1:21],:),[3]))./squeeze(sum(KALL(:,:,:,:),[3]));
            
            %save([TXT_STR '_' DATE_STR '.mat'],'WHAT_TO_RUN','KALL','nALL','nDEATHS','nHOSP_OCC','oHOSP_AD','nHOSP_AD','nCASES','SWAB_PREV','nICU_OCC','nICU_AD','RUN_STARTs','RUN_STOP','byAge*');
            save(['Output_' FILE_STR],'WHAT_TO_RUN','Ratio','Ratio2','nALL','nDEATHS','nHOSP_OCC','oHOSP_AD','nHOSP_AD');
            %save([TXT_STR '_' DATE_STR '.mat'],'RUN_STARTs','RUN_STOP','byAge*');
            
            clear nALL nDEATHS nHOSP_OCC nHOSP_AD nCASES SWAB_PREV nICU_OCC nICU_AD
            
        end
        
        
    end
    
end