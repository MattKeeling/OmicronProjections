clear

number=1;

Length1=2000;
Length2=1000;

Nw=now; NOW=Nw;
% Re_Fit_Earlier(number, datenum(2021,12,1), Length1, [0 1 0 0 0], 'Alpha');
% T=(now-Nw)*60*60*24; fprintf(1,'Alpha, %d steps in %02d:%02d:%02d\n',Length1,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Earlier(number, datenum(2021,12,1), Length1, [0 0 1 0 0], 'Delta');
!cp ReFit_MCMC_Latest1.mat ReFit_MCMC_Delta1.mat
T=(now-Nw)*60*60*24; fprintf(1,'Delta, %d steps in %02d:%02d:%02d\n',Length1,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2021,12,8), Length2, datenum(2021,12,7), 'NV_BETA(23:33,1)=NV_BETA(12:22)*2.8; NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 8/12/21, %d steps in %02d:%02d:%02d\n',Length1,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2021,12,13), Length2, datenum(2021,12,12), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 13/12/21, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2021,12,18), Length2, datenum(2021,12,17), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 18/12/21, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2021,12,23), Length2, datenum(2021,12,22), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 23/12/21, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;

Re_Fit_Omicron( number , datenum(2021,12,28), Length2, datenum(2021,12,27), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 28/12/21, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2022,1,2), Length2, datenum(2022,1,1), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 02/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2022,1,7), Length2, datenum(2022,01,6), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 07/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2022,1,12), Length2, datenum(2022,01,11), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 12/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron( number , datenum(2022,1,17), Length2, datenum(2022,01,16), 'NV_SPEED(23:33,1)=NV_SPEED(12:22); H_FACTOR([34:44],1)=H_FACTOR([23:33],1); I_FACTOR([34:44],1)=I_FACTOR([23:33],1); D_FACTOR([34:44],1)=D_FACTOR([23:33],1); H_STRETCH([34:44],1)=H_STRETCH([23:33],1); I_STRETCH([34:44],1)=I_STRETCH([23:33],1); ALL_2_PILLAR_POS([34:44],1)=ALL_2_PILLAR_POS([23:33],1);', 'NV_Speed3=NV_Speed2; Re_Inf(4)=1.0;', 'h_factor(4)=h_factor(3); i_factor(4)=i_factor(3); d_factor(4)=d_factor(3); h_stretch(4)=h_stretch(3); i_stretch(4)=i_stretch(3);' , 'Delta_P');
T=(now-Nw)*60*60*24; fprintf(1,'Fixed 17/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;

!cp ReFit_MCMC_Latest21_12_28_Delta_P.mat ReFit2_MCMC_Latest1.mat

Re_Fit_Omicron2( number , datenum(2021,12,23), Length2, datenum(2021,12,22), 'IMPORT_DAY(34:end,1)=5000;', '', '' , 'Free_P' );
T=(now-Nw)*60*60*24; fprintf(1,'Variable 23/12/21, %d steps in %02d:%02d:%02d\n',Length1,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron2( number , datenum(2021,12,28), Length2, datenum(2021,12,27), 'IMPORT_DAY(34:end,1)=5000;', '', '' , 'Free_P' );
T=(now-Nw)*60*60*24; fprintf(1,'Variable 28/12/21, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron2( number , datenum(2022,1,2), Length2, datenum(2022,1,1), 'IMPORT_DAY(34:end,1)=5000;', '', '' , 'Free_P' );
T=(now-Nw)*60*60*24; fprintf(1,'Variable 02/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron2( number , datenum(2022,1,7), Length2, datenum(2022,1,6), 'IMPORT_DAY(34:end,1)=5000;', '', '' , 'Free_P' );
T=(now-Nw)*60*60*24; fprintf(1,'Variable 07/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron2( number , datenum(2022,1,12), Length2, datenum(2022,1,11), 'IMPORT_DAY(34:end,1)=5000;', '', '' , 'Free_P' );
T=(now-Nw)*60*60*24; fprintf(1,'Variable 12/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;
Re_Fit_Omicron2( number , datenum(2022,1,17), Length2, datenum(2022,1,16), 'IMPORT_DAY(34:end,1)=5000;', '', '' , 'Free_P' );
T=(now-Nw)*60*60*24; fprintf(1,'Variable 17/01/22, %d steps in %02d:%02d:%02d\n',Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60))); Nw=now;

T=(now-NOW)*60*60*24; fprintf(1,'ALL runs, using %d and %d steps in %02d:%02d:%02d\n',Length1,Length2,floor(T/3600),floor((T-floor(T/3600)*3600)/60),round(mod(T,60)));

%%
clear Names;
D=dir('SALL4*_Delta_P_*.mat'); for i=1:length(D), Names{i}=D(i).name; end; k=i;
D=dir('SALL4*_Free_P_*.mat'); for i=1:length(D), Names{i+k}=D(i).name; end;

%Names{1}='SALL4_MCMC_01_12_21_Alpha_1.mat';
%Names{1}='SALL4_MCMC_01_12_21_Delta_1.mat';

Sort_Data_to_Match_28_10_20;

N=1;
for i=1:length(Names)
    load(Names{i},'HOSP_AD');
    subplot(5,4,i);
    plot([1:size(True_AdHosp,2)],sum(True_AdHosp(2:8,:)),'.b','MarkerSize',15);
    hold on
    plot([1:size(HOSP_AD,2)],squeeze(sum(HOSP_AD(2:8,:,:))),'-r');
    hold off
    axis([650 820 0 5000]);
    txt=Names{i}; txt(txt=='_')=' ';
    title(txt);
end

%%

%%
clear Names;
%D=dir('SALL4*_Delta_P_*.mat'); for i=1:length(D), Names{i}=D(i).name; end; k=i;
D=dir('SALL4*_Free_P_*.mat'); for i=1:length(D), Names{i}=D(i).name; end;

N=1;
for i=1:length(Names)
    load(Names{i});
    n=find(LOGLIKE(2,:)<0,1,'last'); m=1:n; m=n;
    M=1:length(m);  OLD_VAR(11,end,end)=0; NEW_VAR(11,end,end)=0;
    
    nCASES=CASES(:,:,m);
    nACCEPT=ACCEPTED(:,m);
    nACCEPT_N=ACCEPTED_N(:,m);
    nALL=ALL(:,:,m);
    nOLD_VAR=OLD_VAR(:,:,m);
    nNEW_VAR=NEW_VAR(:,:,m);
    nALPHA=ALPHA(m);
    nCASES=CASES(:,:,m);
    nCOMPLIANCE=COMPLIANCE(:,:,m);
    nCOMPLIANCEO=COMPLIANCEO(:,:,m);
    nALL2PP=ALL_2_PILLAR_POS(:,m);
    nDEATHS=DEATHS(:,:,m);
    nFACTOR=FACTOR(:,m);
    nH_FACTOR=H_FACTOR(:,m);
    nI_FACTOR=I_FACTOR(:,m);
    nD_FACTOR=D_FACTOR(:,m);
    nH_STRETCH=H_STRETCH(:,m);
    nI_STRETCH=I_STRETCH(:,m);
    nHOSP_TOTAL=HOSP_TOTAL(:,:,m);
    nICU_TOTAL=ICU_TOTAL(:,:,m);
    nHOSP_AD=HOSP_AD(:,:,m);
    nICU_AD=ICU_AD(:,:,m);
    nLOGLIKE=LOGLIKE(:,m);
    nSCALING=SCALING(:,m);
    %nPOS=POS(:,:,m);
    nSENS=SENS(m);
    %nTPOS=TPOS(:,:,m);
    nTAU=TAU(m);
    nINC_P=INC_P(m);
    nLAG=LAG(:,m);
    nGAMMA=GAMMA(m);
    nNV_BETA=NV_BETA(:,m);
    nNV_SPEED=NV_SPEED(:,m);
    nIMPORT_DAY=IMPORT_DAY(:,m);
    nIMPORT_LEVEL=IMPORT_LEVEL(:,m);
    nSDEL_BASE=SDEL_BASE(:,m);
    nDH_SCALING=DH_SCALING(:,m); DH_flag=1;
    nRE_INF=RE_INF(:,m);
    
    nSERO_DK=SERO_DECAY(m);
    nNumber=N+0*m;
    niter=m;
    nNum_Comp=Num_Comp;
    
    txt=Names{i}; txt(1:11)=[]; txt(10:end)=[];
    %txt=[txt 'Delta_Like.mat'];
    txt=[txt 'Free_Params.mat'];
    eval(['save MCMC_' txt ' n* START_DATE Num_Comp RUN_STOP']);
end

%%
clear Names;
D=dir('SALL4*_Delta_P_*.mat'); for i=1:length(D), Names{i}=D(i).name; end; k=i;
D=dir('SALL4*_Free_P_*.mat'); for i=1:length(D), Names{i+k}=D(i).name; end;

%Names{1}='SALL4_MCMC_01_12_21_Alpha_1.mat';
%Names{1}='SALL4_MCMC_01_12_21_Delta_1.mat';

Sort_Data_to_Match_28_10_20;

N=1;
for i=1:length(Names)
    load(Names{i});
    n=find(LOGLIKE(2,:)<0,1,'last'); m=1:n; m=ceil(n*0.5):n; m=[10:20:200];
    M=1:length(m);  OLD_VAR(11,end,end)=0; NEW_VAR(11,end,end)=0;
    
    nCASES=CASES(:,:,m);
    nACCEPT=ACCEPTED(:,m);
    nACCEPT_N=ACCEPTED_N(:,m);
    nALL=ALL(:,:,m);
    nOLD_VAR=OLD_VAR(:,:,m);
    nNEW_VAR=NEW_VAR(:,:,m);
    nALPHA=ALPHA(m);
    nCASES=CASES(:,:,m);
    nCOMPLIANCE=COMPLIANCE(:,:,m);
    nCOMPLIANCEO=COMPLIANCEO(:,:,m);
    nALL2PP=ALL_2_PILLAR_POS(:,m);
    nDEATHS=DEATHS(:,:,m);
    nFACTOR=FACTOR(:,m);
    nH_FACTOR=H_FACTOR(:,m);
    nI_FACTOR=I_FACTOR(:,m);
    nD_FACTOR=D_FACTOR(:,m);
    nH_STRETCH=H_STRETCH(:,m);
    nI_STRETCH=I_STRETCH(:,m);
    nHOSP_TOTAL=HOSP_TOTAL(:,:,m);
    nICU_TOTAL=ICU_TOTAL(:,:,m);
    nHOSP_AD=HOSP_AD(:,:,m);
    nICU_AD=ICU_AD(:,:,m);
    nLOGLIKE=LOGLIKE(:,m);
    nSCALING=SCALING(:,m);
    %nPOS=POS(:,:,m);
    nSENS=SENS(m);
    %nTPOS=TPOS(:,:,m);
    nTAU=TAU(m);
    nINC_P=INC_P(m);
    nLAG=LAG(:,m);
    nGAMMA=GAMMA(m);
    nNV_BETA=NV_BETA(:,m);
    nNV_SPEED=NV_SPEED(:,m);
    nIMPORT_DAY=IMPORT_DAY(:,m);
    nIMPORT_LEVEL=IMPORT_LEVEL(:,m);
    nSDEL_BASE=SDEL_BASE(:,m);
    nDH_SCALING=DH_SCALING(:,m); DH_flag=1;
    nRE_INF=RE_INF(:,m);
    
    nSERO_DK=SERO_DECAY(m);
    nNumber=N+0*m;
    niter=m;
    nNum_Comp=Num_Comp;
    
    txt=Names{i}; txt(1:11)=[];
    eval(['save SERO4_' txt ' n* START_DATE Num_Comp RUN_STOP']);
    
    %JUST_GENERATE_ANYTHING(txt, 0, [], [0 1], 10);
 
    load(['Output_' txt],'nHOSP_AD','oHOSP_AD','Ratio');
    
    subplot(5,3,i);
    plot([1:size(True_AdHosp,2)],sum(True_AdHosp(2:8,:)),'.','MarkerSize',15);
    hold on
    plot([1:size(nHOSP_AD,2)],squeeze(sum(nHOSP_AD(2:8,:,:))),'-');
    hold off
    set(gca,'XLim',[650 820]);
    
    clear n*
end
   

%%
clf;
D=datenum(2021,12,8)+[0:5:40]; d=D+1-datenum(2020,1,1);

for i=1:length(D)
    try
        load(['Output_' datestr(D(i),'yy_mm_dd') '_Delta_P_1'],'nHOSP_AD','oHOSP_AD','Ratio');
        subplot(4,9,i);
        plot([1:d(i)],sum(True_AdHosp(2:8,1:d(i))),'.c','MarkerSize',10);
        hold on
        plot([(d(i)+1):size(True_AdHosp,2)],sum(True_AdHosp(2:8,(d(i)+1):end)),'ok','MarkerSize',3);
        plot([1:size(nHOSP_AD,2)],squeeze(sum(nHOSP_AD(2:8,:,:))),'-r');
        hold off
        title(datestr(D(i)-1,'dd mmm yy'));
        axis([650 800 0 4000]);
        
        subplot(4,9,i+18)
        plot([1:size(Ratio,1)],squeeze(mean(Ratio(:,2:8),2)),'-r');
        axis([650 800 0 1]);
    end
    
    try
        load(['Output_' datestr(D(i),'yy_mm_dd') '_Free_P_1'],'nHOSP_AD','oHOSP_AD','Ratio');
        subplot(4,9,i+9);
        plot([1:d(i)],sum(True_AdHosp(2:8,1:d(i))),'c.','MarkerSize',10);
        hold on
        plot([(d(i)+1):size(True_AdHosp,2)],sum(True_AdHosp(2:8,(d(i)+1):end)),'ok','MarkerSize',3);
        plot([1:size(nHOSP_AD,2)],squeeze(sum(nHOSP_AD(2:8,:,:))),'-r');
        hold off
        axis([650 800 0 4000]);
        
        subplot(4,9,i+27)
        plot([1:size(Ratio,1)],squeeze(mean(Ratio(:,2:8),2)),'-r');
        axis([650 800 0 1]);
    end
end
    
    
    
    
