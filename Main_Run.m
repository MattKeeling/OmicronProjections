%%
%   RUN THE FULL OMICRON CODE TO GENERAL SIMILAR RESULTS TO FIGURE 5
%%
clear

SAMPLES=1; %up to a maximum of 10.

load OBSERVED_HOSPITALISATIONS.mat

D=dir('MCMC*Delta.mat'); 
for i=1:length(D), Names{i}=D(i).name; end; k=i;

L=length(Names);

figure(1); set(gcf,'position',[0 730 1700 600]); clf;
figure(2); set(gcf,'position',[100 630 1700 600]); clf;
figure(3); set(gcf,'position',[200 530 1700 600]); clf;
drawnow;

X=datenum(2021,1:36,1); XT=X+1-datenum(2020,1,1); XTL=datestr(X,'mmm''yy'); XTL(2:2:18,:)=0;

for i=1:length(Names)
    txt=Names{i}; txt(1:5)=[];
    JUST_GENERATE_ANYTHING(txt, 0, [], [0 1], SAMPLES);

    load(['Output_' txt],'nHOSP_AD','Ratio2','nALL');
    
    txt(9:end)=[]; Date=datenum(txt,'yy_mm_dd');  Date2=Date+1-datenum(2020,1,1);
    figure(1); subplot(2,L,i)
    plot([1:Date2],OBSERVED_HOSPITALISATIONS(1:Date2)/1e3,'ok','MarkerFaceColor','k','MarkerSize',3);
    hold on
    plot([(Date2+1):length(OBSERVED_HOSPITALISATIONS)],OBSERVED_HOSPITALISATIONS((Date2+1):end)/1e3,'ok','MarkerFaceColor','y','MarkerSize',3);
    plot([1:size(nHOSP_AD,2)],squeeze(sum(nHOSP_AD(2:8,:,:),1))/1e3,'-b','LineWidth',2);
    hold off
    title(datestr(Date)); set(gca,'XLim',[650 820],'XTick',XT,'XTickLabel',XTL);
    if i==1 ylabel('Hospital Admissions (thousands)'); end
    
    figure(2); subplot(2,L,i);
    plot([1:size(Ratio2,2)],squeeze(mean(Ratio2(2:8,:,:),1))*100,'-b','LineWidth',2); 
    title(datestr(Date)); set(gca,'XLim',[650 820],'XTick',XT,'XTickLabel',XTL);
    if i==1 ylabel('Proportion Omicron (%)'); end
    
    figure(3); subplot(2,L,i);
    plot([1:size(nALL,2)],squeeze(sum(nALL(2:8,:,:),1))/1e6,'-b','LineWidth',2); 
    title(datestr(Date)); set(gca,'XLim',[650 820],'XTick',XT,'XTickLabel',XTL);
    if i==1 ylabel('Infections (millions)'); end
    drawnow;
end

clear Names;
D=dir('MCMC*Free.mat'); for i=1:length(D), Names{i}=D(i).name; end;
l=length(Names);

for i=1:length(Names)
    txt=Names{i}; txt(1:5)=[];
    JUST_GENERATE_ANYTHING(txt, 0, [], [0 1], SAMPLES);

    load(['Output_' txt],'nHOSP_AD','Ratio2','nALL');
    
    txt(9:end)=[]; Date=datenum(txt,'yy_mm_dd');  Date2=Date+1-datenum(2020,1,1);
    figure(1); subplot(2,L,i+2*L-l);
    plot([1:Date2],OBSERVED_HOSPITALISATIONS(1:Date2)/1e3,'ok','MarkerFaceColor','k','MarkerSize',3);
    hold on
    plot([(Date2+1):length(OBSERVED_HOSPITALISATIONS)],OBSERVED_HOSPITALISATIONS((Date2+1):end)/1e3,'ok','MarkerFaceColor','y','MarkerSize',3);
    plot([1:size(nHOSP_AD,2)],squeeze(sum(nHOSP_AD(2:8,:,:),1))/1e3,'-r','LineWidth',2);
    hold off
    title(datestr(Date));
    set(gca,'XLim',[650 820],'XTick',XT,'XTickLabel',XTL);
    if i==1 ylabel('Hospital Admissions (thousands)'); end
    
    figure(2); subplot(2,L,i+2*L-l);
    plot([1:size(Ratio2,2)],squeeze(mean(Ratio2(2:8,:,:),1))*100,'-r','LineWidth',2); 
    title(datestr(Date)); set(gca,'XLim',[650 820],'XTick',XT,'XTickLabel',XTL);
    if i==1 ylabel('Proportion Omicron (%)'); end
    
    figure(3); subplot(2,L,i+2*L-l);
    plot([1:size(nALL,2)],squeeze(sum(nALL(2:8,:,:),1))/1e6,'-r','LineWidth',2); 
    title(datestr(Date)); set(gca,'XLim',[650 820],'XTick',XT,'XTickLabel',XTL);
    if i==1 ylabel('Infections (millions)'); end
    drawnow;
end


