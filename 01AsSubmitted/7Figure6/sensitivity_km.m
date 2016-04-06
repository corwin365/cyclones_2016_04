clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
AllMLS    = load([LocalDataDir,'/corwin/storms/results/storms_MLS_atz.mat']);    
AllSABER  = load([LocalDataDir,'/corwin/storms/results/storms_SABER_atz.mat']);  
AllHIRDLS = load([LocalDataDir,'/corwin/storms/results/storms_HIRDLS_atz.mat']); 

SABER = AllSABER.SABER.Results;
HIRDLS = AllHIRDLS.HIRDLS.Results;
MLS = AllMLS.MLS.Results;


disp('Data loaded')

%remove anything more than N days from the peak

for iInst=1:1:3;
  
  switch iInst;
    case 1; Data = HIRDLS;
    case 2; Data = SABER;
    case 3; Data = MLS;
  end
  
  Good = find(Data.StormIDs(5,:) >= -1 &  Data.StormIDs(5,:) <= 2);
  Data.StormIDs = Data.StormIDs(:,Good);
  Data.Results  = Data.Results( :,Good);
  
  Data.Results(2,:) = abs(Data.Results(2,:)).*1000; %Pa -> mPa
 

  switch iInst;
    case 1; HIRDLS = Data;
    case 2; SABER = Data;
    case 3; MLS = Data;
  end  
  
end

disp('Data trimmed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bootstrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%loop over distances and find mean T' at this range
Distances = 0:10:1200; %km
Percentiles  = [ 2.5,18,50,82,97.5];
Results   = NaN(3,numel(Distances),numel(Percentiles));
NPoints   = NaN(3,numel(Distances));


NStraps = 1000;

for iInst = 1:1:3;
  iInst
  switch iInst;
    case 1; Data = HIRDLS;
    case 2; Data = SABER;
    case 3; Data = MLS;
  end  
  
  for iDist = 1:1:numel(Distances);
    Good = find(Data.StormIDs(6,:)./1000 >= Distances(iDist)-50 ...
              & Data.StormIDs(6,:)./1000 <= Distances(iDist)+50 ...
              & ~isnan(Data.Results(2,:)));
            
%             stop
    NPoints(iInst,iDist) = numel(Good);         
    if numel(Good) == 0; continue; end
    
    %bootstrap results
    Straps = NaN(NStraps,numel(Good));

    
    for iStrap=1:1:NStraps;
      TheseStraps = randi(numel(Good),numel(Good),1);
      Straps(iStrap,:) = Data.Results(2,Good(TheseStraps));            
    end
    Straps = nanmedian(Straps,2);
    Results(iInst,iDist,:) = prctile(Straps,Percentiles);
  end

end


RawResults = Results;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
set(gcf,'color','w')

Colours = 'bkr';

Results = RawResults;
hold on

Basis = 300;

for i=0.7:0.1:1.1; plot([0 1500],[i i],'--','linewidth',1,'color',[1,1,1].*0.7); end
for i=100:100:1100;  plot([i i],[0 10],'--','linewidth',1,'color',[1,1,1].*0.7); end

for iInst=1:1:3;


  %scale to basis value
  [~,idx] = min(abs(Distances-Basis));
  Results(iInst,:,:) = Results(iInst,:,:)./Results(iInst,idx,3);

  %5th - 95th percentile
  X = [Distances,reverse(Distances),Distances(1)];
  Y = [Results(iInst,:,1),reverse(Results(iInst,:,5)),Results(iInst,1,1)];
  fill(X,Y,Colours(iInst),'edgecolor','none','facealpha',0.3);

  %32nd-68th percentile
  X = [Distances,reverse(Distances),Distances(1)];
  Y = [Results(iInst,:,2),reverse(Results(iInst,:,4)),Results(iInst,1,1)];
  fill(X,Y,Colours(iInst),'edgecolor','none','facealpha',0.3);
  
  %50th percentile
  plot(Distances,Results(iInst,:,3),'color',Colours(iInst),'linewidth',2)
  
end
  
 plot([0 1500],[1,1],'k-','linewidth',2)
 
  
box on
set(gca,'fontsize',12,'YAxisLocation','right')
xlim([0 1200]); 
ylim([0.80 1.15])


set(gca,'ytick',0:0.1:3)
xlabel('Distance [km]')
ylabel(['MF [x ',num2str(Basis),'km value]'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%key
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

text(1010,1.08,'MLS','color','r','fontsize',14,'fontweight','bold')
text(1010,1.05,'SABER','color','k','fontsize',14,'fontweight','bold')
text(1010,1.02,'HIRDLS','color','b','fontsize',14,'fontweight','bold')


text(910,1.08,'RMW','color',[102,204,000]./255,'fontsize',14,'fontweight','bold')
text(910,1.05,'ROCI','color',[255,155,053]./255,'fontsize',14,'fontweight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%overplot histograms of storm sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StormData = load('../09Lifecycle/storms_enhanced.mat');

Axis = 0:2:800; Bandwidth = 20.*mean(diff(Axis(:)));

%produce KDP of ROCI
ROCI = StormData.Storms.OCI; ROCI = ROCI(~isnan(ROCI) & ROCI > 0);
pd.roci  = fitdist(ROCI,'kernel','Kernel','epanechnikov','bandwidth',Bandwidth);
kdp.roci = pdf(pd.roci,Axis); clear ROCI

%produce KDP of MWR
MWR = StormData.Storms.MWR; MWR = MWR(~isnan(MWR) & MWR > 0);
pd.mwr  = fitdist(MWR,'kernel','Kernel','epanechnikov','bandwidth',Bandwidth);
kdp.mwr = pdf(pd.mwr,Axis); clear MWR

%scale results to data
basis = 0.8; top = 0.937; 
kdp.roci = (kdp.roci./max(kdp.roci)).*(top-basis)+basis;
kdp.mwr = (kdp.mwr./max(kdp.mwr)).*(top-basis)+basis;

%plot as pretty plots
x  = [Axis,min(Axis)];
y1 = [kdp.roci,0]; y2 = [kdp.mwr,0];
patch(x,y1,[255,155,053]./255,'facealpha',0.66,'edgecolor',[255,155,053]./255,'linewidth',2)
patch(x,y2,[102,204,000]./255,'facealpha',0.66,'edgecolor',[102,204,000]./255,'linewidth',2)


%and binsizes
text(0,0.65,['\it{KDPs: binsize ',num2str(mean(diff(Axis))),'km, bandwidth ',num2str(Bandwidth),'km}'])


%limit line over top, to make it blacker
plot([Basis Basis],[0 10],'k-','linewidth',2)

%and add a second axis for these data
ax = gca;
ax1_pos = ax.Position; % position of first axes
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
yticklabs = [0,0.33,0.67,1];
set(ax2,'ytick',yticklabs,'yticklabel',yticklabs,'xtick',[],'fontsize',12)
ylim([0 0.5./0.2])
t = text(-0.06,0.5,'Rel. Prop. Obs.','fontsize',12); 
set(t,'rotation',90,'horizontalalignment','center')





export_fig('distance_km','-pdf','-native')