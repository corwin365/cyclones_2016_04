clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time series of composite-storm effects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
MLS    = load([LocalDataDir,'/corwin/storms/results/storms_MLS_close.mat']);    MLS    = MLS.MLS.Results;
SABER  = load([LocalDataDir,'/corwin/storms/results/storms_SABER_close.mat']);  SABER  = SABER.SABER.Results;
HIRDLS = load([LocalDataDir,'/corwin/storms/results/storms_HIRDLS_close.mat']); HIRDLS = HIRDLS.HIRDLS.Results;


%select a height level
idx = 6; 
SABER.z(idx)
% stop
MLS.Results    = squeeze(MLS.Results(    idx,:,:));
SABER.Results  = squeeze(SABER.Results(  idx,:,:));
HIRDLS.Results = squeeze(HIRDLS.Results( idx,:,:));
clear idx

%other settings
MaxDistance  = 300e3; %m
MaxTimeRange = 15; %days
TimeStep     = 1./12; %days
TimeWindow   = 1.; %days
TimeScale    = -MaxTimeRange:TimeStep:MaxTimeRange;

disp('Data loaded and height level selected')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the category of all "storms", so we can exclude the ones so small
%they don't really count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load storm data
load('../../02PhaseTwo/storm_info.mat'); clear Settings

% convert to a category.
%use knots, as that's what the data are in 
%cat5: > 137 kts
%cat4: > 113 kts
%cat3: >  96 kts
%cat2: >  83 kts
%cat1: >  64 kts
%tropical storm: > 34kts

Storms.Category = Storms.Wind; Storms.Category(:) = NaN;
Storms.Category(Storms.Wind <  34) = -1; 
Storms.Category(Storms.Wind >  34) = 0;
Storms.Category(Storms.Wind >  64) = 1;
Storms.Category(Storms.Wind >  83) = 2;
Storms.Category(Storms.Wind >  96) = 3;
Storms.Category(Storms.Wind > 113) = 4;
Storms.Category(Storms.Wind > 137) = 5;





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bootstrapping and noise floor computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bootstrapping
Percentiles = [5,32,50,68,95];
NStraps     = 1000; %number of random trials. main speed bottleneck, set to 1000 for publication figures
NVars       = 2;%Tp,MF
NInsts      = 3; %MLS,SABER,HIRDLS
NTimes      = numel(TimeScale);
Results = NaN(NVars,NInsts,NTimes,numel(Percentiles));

%noise floor
Floors = NaN(NVars,NInsts,2); %mean, stdev


for iInst=1:1:NInsts;
%   for iInst=3:1:3;
  
  switch iInst
    case 1; AllData = MLS;    disp('Bootstrapping MLS');
    case 2; AllData = SABER;  disp('Bootstrapping SABER');
    case 3; AllData = HIRDLS; disp('Bootstrapping HIRDLS');
    otherwise; stop
  end
  
  %remove 2008 from HIRDLS - bad data
  if iInst == 3;
    Last2007Storm = find(nanmean(Storms.Time,1) < datenum(2008,1,1));
    Last2007Storm = nanmax(Last2007Storm);
    AllData.Results(:,AllData.StormIDs(1,:) > Last2007Storm) = NaN;
  end  

  %discard out-of-range profiles
  Good = find(abs(AllData.StormIDs(6,:)) <= MaxDistance);
  if numel(Good) == 0; continue; end
  AllData.StormIDs = AllData.StormIDs(:,Good);
  AllData.Results  = AllData.Results( :,Good);
  
  %discard out-of-timerange profiles
  Good = find(abs(AllData.StormIDs(5,:)) <= MaxTimeRange);
  if numel(Good) == 0; continue; end
  AllData.StormIDs = AllData.StormIDs(:,Good);
  AllData.Results  = AllData.Results( :,Good);

  for iVar = 1:1:NVars;
    
    %scale MF
    if iVar == 2;
      AllData.Results(iVar,:) = AllData.Results(iVar,:).*1000;
    end
    
    %remove any storms too small to be worthy of the name
    Cats = AllData.StormIDs(1,:).*NaN;
    for iCat=1:1:numel(Cats);
      Cats(iCat) = Storms.Category(AllData.StormIDs(2,iCat),AllData.StormIDs(1,iCat));
    end; clear iCat
    AllData.Results(:,Cats == -1) = NaN;
          
    %remove any columns with NaNs in
    Good = sum(AllData.Results,1);
    Good = find(~isnan(Good) == 1);
    AllData.Results = AllData.Results(:,Good);
    AllData.StormIDs = AllData.StormIDs(:,Good);

    %store all daily means, for statistics later
    Means = NaN(NTimes,1);
    
    for iTime=1:1:NTimes;
      
      ThisTime = TimeScale(iTime);
      dt       = mean(diff(TimeScale));

      InTime = find(AllData.StormIDs(5,:) >= ThisTime - TimeWindow./2 ...
                  & AllData.StormIDs(5,:) <  ThisTime + TimeWindow./2 );
      clear ThisTime dt
      if numel(InTime) == 0; continue; end     
      
      Data = squeeze(AllData.Results(iVar,InTime));
      Data  = Data(~isnan(Data));      
      if nansum(Data(:)) == 0; continue; end
            
      Means(iTime) = nanmean(Data(:));
      
      NPoints = numel(Data);
      RandData = NaN(NStraps,1);
      for iStrap = 1:1:NStraps;
        Indices = randi(NPoints,NPoints,1);
        RandData(iStrap) = mean(Data(Indices));
        
      end; clear iStrap Indices 
      
      Results(iVar,iInst,iTime,:) = prctile(RandData,Percentiles);
      
      clear RandData
      
    end; clear iTime
    
    Floors(iVar,iInst,1) = nanmean(Means);
    Floors(iVar,iInst,2) = nanstd(Means);
    
  end; clear iVar
  clear AllData
end; clear iInst

clear NStraps NTimes


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.1], [0.1 0.1], [0.1 0.1]);

Colours = 'rkb';

Instruments = {'MLS','SABER','HIRDLS'};
Letters = 'adbecf';
figcount = 0;
for iInst=1:1:3;
  for iVar=1:1:2;
  
  figcount = figcount+1;
  subplot(3,2,figcount)
  box on
  Colour = Colours(iInst);
  hold on
  
  %linearly detrend
  x = squeeze(Results(iVar,iInst,:,3))';
  t = TimeScale;
  p= polyfit(t,x,1); p(2) = 0;
%   stop
  for i=1:1:5;
    Results(iVar,iInst,:,i) = squeeze(Results(iVar,iInst,:,i))-squeeze(polyval(p,t))';
  end
  clear x t p i
   
  
  %mean properties
  Data   = squeeze(Results(iVar,iInst,:,3));
  TheMedian = prctile(Data,50);
  SE66 = prctile(Data,[18,82]);
  SE95 = prctile(Data,[0,95]);
  B90 = prctile(Data,90);
  
  %shade places where the estimate is outside the confidence bound
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  a = squeeze(Results(iVar,iInst,:,3))';
  Conf  = B90;
  Conf2 = 0;
  Above = find(a > Conf | a < Conf2);
  if numel(Above) >0;
    
    for i=1:1:numel(Above)-1;
      if Above(i) > numel(TimeScale); continue; end
      x1 = TimeScale(Above(i))-1*TimeStep;
      x2 = TimeScale(Above(i))+1*TimeStep;
      x  = [x1,x2,x2,x1,x1];
      y = [-1,-1,1,1,-1].*999;
      fill(x,y,[204,255,209]./255,'edgecolor','none')
    end; clear i
  end
  
  
  %lines etc
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for iDay=-10:5:10; plot([1,1].*iDay,[-999,999],'--','linewidth',1,'color',[1,1,1].*0.7); end  
  plot([0,0],[-999,999],'k-','linewidth',1)  
  
  %significance thresholds
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SigColour = [1,1,1].*0;
  plot([min(TimeScale),max(TimeScale)],[1,1].*TheMedian,'-', 'linewidth',2,'color',SigColour)
  plot([min(TimeScale),max(TimeScale)],[1,1].*SE66(1),'--','linewidth',1,'color',SigColour)
  plot([min(TimeScale),max(TimeScale)],[1,1].*SE95(1),'-.','linewidth',1,'color',SigColour)
  plot([min(TimeScale),max(TimeScale)],[1,1].*SE66(2),'--','linewidth',1,'color',SigColour)
  plot([min(TimeScale),max(TimeScale)],[1,1].*SE95(2),'-.','linewidth',1,'color',SigColour)
  
  %bootstrap confidence intervals and data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %outer
  a = squeeze(Results(iVar,iInst,:,1))';
  b = squeeze(Results(iVar,iInst,:,5))';
  x = [TimeScale,reverse(TimeScale),TimeScale(1)];
  y = [a,reverse(b),a(1)]; %y = smooth(y,5);
  fill(x,y,Colour,'Facealpha',0.2,'edgecolor','none')
  %inner
  a = squeeze(Results(iVar,iInst,:,2))';
  b = squeeze(Results(iVar,iInst,:,4))';
  x = [TimeScale,reverse(TimeScale),TimeScale(1)];
  y = [a,reverse(b),a(1)]; %y = smooth(y,5);
  fill(x,y,Colour,'Facealpha',0.4,'edgecolor','none')
  
  % %centre
  Data = squeeze(Results(iVar,iInst,:,3));
  plot(TimeScale,Data,'-' ,'color',Colour,'linewidth',1)
  
  
  %percent increase from median to max?
%   -100+(max(Data(:))./TheMedian).*100

  %text labelling
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  YLimits = [min(Data(:))-0.005 max(Data(:)).*1.01+0.02*(iVar-1)];
  YText = (min(YLimits)+0.9*range(YLimits));
  text(-14.5,YText,['(',Letters(figcount),') ',Instruments{iInst}],'FontSize',18)
  
  %finish up
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if     iVar == 1; ylabel('T'' [K]')
  elseif iVar == 2; ylabel('MF [mPa]')
  end
  
  set(gca,'xminortick','on')
  if iInst == 3;
    xlabel('Days since storm centre');
  elseif iInst == 1;
    set(gca,'xaxis','top');xlabel('Days since storm centre');
  else set(gca,'xtick',[-999 999]); end
  
  xlim([min(TimeScale) max(TimeScale)])
  ylim(YLimits)
  set(gca,'FontSize',14,'tickdir','out')
  drawnow
  end
end


% export_fig('timeseries_tpmf','-pdf','-native')