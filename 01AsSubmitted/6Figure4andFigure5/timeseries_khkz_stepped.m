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
disp(SABER.z(idx))
MLS.Results    = squeeze(MLS.Results(    idx,:,:));
SABER.Results  = squeeze(SABER.Results(  idx,:,:));
HIRDLS.Results = squeeze(HIRDLS.Results( idx,:,:));
clear idx

%discard tp and kz
MLS.Results    = MLS.Results(   3:4,:);
SABER.Results  = SABER.Results( 3:4,:);
HIRDLS.Results = HIRDLS.Results(3:4,:);

%other settings
MaxDistance  = 300e3; %m
MaxTimeRange = 15; %days
TimeStep     = 1./12; %days
TimeWindow   = 2; %days
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
NVars       = 2;%kh,kz
NInsts      = 3; %MLS,SABER,HIRDLS
NTimes      = numel(TimeScale);
Results = NaN(NVars,NInsts,NTimes,numel(Percentiles));



for iInst=1:1:NInsts;
  
  
  switch iInst
    case 1; AllData = MLS;    disp('Bootstrapping MLS');
    case 2; AllData = SABER;  disp('Bootstrapping SABER');
    case 3; AllData = HIRDLS; disp('Bootstrapping HIRDLS');
    otherwise stop
  end
  
  %remove 2008 from HIRDLS - bad data
  if iInst == 3;
    Last2007Storm = find(nanmean(Storms.Time,1) < datenum(2008,1,1));
    Last2007Storm = nanmax(Last2007Storm);
    AllData.Results(:,AllData.StormIDs(1,:) > Last2007Storm) = NaN;
  end
  
  %remove any storms too small to be worthy of the name
  Cats = AllData.StormIDs(1,:).*NaN;
  for iCat=1:1:numel(Cats);
    Cats(iCat) = Storms.Category(AllData.StormIDs(2,iCat),AllData.StormIDs(1,iCat));
  end; clear iCat
  AllData.Results(:,Cats <= -1) = NaN;
  
  
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
    
    %remove outliers
    CutOff = prctile(AllData.Results(iVar,:),[5,95]);
    Bad = find(AllData.Results(iVar,:) < CutOff(1));
    AllData.Results(iVar,Bad) = NaN;
    Bad = find(AllData.Results(iVar,:) > CutOff(2));
    AllData.Results(iVar,Bad) = NaN;  
    
    %convert to log-km^-1
    AllData.Results(iVar,:) = log10(AllData.Results(iVar,:));
    AllData.Results(~isfinite(AllData.Results)) = NaN; %infinities are a problem
    
    %store all daily means, for statistics later
    Means = NaN(NTimes,1);
        
    for iTime=1:1:NTimes;
      
      ThisTime = TimeScale(iTime);

      InTime = find(AllData.StormIDs(5,:) >= ThisTime - TimeWindow./2 ...
                  & AllData.StormIDs(5,:) <  ThisTime + TimeWindow./2 );
      clear ThisTime dt
      if numel(InTime) == 0; continue; end     
      
      Data = squeeze(AllData.Results(iVar,InTime));
      Data  = Data(~isnan(Data));      
      if nansum(Data(:)) == 0; continue; end      
      
      NPoints = numel(Data);
      RandData = NaN(NStraps,1);
      for iStrap = 1:1:NStraps;
        Indices = randi(NPoints,NPoints,1);
        RandData(iStrap) = mean(Data(Indices));
        
      end; clear iStrap Indices 
      
      Results(iVar,iInst,iTime,:) = prctile(RandData,Percentiles);
      
      clear RandData
      
    end; clear iTime

    
  end; clear iVar
  clear AllData Data InTime
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
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.15], [0.1 0.1], [0.1 0.1]);

Colours = 'rkb';

Instruments = {'MLS','SABER','HIRDLS'};
Letters = 'adbecf';
figcount = 0;



disp('Plotting time');
for iInst=1:1:3;
  for iVar=1:1:2;
    
    %figure positioning etc
    figcount = figcount+1;
    subplot(3,2,figcount)
    box on
    Colour = Colours(iInst);
    hold on
  
    %linearly detrend
    x = squeeze(Results(iVar,iInst,:,3))';
    t = TimeScale;
    p= polyfit(t,x,1); p(2) = 0;
    for i=1:1:5;
      Results(iVar,iInst,:,i) = squeeze(Results(iVar,iInst,:,i))-squeeze(polyval(p,t))';
    end
    clear x t p i
  
   %mean properties
   Data   = squeeze(Results(iVar,iInst,:,3));
   TheMedian = prctile(Data,50);
   SE66 = prctile(Data,[18,82]);
   SE95 = prctile(Data,[2.5,97.5]);
   B90  = prctile(Data,[10,95]);
   
   %shade places where the estimate is outside the confidence bound
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   a = squeeze(Results(iVar,iInst,:,3))';
   Conf  = B90(2);%max(SE95);%;TheMean+2.*TheStd;
   Conf2 = -99e99;%min(SE95);%TheMean-2.*TheStd;
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
   
   
   
   
   
   
   %finish up
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if     iVar == 1; ylabel('log_{10}(k_h [m^{-1}])')
   elseif iVar == 2; ylabel('log_{10}(k_z [m^{-1}])')
   end
   
   if iInst == 3;
     xlabel('Days since storm centre');
   elseif iInst == 1;
     set(gca,'xaxis','top');xlabel('Days since storm centre');
   else set(gca,'xtick',[]); end
   
   xlim([min(TimeScale) max(TimeScale)])
%    axis tight
   
   YBottom = min(squeeze(Results(iVar,iInst,:,1)))-0.001./iVar;
   YTop    = max(squeeze(Results(iVar,iInst,:,5)))+0.001./iVar;
   ylim([YBottom YTop])
   
   set(gca,'FontSize',14,'tickdir','out')
   
   
   
   %text labelling
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   text(-14.5,YBottom+(YTop-YBottom).*0.9,['(',Letters(figcount),') ',Instruments{iInst}],'FontSize',18)
   
   
   %wavelength axis
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %figure out values
   ax = gca;
   ticklocs = ax.YTick;
   
   if iVar == 1; tickvals = round((1e-3./10.^(ax.YTick))./10).*10;
   else          tickvals = round((1e-3./10.^(ax.YTick)).*10)./10;
   end
   
   %overplot second axis set
   ax1_pos = ax.Position; % position of first axes
   ax2 = axes('Position',ax1_pos,...
     'XAxisLocation','top',...
     'YAxisLocation','right',...
     'Color','none');
   set(ax2,'ytick',ticklocs)
   set(ax2,'ytick',ticklocs,'yticklabel',tickvals)
   ylim([YBottom YTop])
   
   if iVar == 1; ylabel('Horizontal \lambda [km]');
   else          ylabel('Vertical \lambda [km]');
   end
   
   set(ax2,'xtick',[])
   %   stop
   
   
   
   
   
   drawnow
  end
end

export_fig('timeseries_khkz','-pdf','-native')