clear all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
MLS    = load([LocalDataDir,'/corwin/storms/results/storms_MLS_close.mat']);    MLS    = MLS.MLS.Results;
SABER  = load([LocalDataDir,'/corwin/storms/results/storms_SABER_close.mat']);  SABER  = SABER.SABER.Results;
HIRDLS = load([LocalDataDir,'/corwin/storms/results/storms_HIRDLS_close.mat']); HIRDLS = HIRDLS.HIRDLS.Results;

% MLS = SABER;

%discard profiles outside our chosen range
MaxDistance = 300e3;
MaxTimeRange = 15; %days
z = SABER.z./1000; %same for all


%time windowing
InnerWindow = [-1,2];
OuterWindow = 5;% minimum separation - we'll random-sample the same number of point as the inner windpw


disp('Data loaded and height level selected')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the category of all "storms", so we can exclude the ones so small
%they don't really count
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load storm data
load('../2MainProcessing/storm_info.mat'); clear Settings

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
%find the distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bootstrapping
NLevels     = numel(z); 
Percentiles = [2.5,18,50,82,97.5];
NVars       = 4;
NInsts      = 3; %MLS,SABER,HIRDLS
Results     = NaN(NVars,NInsts,NLevels,numel(Percentiles)); 
NStraps     = 1000; %bootstrap samples



for iInst=1:1:NInsts;
  
  switch iInst
    case 1; AllData = MLS;    disp('Processing MLS');
    case 2; AllData = SABER;  disp('Processing SABER');
    case 3; AllData = HIRDLS; disp('Processing HIRDLS');
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
  AllData.StormIDs = AllData.StormIDs( :,Good);
  AllData.Results  = AllData.Results(:,:,Good);
  
  %discard out-of-timerange profiles
  Good = find(abs(AllData.StormIDs(5,:)) <= MaxTimeRange);
  if numel(Good) == 0; continue; end
  AllData.StormIDs = AllData.StormIDs( :,Good);
  AllData.Results  = AllData.Results(:,:,Good);
  
  for iVar = 1:1:NVars;
    
    %scale MF
    if iVar == 2;
      AllData.Results(:,iVar,:) = AllData.Results(:,iVar,:).*1000;
    end
    
    %remove any storms too small to be worthy of the name
    Cats = AllData.StormIDs(1,:).*NaN;
    for iCat=1:1:numel(Cats);
      Cats(iCat) = Storms.Category(AllData.StormIDs(2,iCat),AllData.StormIDs(1,iCat));
    end; 
    AllData.Results(:,:,Cats == -1) = NaN;
        

    
    %and the storm period
    Period.Storm = find(AllData.StormIDs(5,:) >= InnerWindow(1) ...
                      & AllData.StormIDs(5,:) <= InnerWindow(2));
                    
    %find the background period
    Period.BG    = find(abs(AllData.StormIDs(5,:)) >= OuterWindow);
    %and random sample a subset of it
    Period.BG = Period.BG(randi(numel(Period.BG),numel(Period.Storm),1));
    
    %loop over levels
    for iLevel=1:1:NLevels;
      Data = squeeze(AllData.Results(iLevel,iVar,:));
      
      %remove any NaNs
      Good = find(~isnan(Data));
      if numel(Good) == 0; continue; end
      
      %create storage array
      Straps = zeros(NStraps,2); %2 is storms and BG
      
      %process storms
      Points = Data(intersect(Period.Storm,Good));
      for i=1:1:NStraps; 
        Random = randi(numel(Points),numel(Points),1);
        Straps(i,1) = nanmean(Points(Random));
      end; clear i Random
      
      %process BG
      Points = Data(intersect(Period.BG,Good));
      for i=1:1:NStraps; 
        Random = randi(numel(Points),numel(Points),1);
        Straps(i,2) = nanmean(Points(Random));
      end; clear i Random      
      clear Points
      
      %divide the one by the other to produce a % deviation
      Straps = ((Straps(:,1)./Straps(:,2)).*100)-100;
      
      %characterise the two datasets
      Results(iVar,iInst,iLevel,:) = prctile(Straps,Percentiles);
    
    end; clear iLevel Straps


  end; clear iVar
  clear AllData
end; clear iInst

clear Cats Data Good HIRDLS InnerWindow MaxLatRange MaxLonRange MaxTimeRange
clear MLS OuterWindow Period SABER Storms 
clear NInsts NLevels NVars



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.09 0.03], [0.1 0.04], [0.14 0.07]);

Count = 0;
Colours = 'rkb';
Letters = 'abcdefghijklmnopqrstuvwxyz';

%fill plot y-axis
y = [z,reverse(z)];

for iInst=1:1:3;
  for iVar=1:1:4;
    Count = Count+1;
    subplot(3,4,Count)
    hold on
    
    %highlight level of other analyses
    plot([-1,1].*999,[1,1].*30,'-','color',[1,1,1].*0,'linewidth',2);    
    
    
    %find data
    St = squeeze(Results(iVar,iInst,:,:));
   
    St = inpaint_nans(St);

    
    %show data bounds
    Colour = Colours(iInst);
    x = [St(:,1);reverse(St(:,5))];
    patch(x,y,Colour,'edgecolor','none','facealpha',0.4)
    x = [St(:,2);reverse(St(:,4))];
    patch(x,y,Colour,'edgecolor','none','facealpha',0.4)
    plot(St(:,3),z,'color',Colour,'linewidth',2)
    hold on
    
    
    %OK. polish panels
    box on
    
    

    %the equality line
    plot([0 0],[-999 999],'k-')
    
    %y-axes and horizontal lines
    if     iVar == 1; set(gca,'yaxis','left');  ylabel('Altitude [km]','FontSize',14);
    elseif iVar == 4; set(gca,'yaxis','right'); ylabel('Altitude [km]','FontSize',14);
    else              set(gca,'ytick',[]);
    end
    for i=min(z):5:max(z); if i==30; continue; end; plot([-1,1].*999,[1,1].*i,'--','color',[1,1,1].*0.5); end
    
    
    %x-axes
    switch iVar
      case 1; xlabel('\DeltaT'' [%]','FontSize',14);
      case 2; xlabel('\DeltaMF [%]','FontSize',14);
      case 3; xlabel('\Deltak_h [%]','FontSize',14);
      case 4; xlabel('\Deltak_z [%]','FontSize',14);
    end

    %limits and labelling
    switch iVar
      case 1; Limits = [98 112]-100;
      case 2; Limits = [95 130]-100;
      case 3; Limits = [92 111]-100;
      case 4; Limits = [97 103]-100;
      otherwise; stop; 
    end
    xlim(Limits);ylim([20 60])
  
    %panel labels
    text(min(Limits)+0.75*range(Limits),57,['(',Letters(Count),')'],'fontsize',22)  
    
    
    %instrument label
    xpos = min(Limits)-0.5*range(Limits);
    if iVar == 1;
      if     iInst == 1; h = text(xpos,mean(z),'MLS',   'FontSize',30,'horizontalalign','center');
      elseif iInst == 2; h = text(xpos,mean(z),'SABER', 'FontSize',30,'horizontalalign','center');
      elseif iInst == 3; h = text(xpos,mean(z),'HIRDLS','FontSize',30,'horizontalalign','center');
      end
    end
    set(h,'rotation',90)
    
    drawnow
    
  end
end


export_fig('height_delta','-pdf','-native')
