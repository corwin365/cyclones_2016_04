clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%violin plots of GW effects by storm category, all datasets
%also store the GW T' and MF associated with each category, for later use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
MLS    = load([LocalDataDir,'/corwin/storms/results/storms_MLS_close.mat']);    MLS    = MLS.MLS.Results;
SABER  = load([LocalDataDir,'/corwin/storms/results/storms_SABER_close.mat']);  SABER  = SABER.SABER.Results;
HIRDLS = load([LocalDataDir,'/corwin/storms/results/storms_HIRDLS_close.mat']); HIRDLS = HIRDLS.HIRDLS.Results;

% MLS = SABER; %temporary

%compress down to a single height level
[~,idx] = min(abs(SABER.z-30000));
MLS.Results    = squeeze(MLS.Results(   idx,:,:));
HIRDLS.Results = squeeze(HIRDLS.Results(idx,:,:));
SABER.Results  = squeeze(SABER.Results( idx,:,:));
clear idx

disp('GW data loaded and prepped')

%discard profiles outside our chosen range
MaxDistance  = 300e3; %m
MaxTimeRange = 15; %days
InnerWindow = [-1,2];
OuterWindow = 5;%minimum time-from-peak - will select same number as in innerwindow

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loaded corresponding storm data
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%link up the two datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Percentiles = [2.5,18,50,82,97.5];%5:1:95;

DistribData = NaN(3,4,9,numel(Percentiles)); %instruments, variables, categories
NPoints   = NaN(3,4,9);
Means = NPoints;

%axes for the kernel distribution functions
for i =1:1:4; ViolinAxis(i,:) = 0.000:0.001:1; end; clear i
ViolinAxis(1,:) = 0+ViolinAxis(1,:).*5; %Tprime
ViolinAxis(2,:) = 0+ViolinAxis(2,:).*20;%MF
ViolinAxis(3,:) = -8.2+ViolinAxis(3,:).*5.0; %log-kh
ViolinAxis(4,:) = -5+ViolinAxis(4,:).*2; %log-kz
ViolinData  = NaN(3,4,9,numel(ViolinAxis(1,:)));

for iDataSet=1:1:3;
  
  switch iDataSet
    case 1; Data = MLS;    disp('Processing MLS')
    case 2; Data = SABER;  disp('Processing SABER')      
    case 3; Data = HIRDLS; disp('Processing HIRDLS')
    otherwise; stop; 
  end
   
  %scale variables appropriately
  Data.Results(2,:) = Data.Results(2,:).*1000; %Pa->mPa
  
  %reduce down to time/space window desired

  %discard out-of-range profiles
  Good = find(abs(Data.StormIDs(6,:)) <= MaxDistance);
  if numel(Good) == 0; continue; end
  Data.StormIDs = Data.StormIDs(:,Good);
  Data.Results  = Data.Results( :,Good);
  
  %discard out-of-timerange profiles
  Good = find(abs(Data.StormIDs(5,:)) <= MaxTimeRange);
  if numel(Good) == 0; continue; end
  Data.StormIDs = Data.StormIDs(:,Good);
  Data.Results  = Data.Results( :,Good);
  
  %find the background period
  Period.BG    = find(abs(Data.StormIDs(5,:)) >= OuterWindow);
  
  %and the storm period
  Period.Storm = find(Data.StormIDs(5,:) >= InnerWindow(1) ...
                    & Data.StormIDs(5,:) <= InnerWindow(2));
  
  %find the category of each point in the storm array
  Cats = Data.StormIDs(1,:).*NaN;
  for iCat=1:1:numel(Cats);
    Cats(iCat) = Storms.Category(Data.StormIDs(2,iCat),Data.StormIDs(1,iCat));
  end; clear iCat
  
  %and produce a histogram of each category from this
  for iCat=1:1:9;
    
    if iCat <= 7; %real category: data for this
      ThisCat = intersect(find(Cats == iCat-2),Period.Storm);
    elseif iCat == 8; %all categories
      ThisCat = intersect(find(Cats > -1),Period.Storm);
    elseif iCat == 9; %background
      ThisCat = intersect(find(Cats > -1),Period.BG);
      %select random number of obs corresponding to real storm obs
      ThisCat = ThisCat(randi(numel(ThisCat),NPoints( iDataSet,iVar,8),1));
    end

    for iVar=1:1:4;
      VarData = Data.Results(iVar,ThisCat);
      VarData = VarData(~isnan(VarData));
      
      if iVar >= 3; VarData = log10(VarData); end %wavelengths logged

      
      %find the numerical distribution
      NPoints( iDataSet,iVar,iCat) = numel(VarData);
      Means(   iDataSet,iVar,iCat) = mean(VarData);
      DistribData(iDataSet,iVar,iCat,:)= prctile(VarData,Percentiles);
      
      %find a kernel density estimate
      bandwidth = mean(diff(ViolinAxis(iVar,:))).*10;
      pd = fitdist(VarData(:),'kernel','Kernel','epanechnikov','bandwidth',bandwidth);
      y = pdf(pd,ViolinAxis(iVar,:));
      ViolinData( iDataSet,iVar,iCat,:) = y;
    end
    
  end; clear iCat iVar VarData ThisCat pd y

  Count = 0;

  clear Data
end; clear iDataSet

clear HIRDLS SABER MLS NTimes Storms
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
subplot = @(m,n,p) subtightplot (m, n, p, [0.18 0.15], 0.1, 0.1);
Vars = {'(a) T''','(b) MF','(c) k_h','(d) k_z'};
Colours = 'rkb';
Cats = -1:1:5;
set(gcf,'color','w')


XLabels = {'TD','TS','Cat1','Cat2','Cat3','Cat4','Cat5','All','BG','All','BG','All','BG'};
vars = {'T'' [K]','MF [mPa]','log_{10}(k_h [m^{-1}])','log_{10}(k_z [m^{-1}])'};


% iDataSet =2; %do SABER first

%labels, plot bounds, tra la la
NumberLoc = [0.,-0.3,-8.27,-4.36]; %count annotation
Tops    = [ 4,4,-5.1,-3.27];
Bottoms = [-0.1,-0.5,-8.4,-4.40]; 
KDFUnits = {'mK','\muPa','x10^{-3}log_{10}(m^{-1})','x10^{-3}log_{10}(m^{-1})'};


for iVar=1:1:4
  subplot(2,2,iVar)

  hold on; box on
  warning off
  set(gca,'xtick',-1:1:11,'xticklabels',XLabels,'XTickLabelRotation',90)
  set(gca,'FontSize',14)
  warning on
  xlim([-0.5 7.5+4]) %note this plots but does not remove TD category
  ylim([Bottoms(iVar),Tops(iVar)])
  ylabel(vars{iVar})
  
  
  %title
  x = -.5;
  y = Bottoms(iVar)+(Tops(iVar)-Bottoms(iVar)).*1.10;
  text(x,y,Vars{iVar},'FontSize',24)  
  
  %setup plots
   %shade special cases region - SABER
  patch([5.5 7.5 7.5 5.5 5.5],      ...
        [-999 -999 999 999 -999], ...
        [-.1 -.1 -.1 -.1 -.1],    ... %negative z, i.e. behind the data
        [255,224,224]./255,'edgecolor','none');
  text(6.5,1.05*(Tops(iVar)-Bottoms(iVar))+Bottoms(iVar),'MLS','color','r','FontWeight','bold','HorizontalAlignment','center')     
  %shade special cases region - MLS
  patch([7.5 9.5 9.5 7.5 7.5],      ...
        [-999 -999 999 999 -999], ...
        [-.1 -.1 -.1 -.1 -.1],    ... %negative z, i.e. behind the data
        [1,1,1].*0.9,'edgecolor','none');
  text(8.5,1.05*(Tops(iVar)-Bottoms(iVar))+Bottoms(iVar),'SABER','color','k','FontWeight','bold','HorizontalAlignment','center')       
  %shade special cases region - HIRDLS
  patch([9.5 11.5 11.5 9.5 9.5],      ...
        [-999 -999 999 999 -999], ...
        [-.1 -.1 -.1 -.1 -.1],    ... %negative z, i.e. behind the data
        [204,229,255]./255,'edgecolor','none');
  text(10.5,1.05*(Tops(iVar)-Bottoms(iVar))+Bottoms(iVar),'HIRDLS','color','b','FontWeight','bold','HorizontalAlignment','center')       

  

  
  for jCat=0:1:11; %including -1 will plot TDs
    
    if jCat <= 7; 
      jDataSet = 1; %MLS
      Colour = 'r';
      Fill   = [255,204,204]./255;
      iCat = jCat;
    elseif jCat <= 9;
      jDataSet = 2; %SABER
      Colour = 'k';  
      Fill   = [192,192,192]./255;      
      iCat = jCat-2;
    elseif jCat <= 11;
      jDataSet = 3; %HIRDLS
      Colour = 'b';
      Fill   = [153,204,255]./255;         
      iCat = jCat-4;
    end
    Data  = squeeze(DistribData(jDataSet,iVar,:,:));
    Viola = squeeze(ViolinData( jDataSet,iVar,:,:));
    Data  = Data( iCat+2,:);
    Viola = Viola(iCat+2,:);
   
    %violin settings
    %%%%%%%%%%%%%%%%%%%%%%%%

    
    %short variable names to keep concise
    MW = 0.45; %maximum width
    
    %compute violin widths
    Width = Viola./max(Viola(:)).*MW;
       
    %hence, x-values
    X =jCat+[Width,0-reverse(Width)];
    
    %y-values are the percentiles
    Y = [ViolinAxis(iVar,:),reverse(ViolinAxis(iVar,:))];
    
    %remove empty ends
    Full = find([Viola,reverse(Viola)] ~= 0);
    X = X(Full); Y =Y(Full);
    


    %plot violin!
    patch(X,Y,Fill,'edgecolor',Colour,'linewidth',0.5);%,'clipping','off')
    
    %now, overplot interesting percentiles as a box plot
    MW = 0.2;
    LW = 1;
    %horizontal lines
    plot(jCat+[-1,1].*MW.*1.0,[1,1].*Data(1),'k-','linewidth',LW)        
    plot(jCat+[-1,1].*MW.*0.8,[1,1].*Data(2),'k-','linewidth',LW)    
  
    plot(jCat+[-1,1].*MW.*0.8,[1,1].*Data(4),'k-','linewidth',LW)    
    plot(jCat+[-1,1].*MW.*1.0,[1,1].*Data(5),'k-','linewidth',LW)    
    %vertical lines
    plot(jCat+[1,1].*MW.*0.8,[Data(2),Data(4)],'k-','linewidth',LW)    
    plot(jCat-[1,1].*MW.*0.8,[Data(2),Data(4)],'k-','linewidth',LW)    
    plot(jCat+[0,0],         [Data(1),Data(2)],'k-','linewidth',LW)       
    plot(jCat+[0,0],         [Data(4),Data(5)],'k-','linewidth',LW)  
    %emphasise median
    plot(jCat,Data(3),'kx','linewidth',2,'markersize',10)
    
    %number of points
    text(jCat,NumberLoc(iVar), ...
         num2str(NPoints(jDataSet,iVar,iCat+2)), ...
         'horizontalalignment','center','fontsize',7)

    
    %line up BG and FG medians
    if iCat == 6
      DS = squeeze(DistribData(jDataSet,iVar,[0,1]+iCat+2,3));
      plot([0,1]+jCat,DS,'k:','linewidth',2)
    end
    
       
    drawnow
       
  end; clear jCat
  

  
  %add a label describing the KDF used
  KDFText1 = ['bin ',num2str(mean(diff(ViolinAxis(iVar,:))).*1000),KDFUnits{iVar}];
  KDFText2 = ['bw ',num2str(20.*mean(diff(ViolinAxis(iVar,:))).*1000),KDFUnits{iVar}]; 
  text(7.5+4,Bottoms(iVar)+(Tops(iVar)-Bottoms(iVar)).*-0.25,...%1.13, ...
      ['\it{KDP properties: ',KDFText1,', ',KDFText2,'}'], ...
      'horizontalalign','right')
  
    set(gca,'tickdir','out')
  if iVar > 2;
 %wavelength axis
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %figure out values
  ax = gca;
  ticklocs = ax.YTick;
  
  if iVar == 3; tickvals = round((1e-3./10.^(ax.YTick))./10).*10;
  else          tickvals = round((1e-3./10.^(ax.YTick)).*10)./10;
  end
  
  %overplot second axis set
  ax1_pos = ax.Position; % position of first axes
  ax2 = axes('Position',ax1_pos,...
             'XAxisLocation','top',...
             'YAxisLocation','right',...
             'Color','none','tickdir','out');
  set(ax2,'ytick',ticklocs)
  set(ax2,'ytick',ticklocs,'yticklabel',tickvals)  
  ylim([Bottoms(iVar) Tops(iVar)])
  
  if iVar == 3; ylabel('Horizontal \lambda [km]');
  else          ylabel('Vertical \lambda [km]');
  end  
  
  set(ax2,'xtick',[])
  end
  
end

% export_fig('categories','-pdf','-native')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store key data for later use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


CatOut.MF   = squeeze(DistribData(:,2,:,:));
CatOut.PCs  = Percentiles;
CatOut.Cats = Cats; 
save('categorised.mat','CatOut');

