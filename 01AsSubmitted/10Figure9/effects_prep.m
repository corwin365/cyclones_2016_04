clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%work out aggregate effect of cyclones in each basin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load cyclone climatology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load storm data - first year
Settings.Storms = [LocalDataDir,'/Miscellany/stormtracks/Year.1990.ibtracs_wmo.v03r08.nc'];
Storms.Wind   = full(ncArray(Settings.Storms,'wind_wmo'));
Storms.Lat    = full(ncArray(Settings.Storms,'lat_wmo'));
Storms.Lon    = full(ncArray(Settings.Storms,'lon_wmo'));
Storms.Time   = full(ncArray(Settings.Storms,'time_wmo'));
Storms.Basin  = full(ncArray(Settings.Storms,'basin'))+1; %shift up 0th
Storms.NYears = 1;

%load other years
for iFile=1991:1:2014;
  Settings.Storms = [LocalDataDir,'/Miscellany/stormtracks/Year.',num2str(iFile),'.ibtracs_wmo.v03r08.nc'];
  Storms.Wind   = cat(2,Storms.Wind,full(ncArray(Settings.Storms,'wind_wmo')));
  Storms.Lat    = cat(2,Storms.Lat ,full(ncArray(Settings.Storms,'lat_wmo')));
  Storms.Lon    = cat(2,Storms.Lon ,full(ncArray(Settings.Storms,'lon_wmo')));
  Storms.Time   = cat(2,Storms.Time,full(ncArray(Settings.Storms,'time_wmo')));
  Storms.Basin  = cat(2,Storms.Basin,full(ncArray(Settings.Storms,'basin'))+1); 
  Storms.NYears = Storms.NYears+1;
end; clear iFile

clear Settings

%convert times to Matlab format
Storms.Time = Storms.Time+datenum(1858,11,17,00,00,00);

%categorise the storms
Storms.Category = Storms.Wind; Storms.Category(:) = NaN;
Storms.Category(Storms.Wind <  34) = NaN;
Storms.Category(Storms.Wind >  34) = 0;
Storms.Category(Storms.Wind >  64) = 1;
Storms.Category(Storms.Wind >  83) = 2;
Storms.Category(Storms.Wind >  96) = 3;
Storms.Category(Storms.Wind > 113) = 4;
Storms.Category(Storms.Wind > 137) = 5;

disp('Storms loaded');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load and apply cyclone climatology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load climatology
load('../8Figure7/categorised.mat');
%instrument order: MLS, SABER, HIRDLS

%produce results arrays
Storms.MF = NaN(3, ...
                numel(Storms.Category(:,1)), ...
                numel(Storms.Category(1,:)), ...
                numel(CatOut.PCs));

%and fill 'em
for iInst = 1:1:3;
  Data = squeeze(CatOut.MF(iInst,:,:));
  
  for iPC = 1:1:numel(CatOut.PCs);
    a = Storms.Category.*NaN;
    for iCat=1:1:numel(CatOut.Cats);
      Cat = CatOut.Cats(iCat);
      ThisCat = find(Storms.Category == Cat);
      a(ThisCat) = Data(iCat,iPC);
    end
    Storms.MF(iInst,:,:,iPC) = a;
  end
end
clear iInst iPC iCat a Cat ThisCat Data


disp('Categories applied');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now produce a time series for each basin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% annualise within each basin
NBasins = nanmax(Storms.Basin(:));
TimeSeries.Time = 0:0.25:365; %day-of-year
TimeSeries.Data = zeros(3,NBasins,numel(TimeSeries.Time),1); %first is instrument, last is percentiles (median only for now)

%find day-of-year of each point
[~,~,~,h,~,~] = datevec(Storms.Time);
doy = datevec2doy(datevec(Storms.Time(:)))+h(:)./24;
disp('Days-of-year computed');
clear h


%store data
for iBasin=1:1:NBasins;
  for iTime=1:1:numel(TimeSeries.Time);
    Useful = find(Storms.Basin(:) == iBasin ...
                & doy >= TimeSeries.Time(iTime) ...
                & doy < TimeSeries.Time(iTime)+mean(diff(TimeSeries.Time)));
    if numel(Useful) == 0; continue; end

    %data in this basin at this time. find sum and store
    for iInst=1:1:3;
      iPC = 3; %median
      jPC = 1; %output location
%       for iPC=1:1:numel(Storms.MF(1,1,1,:));
      Data = squeeze(Storms.MF(iInst,:,:,iPC));
      TimeSeries.Data(iInst,iBasin,iTime,jPC) = nansum(Data(Useful));
%       end
    end
  end
end
clear iBasin iTime iInst iPC Useful Data doy

%divide by number of years. which is at the top of this programme
TimeSeries.Data = TimeSeries.Data./Storms.NYears;

%and scale to daily
TimeSeries.Data = TimeSeries.Data./mean(diff(TimeSeries.Time));

disp('Time series produced');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%separately, produce a map for each basin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LatRange = [-50,50];   LatStep = 2;
LonRange = [-180,180]; LonStep = 2;

[xi,yi] = meshgrid(LonRange(1):LonStep:LonRange(2), ...
                   LatRange(1):LatStep:LatRange(2));
                 
Maps.Data = NaN(3,NBasins,numel(xi(:,1)),numel(xi(1,:)));                 
                 
for iBasin=1:1:NBasins;  
  for iInst=1:1:3;
    InBasin = find(Storms.Basin(:) == iBasin);
    x = Storms.Lon(InBasin); y = Storms.Lat(InBasin);
    z = squeeze(Storms.MF(iInst,:,:,3)); %3 is the median
    z = z(InBasin); %z(:) = 1;
    Good = ~isnan(x+y+z);
    BasinMap = bin2mat(x(Good),y(Good),z(Good),xi,yi);
    Maps.Data(iInst,iBasin,:,:) = BasinMap;
  end
end

Maps.Lon = xi; Maps.Lat = yi;
clear iBasin iInst x y z Good BasinMap

disp('Maps generated');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('effects_data','TimeSeries','Storms','Maps');