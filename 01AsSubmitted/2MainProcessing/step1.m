% function step1()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify which location we want to analyse GWs over
%this includes both storms and backgrounds
%this is STEP TWO of the analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.MinWind = 0; %not used, not removed in case something breaks.
Settings.TimeRange  =[datenum(2002,1,1),datenum(2013,12,31)];


%the following settings are the ones the final analysis will be done on
Settings.LatRange   = [ -50, 50];
Settings.LonRange   = [-180,180];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load storm data - first year
Settings.Storms = [LocalDataDir,'/Miscellany/stormtracks/Year.2002.ibtracs_wmo.v03r08.nc'];
Storms.Wind   = full(ncArray(Settings.Storms,'wind_wmo'));
Storms.Lat    = full(ncArray(Settings.Storms,'lat_wmo'));
Storms.Lon    = full(ncArray(Settings.Storms,'lon_wmo'));
Storms.Time   = full(ncArray(Settings.Storms,'time_wmo'));
Storms.Name   = full(ncArray(Settings.Storms,'name'));
Storms.Basin  = full(ncArray(Settings.Storms,'basin'));
Storms.Serial = full(ncArray(Settings.Storms,'storm_sn'));

%{'storm_sn';'name';'numObs';'season';'track_type';'genesis_basin';'num_basins';'basin';'wind_avg_period';'source';'time_wmo';'lat_wmo';'lon_wmo';'alt';'wind_wmo';'pres_wmo';'sub_basin';'nature_wmo';'source_wmo';'dist2land';'landfall'}

for iFile=2003:1:2014;
  Settings.Storms = [LocalDataDir,'/Miscellany/stormtracks/Year.',num2str(iFile),'.ibtracs_wmo.v03r08.nc'];
  Storms.Wind   = cat(2,Storms.Wind,full(ncArray(Settings.Storms,'wind_wmo')));
  Storms.Lat    = cat(2,Storms.Lat ,full(ncArray(Settings.Storms,'lat_wmo')));
  Storms.Lon    = cat(2,Storms.Lon ,full(ncArray(Settings.Storms,'lon_wmo')));
  Storms.Time   = cat(2,Storms.Time,full(ncArray(Settings.Storms,'time_wmo')));
  Storms.Name   = cat(2,Storms.Name,full(ncArray(Settings.Storms,'name')));  
  Storms.Basin  = cat(2,Storms.Basin,full(ncArray(Settings.Storms,'basin'))); 
  Storms.Serial = cat(2,Storms.Serial,full(ncArray(Settings.Storms,'storm_sn')));    
end; clear iFile

%conversions
%%%%%%%%%%%%%%%%%%%%%%%%

%WMO time -> Matlab time (Storms)
%days since 1858-11-17 00:00:00 (UTC)
Storms.Time = Storms.Time+datenum(1858,11,17,00,00,00);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%temporally/spatially subset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%storms
StormMinTime = nanmin(Storms.Time,[],1); StormMaxTime = nanmax(Storms.Time,[],1);
StormMinLat  = nanmin(Storms.Lat, [],1); StormMaxLat  = nanmax(Storms.Lat ,[],1);
StormMinLon  = nanmin(Storms.Lon, [],1); StormMaxLon  = nanmax(Storms.Lon ,[],1);

%remember, we're looking for anything *overlapping* the ...
%time/space window rather than *entirely contained within* it
InRange = find(StormMaxTime >= min(Settings.TimeRange) ...
             & StormMinTime <= max(Settings.TimeRange) ...
             & StormMaxLat  >= min(Settings.LatRange ) ...     
             & StormMaxLon  >= min(Settings.LonRange ) ...      
             & StormMinLat  <= max(Settings.LatRange ) ...     
             & StormMinLon  <= max(Settings.LonRange ));            
           
Storms.Wind   = Storms.Wind(  :,InRange);
Storms.Lat    = Storms.Lat(   :,InRange);
Storms.Lon    = Storms.Lon(   :,InRange);
Storms.Time   = Storms.Time(  :,InRange);
Storms.Name   = Storms.Name(  :,InRange);
Storms.Basin  = Storms.Basin( :,InRange);
Storms.Serial = Storms.Serial(:,InRange);

clear StormMinTime StormMaxTime StormMinLon StormMaxLon StormMinLat StormMaxLat InRange

save('storm_info.mat','Storms','Settings');
