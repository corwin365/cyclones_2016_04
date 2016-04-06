function [CosmicData,OldFile] = extract_cosmic_data(MatlabDay,DataDir,OldFile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract COSMIC data for a given day
%only runs if calling a new file to previous iteration
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%14/FEB/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found
%2. no data on this day
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the file (if it's not the one from last time) and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the day is before or after the data record starts, give up!
if MatlabDay < datenum(2006,5,1); CosmicData.Error = 1; return; end

if ~isfield(OldFile,'Name'); OldFile.Name = ' '; end; %always load file afresh if old file not specified

%identify file for this month
[y,m,~]   = datevec(MatlabDay);
FileString = strcat('cosmic_atmPrf_data_',num2str(y),'_',cjw_monthname(m,'three'),'.mat');

if strcmp(FileString,OldFile.Name) == 0;
  %new file - load it up
  
  FileName = wildcardsearch(DataDir,FileString);
  if numel(FileName) == 0;  %no file
    CosmicData.Error = 1;
    return;
  end;

  
  %exists! load data
  AllCosmicData = load(FileName{1});
  clear FileName;
  
  %profs x height x [time,lat,lon,TEMP,prs]
  AllCosmicData.MatlabTime = squeeze(AllCosmicData.cosmic_dat(:,:,1))';
  AllCosmicData.Lat        = squeeze(AllCosmicData.cosmic_dat(:,:,2))';
  AllCosmicData.Lon        = squeeze(AllCosmicData.cosmic_dat(:,:,3))';
  AllCosmicData.Temp       = squeeze(AllCosmicData.cosmic_dat(:,:,4))'+273.16;
  AllCosmicData.Prs        = squeeze(AllCosmicData.cosmic_dat(:,:,5))';
  
  AllCosmicData.Height = 10000:500:60000; %from definition of the files used
  
  %duplicate out height scale
  Height2 = AllCosmicData.Temp;
  for i=1:1:numel(Height2(1,:));
    Height2(:,i) = AllCosmicData.Height;
  end
  AllCosmicData.Height = Height2; clear Height2
  
  DayScale = nanmean(AllCosmicData.MatlabTime,1);
  OldFile.Name = FileString;
  
else
  %same as last call - don't reload
  AllCosmicData = OldFile.Data;
  DayScale = nanmean(AllCosmicData.MatlabTime,1);
end
clear m y DataDir FileString

OldFile.Data = AllCosmicData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract the day we actually want
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OnThisDay = find(floor(DayScale) == MatlabDay);
if numel(OnThisDay) ==0; 
  CosmicData.Error = 2;
  return;
end
clear MatlabDay

%subset needed data
CosmicData.Time    = squeeze(nanmean(AllCosmicData.MatlabTime(   :,OnThisDay),1));
CosmicData.Temp    = AllCosmicData.Temp(   :,OnThisDay);
CosmicData.Prs     = nanmean(AllCosmicData.Prs,2);
CosmicData.Lat     = AllCosmicData.Lat(    :,OnThisDay);
CosmicData.Lon     = AllCosmicData.Lon(    :,OnThisDay);
CosmicData.Height  = AllCosmicData.Height;


%use the 30km lat and lon (Neil's suggestion)
[~,z] = min(abs(CosmicData.Height-30000));
CosmicData.Lat = squeeze(CosmicData.Lat(z,:));
CosmicData.Lon = squeeze(CosmicData.Lon(z,:));

%produce lat and lon scale
CosmicData.Lat = nanmean(CosmicData.Lat,1);
CosmicData.Lon = nanmean(CosmicData.Lon,1);

CosmicData.Error = 0; %it worked!

return
end
