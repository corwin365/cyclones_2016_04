function [MlsData] = extract_mls_data(MatlabDay,DataDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract MLS data for a given day
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%08/MAY/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found
%2. no data on this day (not used)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the file and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the day is before or after the data record starts, give up!
if MatlabDay < datenum(2004,1,200); MlsData.Error = 1; return; end

%find data file for this day
[y,m,~] = datevec(MatlabDay);
Day = datevec2doy(datevec(MatlabDay));
y   = sprintf('%04d',y);
m   = sprintf('%02d',m);
Day = sprintf('%03d',Day);


% MLS-Aura_L2GP-Temperature_v03-30-c01_2007d004
DayFile = wildcardsearch(DataDir,['_',y,'d',Day,'.he5']);
if numel(DayFile) ==0;
  MlsData.Error = 1; 
  return;
end
clear Day y m

%check if the file exists, and load it if so
if exist(DayFile{1},'file') == 0;
  MlsData.Error = 1; 
  return;
end

Mls = get_MLS(DayFile{1},'Temperature');

MlsData.Temp   = double(Mls.L2gpValue);
MlsData.Lat    = double(Mls.Latitude)';
MlsData.Lon    = double(Mls.Longitude)';
MlsData.Height = p2h(double(Mls.Pressure)).*1000;
MlsData.Prs    = double(Mls.Pressure);
MlsData.Time   = double(Mls.Time)';

clear Mls

MlsData.Time = (MlsData.Time./86400) + datenum(1993,1,1);

%tidy arrays to match other datasets
NProfs = numel(MlsData.Lat);
MlsData.Height = repmat(MlsData.Height,1,NProfs);


MlsData.Error = 0; %it worked!

return
end
