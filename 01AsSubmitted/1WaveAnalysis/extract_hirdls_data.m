function [HirdlsData] = extract_cosmic_data(MatlabDay,DataDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract HIRDLS data for a given day
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%14/FEB/2014
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
if MatlabDay < datenum(2005,1,1); HirdlsData.Error = 1; return; end
if MatlabDay > datenum(2008,4,1); HirdlsData.Error = 1; return; end

%find data file for this day
[y,m,~] = datevec(MatlabDay);
Day = datevec2doy(datevec(MatlabDay));
y   = sprintf('%04d',y);
m   = sprintf('%02d',m);
Day = sprintf('%03d',Day);

DayFile = [DataDir,'/',y,'/',m,'/HIRDLS-Aura_L2_v07-00-20-c01_',y,'d',Day,'.he5'];
clear Day y m

%check if the file exists, and load it if so
if exist(DayFile,'file') == 0;
  HirdlsData.Error = 1; 
  return;
end

%extract data
HirdlsData.Temp   = get_HIRDLS(DayFile,'Temperature'); HirdlsData.Temp(HirdlsData.Temp == -999) = NaN;
HirdlsData.Lat    = get_HIRDLS(DayFile,'Latitude'   )';
HirdlsData.Lon    = get_HIRDLS(DayFile,'Longitude'  )';
HirdlsData.Height = get_HIRDLS(DayFile,'Altitude'   );
HirdlsData.Prs    = get_HIRDLS(DayFile,'Pressure'   );
HirdlsData.Time   = get_HIRDLS(DayFile,'Time'       )';
clear DayFile


HirdlsData.Time = (HirdlsData.Time./86400) + datenum(1993,1,1);

HirdlsData.Error = 0; %it worked!

return
end
