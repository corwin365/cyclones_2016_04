function [Output] = format_data(SatelliteData,LatRange,LonRange,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove PWs, format data to standard resolution and subset geographically
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%14/FEB/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no profiles in region
%2. no data on this day


%remove planetary waves
[BGTemp,Perturbations] = pwremove(SatelliteData,Settings.PWs.PWLat,Settings.PWs.NPWs);


%reduce down to just region of interest, and discard unwanted variables
InLatRange = find(SatelliteData.Lat >= LatRange(1) ...
                & SatelliteData.Lat <= LatRange(2));
InLonRange = find(SatelliteData.Lon >= LonRange(1) ...
                & SatelliteData.Lon <= LonRange(2));
InRange = intersect(InLatRange,InLonRange); clear InLatRange InLonRange
NProfiles = numel(InRange);


if NProfiles == 0;
  Output.Error = 2;
  return;
end;

SatelliteData.BGTemp  = BGTemp(       :,InRange);
SatelliteData.Perturb = Perturbations(:,InRange);
SatelliteData.Height  = SatelliteData.Height(:,InRange);
if isfield(SatelliteData,'Dens'); SatelliteData.Dens   = SatelliteData.Dens(  :,InRange); end
Output.Lat     = SatelliteData.Lat(InRange);
Output.Lon     = SatelliteData.Lon(InRange);
Output.Time    = SatelliteData.Time(InRange);
clear InRange
  
  
  
%get rid of zeros in the height field and replace them with interpolated values - they break the height interp
%this may lead to some bad values, but only affects heights << tropopause
SatelliteData.Height(SatelliteData.Height == 0) = NaN;
SatelliteData.Height(SatelliteData.Height == -999) = NaN; %bad data
SatelliteData.Height = inpaint_nans(SatelliteData.Height);

%interpolate data onto new height scale
Output.BGTemperature = NaN(NProfiles,Settings.NHeightLevels);
Output.Perturb       = Output.BGTemperature;
for i=1:1:NProfiles;
  if numel(find(isfinite(SatelliteData.Temp(:,i)))) == 0; continue; end;
  warning off %spline interp spits out an error message for NaN columns that will be rejected anyway
  try %bad data won't be processed, leaving a column of NaNs which we can ignore
    Output.BGTemperature(i,:) = interp1(SatelliteData.Height(:,i),SatelliteData.BGTemp( :,i),Settings.HeightScale,'spline');
    Output.Perturb(      i,:) = interp1(SatelliteData.Height(:,i),SatelliteData.Perturb(:,i),Settings.HeightScale,'spline');
  catch;  end
  warning on
end; clear i
Output.Pressure = interp1(nanmean(SatelliteData.Height,2),SatelliteData.Prs,Settings.HeightScale);
Output.BGTemperature = Output.BGTemperature';
Output.Perturb       = Output.Perturb';
Output.Pressure = Output.Pressure';

Output.Error =0; %it worked!
return
end


