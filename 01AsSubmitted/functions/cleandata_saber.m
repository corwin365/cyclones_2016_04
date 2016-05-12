function SaberData = cleandata_saber(SaberData,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clean SABER data, using a set of flags to specify exactly how
%assume all methods if no flags specified
%
%output format is that used by my analysis routines
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%09/APR/2014
%
%inputs
%---------
%
%SaberData - struct containing SABER data, as produced by downsizing routines
%varargin (optional) - flags showing methods to NOT apply, as follows:
%                      - OmitLon360_180 - converts lon from 0-360 to -180-180
%                      - OmitDateConvert - converts date format to matlab
%                      - OmitHeightConvert - puts height onto (lev x prof) grid
%                      - OmitOOR - replaces out-of-range lat and lon with NaNs
%                      - OmitIDNode - identifies asc or desc node profiles
%                      - OmitLSTConvert - converts local solar time from msec since midnight to hours
%                      - OmitTideCompute - computes the local tide phase for each point from GSWM model data
%
%
%
%outputs
%---------
%
% SaberData - cleaned data to return
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert longitude from 0-360 to -180 to +180
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitLon360_180',varargin);
  SaberData.Lon(SaberData.Lon > 180) = SaberData.Lon(SaberData.Lon > 180)-360;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fix out-of-range values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitOOR',varargin);
  SaberData.Lon(SaberData.Lon >  180) = NaN;
  SaberData.Lon(SaberData.Lon < -180) = NaN;  
  SaberData.Lat(SaberData.Lat >   90) = NaN;
  SaberData.Lat(SaberData.Lat <  -90) = NaN;    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%put height on the same grid as the rest of the data (i.e. nlevels x nprofs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if ~ismember('OmitHeightConvert',varargin);
  NProfiles = numel(SaberData.Date);  
  Height = SaberData.Temp .* NaN;
  for iProf=1:1:NProfiles; Height(:,iProf) = SaberData.Height; end; clear iProf
  SaberData.Height = Height;
  clear Height clear NProfiles
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert date format to matlab
%and remove the original time format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitDateConvert',varargin);
  NProfiles = numel(SaberData.Date);
  
  SaberData.MatlabTime = SaberData.Time .* NaN;;
  for iTime=1:1:NProfiles;
      
      Year    = floor(double(SaberData.Date(iTime))/1000.);
      Day     = double(SaberData.Date(iTime))-Year*1000.;
      Seconds = SaberData.Time(:,iTime)/1000.;
      
      SaberData.MatlabTime(:,iTime) = datenum(Year,1,Day,0,0,Seconds);
      clear Year Day Seconds
  end; clear iTime NProfiles
  
  SaberData = rmfield(SaberData,'Time');
  SaberData = rmfield(SaberData,'Date');  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%identify asc or desc modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitIDNode',varargin);

  %work out mean lat of each profile
  MeanLats = nanmean(SaberData.Lat,1);

  %work out difference between each pair of lats 
  NodeCalc = diff(MeanLats);
  NodeCalc(NodeCalc >  0) =  1; %ascending
  NodeCalc(NodeCalc <= 0) = -1; %descending
  NodeCalc(numel(MeanLats)) = NodeCalc(numel(MeanLats)-1); %assume last is same as previous
  clear MeanLats

  %assume the last profile is the same node as the previous
  NodeCalc(end) = NodeCalc(end-1);

  %put into original format
  SaberData.Nodes = NaN.*SaberData.Lat;
  
  NLevels = numel(SaberData.Lat(:,1));
  for iLevel=1:1:NLevels
    SaberData.Nodes(iLevel,:) = NodeCalc;
  end; clear iLevel
  
  clear NodeCalc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert local solar time to hours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember('OmitLSTConvert',varargin);
  SaberData.LST = SaberData.LST./1000./60./60;
end

