function step3(Instrument)

% Instrument = 'SABER'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the seasonal cycle of the satellite data
%this is STEP THREE of the analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Instrument = Instrument;
Settings.InFile     = [LocalDataDir,'/corwin/storms2/background/',Settings.Instrument,'.mat'];
Settings.OutFile    = [LocalDataDir,'/corwin/storms2/climatology/',Settings.Instrument,'.mat'];
Settings.SmoothSize = 30; %days' smoothing to apply to the seasonal trend as inst noise high


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data and prep arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LargeScale = load(Settings.InFile); 

%transpose old settings
Settings.TimeStep = LargeScale.Settings.TimeStep;
Settings.LatRange = LargeScale.Settings.LatRange;
Settings.LonRange = LargeScale.Settings.LonRange;
Settings.LatStep  = 2;
Settings.LonStep  = 2;


Seas.Lon  = Settings.LonRange(1):Settings.LonStep:Settings.LonRange(2);
Seas.Lat  = Settings.LatRange(1):Settings.LatStep:Settings.LatRange(2);
Seas.Time = 1:Settings.TimeStep:365;
Seas.z    = LargeScale.Settings.z;
Seas.Tp   = NaN(numel(Seas.Lon),numel(Seas.Lat),numel(Seas.z),numel(Seas.Time));
Seas.MF   = Seas.Tp;
Seas.Kh   = Seas.Tp;
Seas.Kz   = Seas.Tp;

LargeScale = LargeScale.BGGrid; %simplifies typing later


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the seasonality, and store
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

textprogressbar('Computing seasonality ')
for iLevel=1:1:numel(LargeScale.z);
  textprogressbar(100.*iLevel./numel(LargeScale.z))
  
  %extract data for this level
  AllData.Tp = squeeze(LargeScale.Grid.Tp(:,:,iLevel,:));
  AllData.MF = squeeze(LargeScale.Grid.MF(:,:,iLevel,:));  
  AllData.Kh = squeeze(LargeScale.Grid.Kh(:,:,iLevel,:));  
  AllData.Kz = squeeze(LargeScale.Grid.Kz(:,:,iLevel,:));  

  %OK. now, loop over the grid points and produce gappy time series
  %then use these to produce a seasonal pattern
  %see http://uk.mathworks.com/help/econ/seasonal-adjustment.html
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for iVar=1:1:4;
    
    switch iVar
      case 1; Data = AllData.Tp;
      case 2; Data = AllData.MF;
      case 3; Data = AllData.Kh;
      case 4; Data = AllData.Kz;
      otherwise; stop; 
    end

  
  NPoints = numel(LargeScale.Lon).*numel(LargeScale.Lat);
  Data = reshape(Data,NPoints,numel(LargeScale.Time));
  Cycle = NaN(NPoints,numel(Seas.Time));
  
  for iPoint=1:1:NPoints;

    %1. extract time series
    TimeSeries = squeeze(Data(iPoint,:))';
    if nansum(TimeSeries) == 0; continue; end %no data
    ChunkSize = 365./mean(diff(LargeScale.Time)); %number of time unit
    SmoothSize = ceil(Settings.SmoothSize./mean(diff(LargeScale.Time)));

    
% %     %2. smooth with a one-year moving average to detrend
%commented out as trend is known to be very small and this is the dominant
%bottleneck in the routine
% %     SmoothSeries = smooth(TimeSeries,SmoothSize);
% %     Detrended = TimeSeries-SmoothSeries;
  
    %3. create seasonal indices
    sidx = cell(ChunkSize,1);
    for iSm = 1:ChunkSize
      sidx{iSm,1} = iSm:ChunkSize:numel(TimeSeries);
    end; clear iSm
    
    %4. apply a stable seasonal filter
    sst = cellfun(@(x) nanmean(TimeSeries(x)),sidx);
    
%     %4. interpolate onto our new time scale
%%commented out as we're keeping the old one now
%     %(geographic interpolation will be done later)
%     OldTime = mean(diff(LargeScale.Time)):mean(diff(LargeScale.Time)):365;
%     OldTime = [OldTime-365,OldTime,OldTime+365];
%     sst     = [sst,sst,sst];
%     NewTime = Seas.Time;
%     valid = find(~isnan(sst));
%     sst = interp1(OldTime(valid),sst(valid),NewTime,'spline');
    

    %5. apply smoothing
    sst = smooth(sst,SmoothSize);


    %store
    Cycle(iPoint,:) = sst;

    %tidy
    clear TimeSeries SmoothSize SmoothSeries Detrended sidx sst OldTime NewTime valid
    
    
  end; clear iPoint NPoints
 
  %glue back into a map
  Cycle = reshape(Cycle,numel(LargeScale.Lon),numel(LargeScale.Lat),numel(Seas.Time));  
  
  %interpolate this map onto the output map, and store
  [xi,yi] = meshgrid(LargeScale.Lon,LargeScale.Lat);
  [xq,yq] = meshgrid(Seas.Lon,Seas.Lat);  
  for iTime=1:1:numel(Seas.Time);
    Out = interp2(xi,yi,squeeze(Cycle(:,:,iTime))',xq,yq)';
    switch iVar
      case 1; Seas.Tp(:,:,iLevel,iTime) = Out;
      case 2; Seas.MF(:,:,iLevel,iTime) = Out;
      case 3; Seas.Kh(:,:,iLevel,iTime) = Out;
      case 4; Seas.Kz(:,:,iLevel,iTime) = Out;
      otherwise; stop; 
    end
    clear Out
  end; clear iTime
  clear xi yi xq yq
  clear Cycle


  clear Data
  end; clear iVar
  clear AllData
end; clear iLevel

textprogressbar(' Done!')
save(Settings.OutFile,'Settings','Seas','-v7.3')
