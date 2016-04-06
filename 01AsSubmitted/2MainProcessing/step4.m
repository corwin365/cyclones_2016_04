function step4(Instrument)
% clear all
% Instrument = 'SABER'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do the assessment of each variable over each storm
%this is STEP FOUR of the analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Instrument = Instrument;
Settings.ClimaFile  = [LocalDataDir,'/corwin/storms/climatology/',Settings.Instrument,'.mat'];
Settings.StormFile  = 'storm_info.mat';
Settings.DataDir    = [LocalDataDir,'/corwin/limb_out5/'];

Settings.TimeWindow = 15; %days before and after storm centre to record data

Settings.OutFile    = [LocalDataDir,'/corwin/storms2/results/storms_',Settings.Instrument,'.mat'];

Settings.MaxDistance = 10; 
%degrees - make large, this is just intended as a processing time filter really
%the actual reduction to within-distance will be handled in later routines

%geographic etc binning settings come from earlier files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load climatologial data and prep arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the climatology that we'll be differencing from
Climatology = load(Settings.ClimaFile); Climatology = Climatology.Seas;

%and the storm data
Storms = load(Settings.StormFile); Storms = Storms.Storms;


% %and the grid spacings etc
Params = load(Settings.StormFile);
Settings.TimeRange = Params.Settings.TimeRange;
% clear Params

disp('--> Data loaded');

%main results
NLevs = numel(Climatology.z);
NVars = 4; %tp,mf,kh,kz
NPossiblePoints = 1e6;%trimmed or expanded as needed, this is just to minimise redeclaring arrays which is a major time bottleneck
Results.Results = NaN(NLevs,NVars,NPossiblePoints);
Results.StormIDs = NaN(6,NPossiblePoints); %storm, point, dlat, dLon, dt, distance
Results.z = Climatology.z;
clear NLevs NVars NPossiblePoints


disp('--> Storage arrays prepped');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%store vars for each point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%loop over months, as that's how the data are stored
[FirstYear,~,~] = datevec(min(Settings.TimeRange));
[LastYear, ~,~] = datevec(max(Settings.TimeRange));
%  FirstYear = 2006; LastYear = 2006; %testing
RowCount = 0; %for data storage. counts number of stormpoints 

textprogressbar('Processing data, ')
for iYear=FirstYear:1:LastYear
  for iMonth=1:1:12;
    Total = 12.*(LastYear-FirstYear+1);
    Sum = 12*(iYear-FirstYear)+iMonth;
    Percentage = Sum./Total.*100; 
    textprogressbar(Percentage)
    clear Total Sum Percentage
  
    %identify storage file
    FileName = [Settings.DataDir,'/',Settings.Instrument,'_',sprintf('%04d',iYear),'_',sprintf('%02d',iMonth),'.mat'];
    if exist(FileName) ~= 2; clear FileName; continue; end %no file
    Data = load(FileName); clear FileName; SatData = SatData.Results;

    

    %loop over storms and points within them
    for iStorm=1:1:numel(Storms.Time(1,:));
      for iPoint=1:1:numel(Storms.Time(:,1));

        %geolocation
        Lat  = Storms.Lat( iPoint,iStorm);
        Lon  = Storms.Lon( iPoint,iStorm);
        Time = Storms.Time(iPoint,iStorm);
        
        %if there's no data, skip
        if ~isnan(Time);
          
          %if no satellite time is near the storm time, skip
          dt = Data.Time-Time;
          mindt = dt(1); %this assumption *massively* reduces runtime
          if mindt > Settings.TimeWindow; continue; end
          
          %if there are no geographic points anywhere near the storm, skip
          dLat = abs(Lat-Data.Lat); if min(dLat) > Settings.MaxDistance; continue; end
          dLon = abs(Lon-Data.Lon); if min(dLon) > Settings.MaxDistance; continue; end
          
          %ok. find any satellite points close enough to the storm point
          %this may still be zero as the intersection of the two above is required
          NearStorm = find(abs(dLat) < Settings.MaxDistance ...
                         & abs(dLon) < Settings.MaxDistance ...
                         & abs(dt)   < Settings.TimeWindow);
          if numel(NearStorm) == 0; continue; end
          
          
          
          %found points! loop over them
          for iVal = 1:1:numel(NearStorm);
            %create a new row
            RowCount = RowCount+1;
            
            %find out how far away they are in time from our point, rounded onto our scale
            %           [deltat,idx] = min(abs(dt(NearStorm(iVal))-Results.TimeScale));
            deltat = dt(NearStorm(iVal));
            
            
            
            %geographic distance
            distance_to_storm = distdim(distance(Lat,...
                                                 Lon,...
                                                 Data.Lat(NearStorm(iVal)),...
                                                 Data.Lon(NearStorm(iVal))), ...
                                                 'deg','meters');
            
                                               
                                               
            
            Results.StormIDs(:,RowCount) = [iStorm,iPoint,dLat(NearStorm(iVal)),dLon(NearStorm(iVal)),deltat,distance_to_storm];
            
%              %deseasonalise
%              [~,latidx] = min(abs(Data.Lat(NearStorm(iVal))-Climatology.Lat));
%              [~,lonidx] = min(abs(Data.Lon(NearStorm(iVal))-Climatology.Lon));
%              [~,~,d] = datevec(Data.Time(NearStorm(iVal)));
%              [~,dayidx] = min(abs(d-Climatology.Time)); clear d
%              BG.Tp = squeeze(Climatology.Tp(lonidx,latidx,:,dayidx))';
%              BG.MF = squeeze(Climatology.MF(lonidx,latidx,:,dayidx))';
%  %             BG.Kh = squeeze(Climatology.Kh(lonidx,latidx,:,dayidx))';
%  %             BG.Kz = squeeze(Climatology.Kz(lonidx,latidx,:,dayidx))';
%              
%              clear lonidx latidx
            
            %store the results
            Results.Results(:,1,RowCount) = interp1(Data.HeightScale,Data.Tp(:,NearStorm(iVal)),Climatology.z);%-BG.Tp;
            Results.Results(:,2,RowCount) = interp1(Data.HeightScale,Data.MF(:,NearStorm(iVal)),Climatology.z);%-BG.MF;
            Results.Results(:,3,RowCount) = interp1(Data.HeightScale,Data.Kh(:,NearStorm(iVal)),Climatology.z);%-BG.Kh;
            Results.Results(:,4,RowCount) = interp1(Data.HeightScale,Data.Kz(:,NearStorm(iVal)),Climatology.z);%-BG.Kz;
            
          end; clear iVal
          clear NearStorm idx BG
          

        end; %if-time block ends here
      end;

    end; clear iStorm iPoint
    clear dt dLon dLat Lat Lon Time NearStorm mindt dt1 dt2


    clear Data
  end;clear iMonth
end; clear iYear FirstYear LastYear
textprogressbar(' Done!')


%truncate the data arrays from the possible maximum size to the actual used
%assume if we don't have T' we don't have anything, which I think should always be true
Used = squeeze(Results.Results(:,1,:));
Used = squeeze(nansum(Used,1));
Used = squeeze(nansum(Used,1));
Used = find(Used > 0);
Results.Results = Results.Results(:,:,Used);
Results.StormIDs = Results.StormIDs(:,Used);
clear Used

save(Settings.OutFile,'Results','-v7.3')
disp('--> Saved');