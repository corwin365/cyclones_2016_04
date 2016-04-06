function generate_data(Instrument,LonRange,LatRange,TimeRange,OutPath)
% 

% 
% %temporary for testing
% clear all
% Instrument = 'HIRDLS';
% LonRange  = [-180,180];
% LatRange  = [-90,90];
% TimeRange = [datenum(2006,1,1),datenum(2006,1,2)];
% OutPath = './temp.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute lambda_z vertical climatology in desired region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%


Settings.NWavesToFind = 1; %this is equivalent to Alexander et al (JGR,2008)
Settings.Inst =  Instrument; %e.g. 'HIRDLS';

%coverage
Settings.LonRange = LonRange; clear LonRange
Settings.LatRange = LatRange; clear LatRange
Settings.Days  = TimeRange(1):1:TimeRange(2); clear TimeRange
Settings.NDays = numel(Settings.Days);

%minimum profile sep for MF. From Ern et al (JGR, 2011)
%we take the closest in each case, provided it's below this
Settings.MinSep    = 300.*1000; %m
%minimum time sep
Settings.MinDt     = 30./(24.*60); %fraction of day  


%dataset-specific settings
if strcmp(Settings.Inst,'SABER') == 1;
  Settings.DataDir = [LocalDataDir,'/SABER/reprocessed-v2/'];
  Settings.MinLambda = 4000;
  Settings.dZ        = 2000;
end
if strcmp(Settings.Inst,'COSMIC') == 1;
  Settings.DataDir = [LocalDataDir,'/AM_COSMIC/Cosmic_500m_10_60km/'];
  Settings.MinLambda = 2000;
  Settings.dZ        = 1000; 
end
if strcmp(Settings.Inst,'HIRDLS') == 1;
  Settings.DataDir = [LocalDataDir,'/HIRDLS/HIRDLS-L2/'];
  Settings.MinLambda = 2000;
  Settings.dZ        = 1000;
end
if strcmp(Settings.Inst,'MLS') == 1;
  Settings.DataDir     = [LocalDataDir,'/MLS/L2/'];
  Settings.MinLambda   = 7200; 
  Settings.dZ          = 2000; 
end
%don't do sondes with this routine - PW sections do not apply

%generally applicable settings
Settings.PWs.PWLat = 5;
Settings.PWs.NPWs  = 3;
Settings.MaxLambda   = 18000;
Settings.HeightRange = [18,80].*1000;  
Settings.OutFile  = OutPath; clear OutPath;
  
%heightscale
Settings.NHeightLevels = ceil(range(Settings.HeightRange)./Settings.dZ);
Settings.HeightScale   = Settings.HeightRange(1) + (0:1:Settings.NHeightLevels-1).*Settings.dZ;


%output prep
Results.Tp          = NaN(Settings.NHeightLevels,1); %the 1 will expand out
Results.MF          = Results.Tp;
Results.Kh          = Results.Tp;
Results.Kz          = Results.Tp;
Results.DayScale    = [NaN];
Results.ProfScale   = [NaN];
Results.HeightScale = Settings.HeightScale;
Results.Lat         = [NaN];
Results.Lon         = [NaN];
Results.Time        = [NaN];

%load the density climatology
Density = load('saber_density_new.mat'); Density = Density.Results;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loop over days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OldFile.Name = ' ';

for iDay=1:1:Settings.NDays
  disp(['===== Processing ',Instrument,' ',datestr(Settings.Days(iDay)),' =====']);

  %find profiles for this day, and format them
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if     strcmp(Settings.Inst,'SABER') == 1;
    [SatelliteData,OldFile] = extract_saber_data(Settings.Days(iDay),Settings.DataDir,OldFile);
  elseif strcmp(Settings.Inst,'COSMIC') == 1;
    [SatelliteData,OldFile] = extract_cosmic_data(Settings.Days(iDay),Settings.DataDir,OldFile);
  elseif strcmp(Settings.Inst,'HIRDLS') == 1;
    [SatelliteData] = extract_hirdls_data(Settings.Days(iDay),Settings.DataDir);
  elseif strcmp(Settings.Inst,'MLS') == 1;
    [SatelliteData] = extract_mls_data(Settings.Days(iDay),Settings.DataDir); 
  else
    disp('Error: extract routine not provided for this dataset');
    stop
  end
  
  if SatelliteData.Error ~=0 ; continue; end %problem extracting data

  
  FormattedData = format_data_pw(SatelliteData, ...
                                 Settings.LatRange, ...
                                 Settings.LonRange, ...
                                 Settings);
  if FormattedData.Error ~= 0; continue; end;   
  BackgroundTemp = FormattedData.BGTemperature;
  Perturbations  = FormattedData.Perturb;
  Pressure       = FormattedData.Pressure;
  clear SatelliteData 

  %fill any gaps, as they'll break the ST routine
  %don't do this to the BG temp, that way any bad temps will be removed at end
  Perturbations = inpaint_nans(Perturbations);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %work out the Brunt-Vaisala frequency
  %assume daily mean values for density and hence BVF
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  Dens = BackgroundTemp .* NaN;
  NProfiles = numel(FormattedData.Lat); 

  %work out density
  
  for i=1:1:NProfiles;
    %do every nth (arbitrary) level and interpolate between, for speed
    %loss of information should be very small and it's an approx. anyway
    for iLevel=1:3:Settings.NHeightLevels;
      if isfinite(BackgroundTemp(iLevel,i)) == 0; continue; end;
      Dens(iLevel,i) = cjw_airdensity(Pressure(iLevel),BackgroundTemp(iLevel,i));
    end; clear iLevel
    Dens(:,i) = inpaint_nans(Dens(:,i));
  end; clear i
  Dens = nanmean(Dens,2);
  
  %vertical derivative of density
  DensDeriv = Dens.*NaN;
  DensDeriv(2:Settings.NHeightLevels) = diff(Dens);
  DensDeriv = abs(DensDeriv)./Settings.dZ;
  
  %BVF
  BVF = 0.5.*sqrt((9.81./Dens).*DensDeriv);
  
  %duplicate to allow array multiplication later
  BVF = repmat(BVF,1,NProfiles);

  clear Dens DensDeriv  

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %find all profile pairs on the day
  %for most instruments, this will be the next profile
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %identify instrument - track limb sounders can be sped up significantly
  if     strcmp(Instrument,'HIRDLS') == 1; TrackInst = 1;
  elseif strcmp(Instrument,'SABER')  == 1; TrackInst = 1;
  elseif strcmp(Instrument,'MLS')    == 1; TrackInst = 1;
  else                                     TrackInst = 0;
  end
    
  textprogressbar('Finding Profile Separations  ')
  if TrackInst == 0;
    %full pair-finding routine, slow.
    
    %first, find the separation of each profile pair (in time and space)
    Distance = NaN(NProfiles,NProfiles);
    TimeSep  = Distance;
    
    for i=1:1:NProfiles;
      textprogressbar(i./NProfiles.*100)
      
      for j=1:1:i-1; %anything beyond this point produces duplicates
            
        %check time separation, as it's the most likely to lead to rejection
        dt = abs(FormattedData.Time(i)-FormattedData.Time(j));
        if TimeSep(i,j) > Settings.MinDt;  clear dt; continue; end %no point in going further
        TimeSep(i,j) = dt; clear dt;
        
        %very widely geographically separated profiles obviously won't produce pairs
        dLon = abs(FormattedData.Lon(i)-FormattedData.Lon(j));
        dLat = abs(FormattedData.Lat(i)-FormattedData.Lat(j));
        if abs(FormattedData.Lat(i)) < 70; %not valid near poles
          if dLon > 5; clear dLon dLat; continue; end
          if dLat > 5; clear dLon dLat; continue; end
        end
        clear dLon dLat
        
        %check remaining pairs for actual distance, to be sure
        Distance(i,j) =  distdim(distance(FormattedData.Lat(i),...
                                          FormattedData.Lon(i),...
                                          FormattedData.Lat(j),...
                                          FormattedData.Lon(j)), ...
                                          'deg','meters');
        
      end; clear i
    end; clear j
    
    
    %discard profile pairs too far apart in time or space
    BadTime  = find(TimeSep  > Settings.MinDt);
    BadSpace = find(Distance > Settings.MinSep);
    Bad = union(BadTime,BadSpace);
    Distance(Bad) = NaN;
    TimeSep( Bad) = NaN;
    clear BadTime BadSpace Bad
    
    
    %now, find the *nearest* profile to each other profile
    Nearest = NaN(3,NProfiles);
    for i=2:1:NProfiles;
      [a,b] = min(Distance(i,:));
      dt = TimeSep(i,b);
      if isnan(a); b = NaN; end
      Nearest(:,i) = [a,b,dt];
      clear a b dt;
    end; clear i
    clear Distance TimeSep
    
  else
    
    %much simpler! only consider adjacent profiles
    Nearest = NaN(3,NProfiles);
    for i=2:1:NProfiles;
      textprogressbar(i./NProfiles.*100)      
      dt = abs(FormattedData.Time(i)-FormattedData.Time(i-1));
      if dt > Settings.MinDt; clear dt; continue; end     
      dx = distdim(distance(FormattedData.Lat(i),...
                            FormattedData.Lon(i),...
                            FormattedData.Lat(i-1),...
                            FormattedData.Lon(i-1)), ...
                            'deg','meters');
      if dx > Settings.MinSep; clear dx; continue; end                          
   
      Nearest(:,i) = [dx,i-1,dt];
      clear dx dt
    end; clear i
  end;
  

  textprogressbar(' ')
  

  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %now loop over the profile pairs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  %results storage arrays
  Store.MF   = NaN.*FormattedData.Perturb;
  Store.Kh   = Store.MF;
  Store.Tp   = Store.MF;
  Store.Kz   = Store.MF;
  Store.Lat  = FormattedData.Lat.*NaN;
  Store.Lon  = FormattedData.Lon.*NaN;
  Store.Time = FormattedData.Time.*NaN;
  
  textprogressbar('ST Processing Profile Pairs  ')
  for iPair=1:1:NProfiles;
    textprogressbar(iPair./NProfiles.*100)    
    if ~isfinite(Nearest(2,iPair)); continue; end %no valid pair for this profile

    
    %extract the pair IDs
    ProfA = iPair;
    ProfB = Nearest(2,iPair);
    
    %geolocate their mean
    PairLats = [FormattedData.Lat(ProfA),FormattedData.Lat(ProfB)];
    PairLons = [FormattedData.Lon(ProfA),FormattedData.Lon(ProfB)];
    [Store.Lat(iPair),Store.Lon(iPair)] = meanm(PairLats,PairLons);
    Store.Time(iPair) = nanmean([FormattedData.Time(ProfA),FormattedData.Time(ProfB)]);    
    clear PairLats PairLons

    %zero pad profiles
    ProfileA = vertcat(zeros(20,1),Perturbations(:,ProfA),zeros(20,1));
    ProfileB = vertcat(zeros(20,1),Perturbations(:,ProfB),zeros(20,1));
    clear ProfA ProfB
    
    %S-Transform each profile
    [StOutA,~,StFrequency] = st(ProfileA,0,fix(length(ProfileA)/2),Settings.dZ,1);
    [StOutB,~,~]           = st(ProfileB,0,fix(length(ProfileB)/2),Settings.dZ,1);
    clear ProfileA ProfileB
    
    %unpad
    StOutA = StOutA(:,21:end-20);
    StOutB = StOutB(:,21:end-20);   
    
    %discard unwanted frequencies
    FreqMax = 1./Settings.MaxLambda;
    FreqMin = 1./Settings.MinLambda;
    StOutA = StOutA(StFrequency >= FreqMax & StFrequency <= FreqMin,:);
    StOutB = StOutB(StFrequency >= FreqMax & StFrequency <= FreqMin,:);    
    StFrequency = StFrequency(StFrequency >= FreqMax & StFrequency <= FreqMin);
    clear FreqMax FreqMin
    
    
    %produce cospectrum of the two profiles
    CoSpec = StOutA.*conj(StOutB);
    clear StOutA StOutB
        
    %find the peak amplitude at each level
    [~,PeakIdx] = nanmax(abs(CoSpec),[],1);
    
    %hence, find the amplitude, kz and delta-phase of the covarying signal
    Amp  = NaN(Settings.NHeightLevels,1);
    Dphi = Amp;
    Kz   = Amp;
    for i=1:1:Settings.NHeightLevels;
      Amp(i)  = sqrt(abs(CoSpec(PeakIdx(i),i)));
      Dphi(i) = angle(CoSpec(PeakIdx(i),i));
      Kz(i)   = StFrequency(PeakIdx(i));
    end; clear i 
    clear StFrequency CoSpec
    
    %convert Dphi to wavelength
    Dphi = Dphi ./ (2.*pi); %convert to fraction of circle
    Kh = abs(Dphi./Nearest(1,iPair));
    
    %work out relevant density profile
    [~,DayIdx] = min(abs(Settings.Days(iDay)-Density.DayScale));
    DensityProfile = interp1(Density.HeightScale,squeeze(Density.Density(DayIdx,:)),Settings.HeightScale)';
    clear DayIdx

    %and work out MF! (m^2 s^-2)
    MF = 0.5 .* DensityProfile                   ...
             .* (Kh./Kz)                         ...
             .* (9.81./BVF(:,iPair)).^2          ...
             .* (Amp./BackgroundTemp(:,iPair)).^2;
    %retain results
    Store.MF(:,iPair) = MF;
    Store.Kh(:,iPair) = Kh;
    Store.Tp(:,iPair) = Amp;
    Store.Kz(:,iPair) = Kz;
    
    
    clear Amp Dphi Kh Kz MF PeakIdx DensityProfile

  end; clear Nearest iPair
  textprogressbar(' ')
  clear BVF BackgroundTemp Perturbations FormattedData 


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %store results
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Results.MF = horzcat(Results.MF,Store.MF);
  Results.Kh = horzcat(Results.Kh,Store.Kh);
  Results.Tp = horzcat(Results.Tp,Store.Tp);
  Results.Kz = horzcat(Results.Kz,Store.Kz);  
  Results.Lat = vertcat(Results.Lat,Store.Lat');
  Results.Lon = vertcat(Results.Lon,Store.Lon');  
  Results.Time = vertcat(Results.Time,Store.Time');

  clear NProfiles Store 
  
end
disp('')
save(Settings.OutFile,'Results');
disp('--> Processed and saved')