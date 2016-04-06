function step2(Instrument)
% clear all
% Instrument = 'SABER'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compute the large-scale MF distribution of each instruments, so we
%can compute the deseasonalisng terms
%this is STEP TWO of the analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Settings.Instrument = Instrument;
Settings.MFDir      = [LocalDataDir,'/corwin/limb_out_multi/'];
Settings.OutFile    = [LocalDataDir,'/corwin/storms2/background/',Settings.Instrument,'.mat'];

Settings.TimeRange  =[datenum(2002,1,1),datenum(2013,12,31)];
Settings.LatRange   = [ -50, 50];
Settings.LonRange   = [-180,180];
Settings.z          = [20:2:60].*1000;

Settings.LatStep  = 10;
Settings.LonStep  = 20;
Settings.TimeStep = 5; %days


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create output arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BGGrid.Time   = Settings.TimeRange(1):Settings.TimeStep:Settings.TimeRange(2);
BGGrid.Lon    = Settings.LonRange(1):Settings.LonStep:Settings.LonRange(2);
BGGrid.Lat    = Settings.LatRange(1):Settings.LatStep:Settings.LatRange(2);
BGGrid.z      = Settings.z;
BGGrid.Grid.Tp = NaN(numel(BGGrid.Lon),numel(BGGrid.Lat),numel(BGGrid.z),numel(BGGrid.Time));
BGGrid.Grid.MF = BGGrid.Grid.Tp;
BGGrid.Grid.Kh = BGGrid.Grid.Tp;
BGGrid.Grid.Kz = BGGrid.Grid.Tp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%average data onto grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%simple things outside loop to avoid repetition
[xi,yi] = meshgrid(BGGrid.Lon,BGGrid.Lat);

OldFile = '';
textprogressbar('Computing Background  ')
for iTime=1:1:numel(BGGrid.Time)-1;
  textprogressbar(iTime./(numel(BGGrid.Time)-1).*100);
  try
    %what days do we want?
    DaysToUse = BGGrid.Time(iTime):1:BGGrid.Time(iTime+1);
    
    %temporary storage
    Tp.Data  = BGGrid.Grid.Tp(:,:,:,iTime); Tp.Data(:) = 0; Tp.Count = Tp.Data;
    MF.Data = Tp.Data; MF.Count = Tp.Count;
    Kh.Data = Tp.Data; Kh.Count = Tp.Count;
    Kz.Data = Tp.Data; Kz.Count = Tp.Count;
    
    %fill temporary storage
    for iDay=1:1:numel(DaysToUse);
      
      %find which file to load
      [y,m,~] = datevec(DaysToUse(iDay));
      FileName = [Settings.MFDir,'/',Settings.Instrument,'_',sprintf('%04d',y),'_',sprintf('%02d',m),'.mat'];
      clear y m
      if exist(FileName)~= 2; clear FileName; continue; end %file does not exist
      
      %load the file
      if strcmp(FileName,OldFile) ~=1;
        %load new file
        SatData = load(FileName); SatData = SatData.Results;
        OldFile = FileName;
       
        
      end;
      clear FileName

      %extract just this day and grid onto our output height axis
      OnDay = find(floor(SatData.Time) == DaysToUse(iDay));

      if numel(OnDay) == 0; continue; end %no data this day
      
      %interpolate the data to our background grid
      %we can't use interp as there are NaNs everywhere (T' < noise)
      %so use a nearest-neighbour basis
      DailyData.Lat = SatData.Lat(OnDay); DailyData.Lon = SatData.Lon(OnDay);
      for iVar=1:1:4;
        switch iVar
          case 1; Var = SatData.Tp(:,OnDay); 
          case 2; Var = SatData.MF(:,OnDay); 
          case 3; Var = SatData.Kz(:,OnDay); 
          case 4; Var = SatData.Kh(:,OnDay); 
        end
        VarOut = NaN(numel(BGGrid.z),numel(OnDay));
        

        for iLev=1:1:numel(BGGrid.z);
          [~,idx] = min(abs(BGGrid.z(iLev)-SatData.HeightScale));
          VarOut(iLev,:) = Var(idx,:);
        end
        
        switch iVar
          case 1; DailyData.Tp = VarOut; 
          case 2; DailyData.MF = VarOut; 
          case 3; DailyData.Kz = VarOut; 
          case 4; DailyData.Kh = VarOut; 
        end
        
        
      end; clear iVar


      %append to our grids
      for iLevel=1:1:numel(Settings.z);

        
        %geolocation
        x = DailyData.Lon(:); y = DailyData.Lat(:);
        
        %each variable
        for iVar=1:1:4;
          switch iVar;
            case 1; z = DailyData.Tp(iLevel,:);
            case 2; z = DailyData.MF(iLevel,:);
            case 3; z = DailyData.Kh(iLevel,:);
            case 4; z = DailyData.Kz(iLevel,:);
            otherwise stop;
          end;
          
          %bin
          valid = find(~isnan(x+y+z'));
          a = bin2mat(x(valid),y(valid),z(valid),xi,yi,'@sum');
          c = z; c(valid) = 1;
          b = bin2mat(x(valid),y(valid),c(valid),xi,yi,'@sum');
          clear valid
          
          %store by variable
          switch iVar
            case 1; Tp.Data(:,:,iLevel) = squeeze(Tp.Data(:,:,iLevel))+a';
            case 2; MF.Data(:,:,iLevel) = squeeze(MF.Data(:,:,iLevel))+a';
            case 3; Kh.Data(:,:,iLevel) = squeeze(Kh.Data(:,:,iLevel))+a';
            case 4; Kz.Data(:,:,iLevel) = squeeze(Kz.Data(:,:,iLevel))+a';
            otherwise stop;
          end;
          switch iVar
            case 1; Tp.Count(:,:,iLevel) = squeeze(Tp.Count(:,:,iLevel))+b';
            case 2; MF.Count(:,:,iLevel) = squeeze(MF.Count(:,:,iLevel))+b';
            case 3; Kh.Count(:,:,iLevel) = squeeze(Kh.Count(:,:,iLevel))+b';
            case 4; Kz.Count(:,:,iLevel) = squeeze(Kz.Count(:,:,iLevel))+b';
            otherwise stop;
          end;
          clear a b c z
        end; clear iVar
        clear x y
      end; clear iLevel;
      
      clear DailyData OnDay
    end; clear iDay
    clear DaysToUse
    
    %store the data
    BGGrid.Grid.Tp(:,:,:,iTime) = Tp.Data./Tp.Count;
    BGGrid.Grid.MF(:,:,:,iTime) = MF.Data./MF.Count;
    BGGrid.Grid.Kz(:,:,:,iTime) = Kz.Data./Kz.Count;
    BGGrid.Grid.Kh(:,:,:,iTime) = Kh.Data./Kh.Count;
    
    clear Tp MF Kz Kh
  catch; end
  
end; clear iTime OldFile
textprogressbar(' Done!')
clear xi yi SatData


%save output
save(Settings.OutFile, 'BGGrid','Settings','-v7.3');