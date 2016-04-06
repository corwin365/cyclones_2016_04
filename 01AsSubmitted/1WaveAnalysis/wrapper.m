


LonRange  = [-180,180];
LatRange  = [-90,90];

for iYear=2002:1:2013;
  for iMonth = 1:1:12;
    
    TimeRange = [datenum(iYear,iMonth,1),datenum(iYear,iMonth+1,0)];
    
    OutPath = [LocalDataDir,'/corwin/limb_out5/HIRDLS_',sprintf('%04d',iYear),'_',sprintf('%02d',iMonth),'.mat'];
    if exist(OutPath) == 0; generate_data('HIRDLS',LonRange,LatRange,TimeRange,OutPath); end

    OutPath = [LocalDataDir,'/corwin/limb_out5/SABER_',sprintf('%04d',iYear),'_',sprintf('%02d',iMonth),'.mat'];
    if exist(OutPath) == 0; generate_data('SABER',LonRange,LatRange,TimeRange,OutPath); end
    
    OutPath = [LocalDataDir,'/corwin/limb_out5/MLS_',sprintf('%04d',iYear),'_',sprintf('%02d',iMonth),'.mat'];
    if exist(OutPath) == 0; generate_data('MLS',LonRange,LatRange,TimeRange,OutPath); end
    
    %didn't use this - too little data to be meaningful as only a few
    %hundred profiles globally are pairable
%     if datenum(iYear,iMonth,30) <= datenum(2007,5,1);
%       OutPath = [LocalDataDir,'/corwin/limb_out2/COSMIC_',sprintf('%04d',iYear),'_',sprintf('%02d',iMonth),'.mat'];
%      if exist(OutPath) == 0; generate_data('COSMIC',LonRange,LatRange,TimeRange,OutPath); end
%     end
    
    
  end; 
end;
