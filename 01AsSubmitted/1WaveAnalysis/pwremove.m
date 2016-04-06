function [BGTemp,Perturbations] = pwremove(SatelliteData,PWLat,NPWs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove PWs, by fitting sine waves to lat bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%copy of raw T, for background
SatelliteData.BGTemp = SatelliteData.Temp;

%loop over height
textprogressbar('Removing planetary waves     ')
for iLevel=1:1:numel(SatelliteData.Height(:,1));
  textprogressbar(iLevel./numel(SatelliteData.Height(:,1)).*100)
  
  PWCalc.Lat  = SatelliteData.Lat(       :);
  PWCalc.Lon  = SatelliteData.Lon(       :);
  PWCalc.Temp = SatelliteData.Temp(      iLevel,:)';
  
  
  %planetary wave removal
  for iLat=1:1:ceil(180./PWLat);
    
    %find boundaries of lat range
    LatBinStart = -90+(iLat-1).*PWLat;
    LatBinEnd   = -90+(iLat  ).*PWLat;
    
    %find profiles in bin
    InLatBin = find(PWCalc.Lat >= LatBinStart & PWCalc.Lat < LatBinEnd);
    if numel(InLatBin) == 0; clear InLatBin LatBinStart LatBinEnd; continue; end;
    
    TheLons  = PWCalc.Lon( InLatBin);
    TheTemps = PWCalc.Temp(InLatBin);
    
    %remove NaNs
    Valid = find(~isnan(TheLons + TheTemps));
    TheLons  = TheLons( Valid);
    TheTemps = TheTemps(Valid);

    
    if numel(Valid) == 0; continue; end
    
    %sort by lon
    [~,idx] = sort(TheLons,'ascend');
    TheLons  = TheLons(idx);
    TheTemps = TheTemps(idx);
    
    %fit planetary wave modes
    for iPW=1:1:NPWs;
      %identify frequency of signal we're looking for
      WaveFreq = 1./(360./iPW); %per degree
      
      %compute the wave using sinefit (IEEE-1057 standard fitting)
      %flags: verbose,plot,iterate to improve
      
      [~,Wave] = sinefit(TheTemps,TheLons,WaveFreq,0,0,0);
      TheTemps = TheTemps - Wave; %remove PW
      
      clear Wave WaveFreq
    end; clear iPW
    
    %put the temperatures back where they came from
    SatelliteData.Temp(iLevel,(InLatBin(Valid(idx)))) = TheTemps;
    
    clear LatBinStart LatBinEnd TheLons TheLats TheTemps idx Valid
  end; clear iLat
  
  clear PWCalc
end; clear iLevel
textprogressbar(' ')

%remove any remaining outlier T (e.g. only one point in the latbin will not fit successfully)
SatelliteData.Temp(abs(SatelliteData.Temp) > 100) = NaN;

%compute the background
BGTemp = SatelliteData.BGTemp - SatelliteData.Temp;
Perturbations =  SatelliteData.Temp;

end