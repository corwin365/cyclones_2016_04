clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%georeferenced image of hurricane Nabi, to show analysis on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load image
%stored upside down because I can't work out how to flip in code...
[P,~]=imread('Typhoon_Nabi_06_sep_2005_0205Z_flipped.jpg');

%georeferencing coordinates eyeballed from google maps
%close enough at the scales we're worked on, and checked by eye against overplotted land
m_proj('mercator', ...
       'long',[120 141],...
       'lat',[ 21   40]);

clf;
set(gcf,'color','w')
image(P);
hold on
clear P 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%overplot the path of the storm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../2MainProcessing/storm_info.mat')
Storm = 340; %ID in the dataset as I have it.

%plot storm track
Colour1 = [255,128,0]./255;

%make it an alphaed line around it, for visual effect
Lon = Storms.Lon(:,Storm); 
Lat = Storms.Lat(:,Storm);
Good = find(~isnan(Lon+Lat));
Lon = Lon(Good,:);
Lat = Lat(Good,:);
clear Good

%plot daily points
Store = 0;
a = []; b = [];
for i=1:1:142;
  dt = Storms.Time(i,Storm)-Store;
  if isnan(dt); continue; end
  if dt < 0.25; continue; end %quarter-day
  Store = Storms.Time(i,Storm);
  m_plot(Storms.Lon(i+Shift,Storm),Storms.Lat(i+Shift,Storm),'o','color',Colour1,'linewidth',1,'markerfacecolor',Colour1,'markersize',8)
  a(end+1) = Storms.Lon(i+Shift,Storm);
  b(end+1) = Storms.Lat(i+Shift,Storm);
end;
 m_plot(a,b,'-','color',Colour1,'linewidth',2)
clear dt Store i Shift
  

%highlight nearest time to image
ImageTime = datenum(2005,09,06,02,05,00);
[~,idx] = min(abs(Storms.Time(:,Storm) - ImageTime));
Lat = Storms.Lat(idx,Storm);
Lon = Storms.Lon(idx,Storm);
Storms.Wind(idx,Storm)
datestr(Storms.Time(idx,Storm))

Colour2 = [102,0,204]./255;
m_plot(Lon,Lat,'o','color',Colour2,'linewidth',1,'markerfacecolor',Colour2,'markersize',15)

%plot our averaging area
m_range_ring(Lon,Lat,300,'color',Colour2,'linewi',3);

clear Colour1 Colour2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%overplot the satellites on a given day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%HIRDLS
HIRDLSFile = [LocalDataDir,'/HIRDLS/HIRDLS-L2/2005/09/HIRDLS-Aura_L2_v07-00-20-c01_2005d249.he5'];
Hirdls.Lats = get_HIRDLS(HIRDLSFile,'Latitude');
Hirdls.Lons = get_HIRDLS(HIRDLSFile,'Longitude');
m_plot(Hirdls.Lons,Hirdls.Lats,'b.-','linewidth',2,'markersize',30,'clipping','on')

%MLS
MLSFile = [LocalDataDir,'/MLS/L2/2005/MLS-Aura_L2GP-Temperature_v03-30-c01_2005d249.he5'];
Mls = get_MLS(MLSFile,'Temperature');
m_plot(Mls.Longitude,Mls.Latitude,'r.-','linewidth',2,'markersize',30,'clipping','on')

%SABER
SABERFile = [LocalDataDir,'/SABER/reprocessed-v2/September2005_v2.0.mat'];
SABER = load(SABERFile); SABER = cleandata_saber(SABER.SaberData);
SABER.Lat = SABER.Lat(14,SABER.MatlabTime(14,:) >= datenum(2005,9,6) & SABER.MatlabTime(14,:) < datenum(2005,9,7));
SABER.Lon = SABER.Lon(14,SABER.MatlabTime(14,:) >= datenum(2005,9,6) & SABER.MatlabTime(14,:) < datenum(2005,9,7));
%odd jumps due to how m_map works, fix
Deriv = [abs(diff(SABER.Lon)),9999]; SABER.Lon(Deriv > 5) = NaN;
%first line: all points
m_plot(SABER.Lon,SABER.Lat,'k.:','linewidth',2,'markersize',30,'clipping','on')
%second line: only close pairs
for i=1:1:numel(SABER.Lon)-1;
  if mod(i,2) == 1; continue; end
  m_plot([SABER.Lon(i),SABER.Lon(i+1)], ...
         [SABER.Lat(i),SABER.Lat(i+1)],'k.-','linewidth',2,'markersize',30,'clipping','on')
end

clear MLS HIRDLS SABER HIRDLSFile MLSFile SABERFile Deriv

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%finish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'ydir','normal');
m_grid('linest',':','FontSize',14)
title('Typhoon Nabi, 2005/09/06 02:05Z','FontSize',16)

export_fig('storm_image','-png','-m1.5')