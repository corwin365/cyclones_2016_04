clear all

clf; set(gcf,'color','w')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%map/time series of IBTRacS storm climatology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


StormData = load('../2MainProcessing/storm_info.mat'); 
StormData = StormData.Storms;

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first panel: binned map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%prepare topography
%%%%%%%%%%%%%%%%%%%%
load([LocalDataDir,'/Topography/easy_tenth_degree_topography/easy_topo.mat']);
topo.elev(topo.elev <= -1) =NaN;

%interpolate to a less intensive scale
lons = -180:0.5:180;
lats = -90:0.5:90;
[lons,lats] = meshgrid(lons,lats);
elev = interp2(topo.lons,topo.lats,topo.elev,lons,lats);
topo.elev = elev;
topo.lons = lons;
topo.lats = lats;
clear elev lons lats

% %shift onto a 0-360 range
% topo.lons = wrapTo360(topo.lons);
% clear a b c
disp('go!')



subplot(2,1,1)

Lat = -90:1:90; Lon = -180:1:180;
[xi,yi] = meshgrid(Lon,Lat);
x = StormData.Lon(:);
y = StormData.Lat(:);
z = ones(numel(StormData.Lon(:)),1);
valid = find(~isnan(x+y+z));
zz = bin2mat(x(valid),y(valid),z(valid),xi,yi,'@sum');
clear xi yi x y z valid


m_proj('Robinson','long',[-180 180],'lat',[-1,1].*60)
CT = cbrewer('seq','Greys',36);
CT = CT(6:36,:);
colormap(CT);


m_pcolor(topo.lons,topo.lats,topo.elev); shading flat
caxis([-200 6000])
freezeColors
hold on




colormap(flipud(cbrewer('seq','Blues',32)));
m_pcolor(Lon,Lat,zz); shading flat; 
c = colorbar; 
c.Label.String = 'Number of storm points recorded'; 
c.Label.FontSize = 12;

caxis([0 nanmax(zz(:))])

m_coast('color','k', 'linewi',1   );
m_grid('linest','--')

clear c CT Lat Lon topo zz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second panel: temporal histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)
TimeScale = datenum(2002,(1:1:(12.*12)),1);

%create axes
axis([min(TimeScale) max(TimeScale) 0 800]); hold on; box on;

%background fill
for iYear=2002:2:2013;
  x = [datenum(iYear,1,1),datenum(iYear+1,1,0),datenum(iYear+1,1,0),datenum(iYear,1,1),datenum(iYear,1,1)];
  y = [0 0 800 800 0];
  
  fill(x,y,[1,1,1].*0.9,'edgecolor','none')
  
end

a = hist(StormData.Time(:),TimeScale);

bar(TimeScale,a,'barwidth',1,'facecolor',[204,0,0]./255,'edgecolor',[153,0,0]./255)
datetick;
ylabel('Number of storm points recorded')
% 

text(datenum(2001,1,1),1900,'(a)','FontSize',18,'clipping','off')
text(datenum(2001,1,1),900,'(b)','FontSize',18,'clipping','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawnow
export_fig('storm_map','-pdf','-native')

