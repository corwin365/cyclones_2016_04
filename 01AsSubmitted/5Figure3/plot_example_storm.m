clear all

clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.02], [0.1 0.1], [0.1 0.1]);



%%ronseal....


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data - T'
load('../../02PhaseTwo/example_storm_1.mat')

Storm   = EXAMPLE.Storm;
ResultsTp = EXAMPLE.Results;
clear EXAMPLE

%load data - MF
load('../../02PhaseTwo/example_storm_2.mat')
ResultsMF = EXAMPLE.Results;
clear EXAMPLE

% load and interpolate topo data

%prepare topography
%%%%%%%%%%%%%%%%%%%%
load([LocalDataDir,'/Topography/easy_tenth_degree_topography/easy_topo.mat']);
topo.elev(topo.elev <= -1) =NaN;%-300; %nicer blue for the sea

% %interpolate to a less intensive scale
% lons = -180:2.5:180;
% lats = -90:2.5:90;
% [lons,lats] = meshgrid(lons,lats);
% elev = interp2(topo.lons,topo.lats,topo.elev,lons,lats);
% topo.elev = elev;
% topo.lons = lons;
% topo.lats = lats;
% clear elev lons lats



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first panel: map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,4)

m_proj('satellite','lat',32,'lon',140,'radius',0.25);
CT = cbrewer('seq','Greys',36); CT = CT(6:36,:); colormap(CT);
m_pcolor(topo.lons,topo.lats,topo.elev); shading flat
caxis([0 3000]);freezeColors
hold on
m_coast('color','k','linewidth',1 );
% m_grid
m_grid('color','k','xticklabels',[],'yticklabels',[]);

%plot storm track
Colour1 = [153,0,0]./255;
%make it an alphaed line around it, for visual effect
Lon = Storm.Lon; 
Lat = Storm.Lat;
Good = find(~isnan(Lon+Lat));
Lon = Lon(Good,:);
Lat = Lat(Good,:);
clear Good

%overinterpolate and smooth, for pretty
NewScale = 1:0.01:numel(Lat);
Lon = interp1(1:1:numel(Lat),Lon,NewScale,'spline')'; Lon = smooth(Lon,400);
Lat = interp1(1:1:numel(Lat),Lat,NewScale,'spline')'; Lat = smooth(Lat,400);
clear NewScale

%plot as a patch
Lon(:,2) = Lon+1; Lon(:,1) = Lon(:,2)-2;
Lon = [Lon(:,1);flipud(Lon(:,2));Lon(1,1)];
Lat = [Lat(:,1);flipud(Lat(:,1));Lat(1,1)];
m_patch(Lon,Lat,Colour1,'facealpha',0.6,'edgealpha',0.7)
clear Lon lat

title('Typhoon Nabi, 2005','FontSize',14)

%%%%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second plot: wind time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,[1,2,3])
box on
set(gca,'FontSize',12)

%get data
Time = Storm.Time;
Wind = Storm.Wind./1.94; %kts -> m/s


%shade days
for iDay=floor(min(Time)):2:ceil(max(Time));
  h = fill([iDay-1,iDay,iDay,iDay-1,iDay-1],[0 0 100 100 0],[1,1,1].*0.9,'edgecolor','none');
  hold on
end



%overinterpolate and spline
Time2 = min(Time):0.01:max(Time);
Wind2 = interp1(Time,Wind,Time2,'spline');



plot(Time2,Wind2,'k-','linewidth',2)
datetick('x','dd/mm'); set(gca,'xtick',min(Time)-1:1:max(Time)+1)
clear Time2 Wind2

ylabel('Wind Speed [m/s]')
%create axes
axis([datenum(2005,08,30),datenum(2005,09,08),0,100./1.94])


%position figure label and map figure label
text(datenum(2005,08,30,04,00,00),15,'(a)','FontSize',25)
text(datenum(2005,09,10,12,00,00),10./1.94,'(b)','FontSize',25)


%add lat/lon information
for iDay=datenum(2005,08,30,12,00,00):1:datenum(2005,09,07,12,00,00);
  [~,idx] = min(abs(iDay-Storm.Time+1));
  Lat = Storm.Lat(idx); Lat = [num2str(Lat),'N'];
  text(iDay,113./1.94,Lat,'clipping','off','horizontalalignment','center','fontsize',11)
  Lon = Storm.Lon(idx); Lon = [num2str(Lon),'E'];  
  text(iDay,105./1.94,Lon,'clipping','off','horizontalalignment','center','fontsize',11)
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%third plot: observed GW T'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(3,4,[1,2,3]+4)

hold on; box on
set(gca,'FontSize',12)
% datetick('x',19)

%shade days
for iDay=floor(min(Time)):2:ceil(max(Time));
  fill([iDay,iDay+1,iDay+1,iDay,iDay],[-15 -15 15 15 -15],[1,1,1].*0.95,'edgecolor','none')
end; clear iDay

%zero axis
plot([min(Time)-1 max(Time)+1],[0 0],'-','linewidth',2,'color',[1,1,1].*0.4)
ylabel('Days since typhoon centre')

%daily axes
for iDay=-14:1:14;
  if iDay == 0; continue; end
  plot([min(Time)-1 max(Time)+1],[iDay iDay],'-','linewidth',1,'color',[1,1,1].*0.8)
end ; clear iDay

Colours = flipud(cbrewer('div','RdYlBu',32)); %cbrewer('seq','YlOrBr',32);
CAx     = 0+(1:1:32)./32.*4;
Symbols = 'ox^s';


%plot individual days
for iPoint=1:1:numel(ResultsTp(:,1));
  
  t  = ResultsTp(iPoint,1)+2;
  dt = ResultsTp(iPoint,2);
  v  = ResultsTp(iPoint,3);
  ds = ResultsTp(iPoint,4);
  if isnan(t+dt+v+ds); continue; end
  
  [~,colour] = min(abs(CAx-v));
  colour = Colours(colour,:);
  
  plot(t,dt,Symbols(ds),'color',colour,'markerfacecolor',colour,'markersize',8)
  
end; clear iPoint
datetick('x','dd/mm'); set(gca,'xtick',min(Time)-1:1:max(Time)+1)
axis([datenum(2005,08,31),datenum(2005,09,09),-15,15])


colormap(flipud(cbrewer('div','RdYlBu',32)));
c = colorbar('manual','position',[0.04 0.4 0.02 0.2]);
caxis([min(CAx) max(CAx)])
c.Label.String = 'T'' [\timesdataset mean]';
c.Label.FontSize = 12;

clear t dt v ds colour CAx Colours Symbols

text(datenum(2005,08,31,1,00,00),11,'(c)','FontSize',25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fourth plot: observed GW MF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,[1,2,3]+8)
hold on; box on
set(gca,'FontSize',12)


%shade days
for iDay=floor(min(Time)):2:ceil(max(Time));
  fill([iDay,iDay+1,iDay+1,iDay,iDay],[-15 -15 15 15 -15],[1,1,1].*0.95,'edgecolor','none')
end; clear iDay

%zero axis
plot([min(Time)-1 max(Time)+1],[0 0],'-','linewidth',2,'color',[1,1,1].*0.4)
ylabel('Days since typhoon centre')

%daily axes
for iDay=-14:1:14;
  if iDay == 0; continue; end
  plot([min(Time)-1 max(Time)+1],[iDay iDay],'-','linewidth',1,'color',[1,1,1].*0.8)
end ; clear iDay

Colours = flipud(cbrewer('div','RdYlBu',32)); %cbrewer('seq','YlOrBr',32);
CAx     = 0+(1:1:32)./32.*10;
Symbols = 'ox^s';

%plot individual days
for iPoint=1:1:numel(ResultsMF(:,1));
  
  t  = ResultsMF(iPoint,1)+2;
  dt = ResultsMF(iPoint,2);
  v  = ResultsMF(iPoint,3);
  ds = ResultsMF(iPoint,4);
  if isnan(t+dt+v+ds); continue; end
  
  [~,colour] = min(abs(CAx-v));
  colour = Colours(colour,:);
  
  plot(t,dt,Symbols(ds),'color',colour,'markerfacecolor',colour,'markersize',8)
  
end; clear iPoint

datetick('x','dd/mm'); set(gca,'xtick',min(Time)-1:1:max(Time)+1)
axis([datenum(2005,08,31),datenum(2005,09,09),-15,15])


colormap(flipud(cbrewer('div','RdYlBu',32)));
c = colorbar('manual','position',[0.04 0.12 0.02 0.2]);
caxis([min(CAx) max(CAx)])
c.Label.String = 'MF [\timesdataset mean]';
c.Label.FontSize = 12;

clear t dt v ds colour CAx Colours Symbols

text(datenum(2005,08,31,1,00,00),11,'(c)','FontSize',25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%let's handle the key here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(datenum(2005,09,02,-3,00,00),-23,'ko','clipping','off','markerfacecolor','k')
text(datenum(2005,09,02,00,00,00),-23,'HIRDLS','FontSize',14,'clipping','off')

plot(datenum(2005,09,03,01,00,00),-23,'k^','clipping','off','markerfacecolor','k')
text(datenum(2005,09,03,04,00,00),-23,'MLS','FontSize',14,'clipping','off')

plot(datenum(2005,09,04,00,00,00),-23,'ks','clipping','off','markerfacecolor','k')
text(datenum(2005,09,04,03,00,00),-23,'SABER','FontSize',14,'clipping','off')

%no COSMIC data present in figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fifth plot: histogram of T'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,8)
axis([0 15 -15 15])
hold on; box on; set (gca,'fontsize',12)
set(gca,'ytick',{})
plot([0 25 ],[0 0],'-','linewidth',2,'color',[1,1,1].*0.4)
for iDay=-15:1:15; if iDay == 0; continue; end; plot([0 25],[1,1].*iDay,'-','linewidth',1,'color',[1,1,1].*0.8); end

dt = ResultsTp(:,2);
MF = ResultsTp(:,3);

plot(MF,dt,'ko','markerfacecolor','k','markersize',4)
ylim([-15,15])
xlabel('T'' [\timesdataset mean]')
text(15,1-10,'(d)','FontSize',25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sixth plot: histogram of MF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,4,12)
axis([0 24 -15 15])
hold on; box on; set (gca,'fontsize',12)
set(gca,'ytick',{})
plot([0 25 ],[0 0],'-','linewidth',2,'color',[1,1,1].*0.4)
for iDay=-15:1:15; if iDay == 0; continue; end; plot([0 25],[1,1].*iDay,'-','linewidth',1,'color',[1,1,1].*0.8); end

dt = ResultsMF(:,2);
MF = ResultsMF(:,3);



plot(MF,dt,'ko','markerfacecolor','k','markersize',4)
ylim([-15,15])
xlabel('MF [\timesdataset mean]')

text(20,1-10,'(f)','FontSize',25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%done!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

export_fig('example_storm','-pdf','-native')