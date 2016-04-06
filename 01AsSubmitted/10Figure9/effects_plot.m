clear all
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.05], [0.05 0.05], [0.05 0.05]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the aggregate effect of cyclones in each basin
%take outermost bounds of each instrument pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load('effects_data.mat');

%days smoothing to apply
SmoothDays = 31;

%colours for instruments, in order
InstColours = 'rkb';

%colour maps for regions
Colours = {'Reds','','Blues','Greens','Oranges','Purples','YlOrBr'};

%plotting order, to put time series near regions in output figure
PlotPos     = [2,3,1,10,11,12];
PlotLetters = 'bcadef'; 
BarY        = [2,3,1,4,5,6];

%basins
%all +1 in our arrays
%0 = NA - North Atlantic 
%1 = SA - South Atlantic 
%2 = WP - West Pacific 
%3 = EP - East Pacific 
%4 = SP - South Pacific 
%5 = NI - North Indian 
%6 = SI - South Indian 
Basins = {'North Atlantic','South Atlantic','West Pacific', ...
          'East Pacific','South Pacific','North Indian','South Indian'};
%omit south atlantic from plotting below - almost no data
       



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare topography
%%%%%%%%%%%%%%%%%%%%
load([LocalDataDir,'/Topography/easy_tenth_degree_topography/easy_topo.mat']);
topo.elev(topo.elev <= -1) =NaN;

%interpolate to a less intensive scale
lons = -180:.5:180;
lats = -90:.5:90;
[lons,lats] = meshgrid(lons,lats);
elev = interp2(topo.lons,topo.lats,topo.elev,lons,lats);
topo.elev = elev;
topo.lons = lons;
topo.lats = lats;
clear elev lons lats

%plot map and data
%%%%%%%%%%%%%%%%%%%%

subplot(4,4,[5,6,7,9,10,11])

%produce the projection
m_proj('Robinson','long',[-180 180],'lat',[-1,1].*60)

% %first, plot topography in greyscale
CT = cbrewer('seq','Greys',36); CT = CT(6:36,:); colormap(CT);
m_pcolor(topo.lons,topo.lats,topo.elev); shading flat
caxis([-200 6000])
freezeColors; hold on

%now, loop over the basins and plot the median MF per gridbox
%use HIRDLS, as it's the median of the three

iInst = 2;
for iBasin=1:1:7;
  if iBasin == 2; continue; end %no point in South Atlantic
  
  %load regional colour map
  CT = cbrewer('seq',Colours{iBasin},128); CT = CT(16:128,:); 
  colormap(CT);
  
  %plot data
  Data = squeeze(Maps.Data(iInst,iBasin,:,:));
  m_pcolor(Maps.Lon,Maps.Lat,Data); shading flat;

  %store the ends of the bars, for the key
  MinStore(iBasin) = min(Data(:));%1.2;
  MaxStore(iBasin) = max(Data(:));
%   %set all to same range, more comparable
%   MinStore(iBasin) = 0;%min(Data(:));%1.2;
%   MaxStore(iBasin) = 0.2;%max(Data(:));
  
  
  %freeze colours
  caxis([MinStore(iBasin) MaxStore(iBasin)])  
  freezeColors;
  drawnow
end
  
clear Data iBasin iInst CT 



%now, finish up the plot!
m_coast('color','k', 'linewi',1   );
m_grid('linest','--')
hold on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop over basins

jBasin = 0;
for iBasin = 1:1:7;
  if iBasin == 2; continue; end
  jBasin = jBasin +1;
 
  %produce a panel, with an appropriately-shaded background
  subplot(4,3,PlotPos(jBasin));
  hold on; box on
  
  %produce axes with shading
  axis([0 365 0 7.5])
  Colour = cbrewer('seq',Colours{iBasin},36); Colour = Colour(2,:);
  patch([0 365 365 0 0],[0 0 10 10 0],Colour,'edgecolor','none')
    
  %produce a grid
  MonthStarts = zeros(12,1);
  for iMonth = 2:1:12; MonthStarts(iMonth) = MonthStarts(iMonth-1)+cjw_nmonthdays(iMonth,1990); end
  for iMonth=1:1:12; plot([1,1].*MonthStarts(iMonth),[-999,999],':','color',[1,1,1].*0.3); end
  for i=1:1:50; plot([0 365],[1,1].*i,':','color',[1,1,1].*0.3); end
  
  
  %plot medians for each instrument
  Max = 0;
  for iInst=1:1:3;
    Data = squeeze(TimeSeries.Data(iInst,iBasin,:));
    Data = smooth([Data,Data,Data],SmoothDays./mean(diff(TimeSeries.Time)));
    Data = Data(numel(TimeSeries.Time)+1:numel(TimeSeries.Time).*2);
    
    plot(TimeSeries.Time,Data,'color',InstColours(iInst),'linewidth',1);
    
    if max(Data(:)) > Max; Max = max(Data(:)); end
    
    %work out and plot the mean across the basin, for discussion
    TheMean = mean(Data(:));
    plot([0,365],[1,1]*TheMean,'--','color',InstColours(iInst),'linewidth',1);
    
  end
  
    
  
  
  
  %rescale axis to the data
  axis([0 365 0 Max.*1.05])
  
  %replace the time axis with months
  MonthNames = {'J','F','M','A','M','J','J','A','S','O','N','D'};
  set(gca,'xtick',MonthStarts+15,'xticklabel',MonthNames)

  
  %finish up
  ylabel('Total MF [mPa]')
  title(['(',PlotLetters(jBasin),')  ',Basins{iBasin}])
  freezeColors;
  drawnow
end

clear iBasin iInst MonthNames Data Max MonthStarts Colour
clear i iMonth PlotLetters PlotPos SmoothDays
clear TimeSeries Storms jBasin Maps

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%colour bars. bit complicated and cheaty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.05], [0.05 0.1], [0.01 0.05]);
subplot(6,4,[12])


hold on


iBasin = 0;
for jBasin=1:1:7;
  if jBasin == 2; continue; end
  iBasin = iBasin+1;
  
  
  %load regional colour map
  CT = cbrewer('seq',Colours{jBasin},128); CT = CT(16:numel(CT(:,1)),:); 
  colormap(CT); 
  
  %produce a bar for this region
  x = MinStore(jBasin):(MaxStore(jBasin)-MinStore(jBasin))./128:MaxStore(jBasin);
  y = [-1,1].*0.4+BarY(iBasin); 
  Bar = [x;x];
  
  %plot the bar
  pcolor(x,y,Bar); shading flat

  %plot the lines around it
  plot([1,1].*MinStore(jBasin),y,'k-')
  plot([1,1].*MaxStore(jBasin),y,'k-')
  plot([MinStore(jBasin),MaxStore(jBasin)],[1,1].*y(1),'k-')
  plot([MinStore(jBasin),MaxStore(jBasin)],[1,1].*y(2),'k-')
  
  %label with region
  
  %finish up bar
  caxis([MinStore(jBasin) MaxStore(jBasin)])
  freezeColors;
  drawnow  
  
end

% 
% %lines on axis
% plot([1,1].*0.0,[-999,999],'k-');
% plot([1,1].*1.3,[-999,999],'k-');
for x=0.1:0.1:0.6; plot([1,1].*x,[-999,999],'k--');end

%bottom line
plot([-999 999],[6.5,6.5],'k-')

%finish up
ylim([0.5 6.5])
xlim([0.3 0.6])
set(gca,'ytick',[],'ycolor','w','ydir','reverse')
set(gca,'xaxis','top','xtick',0.3:0.1:0.6)
xlabel('Gridbox Mean MF [mPa]')
  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%instrument key
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MLS
plot([0.35,0.40],[1,1].*8,'r-','clipping','off','linewidth',2)
text(0.42,8,'MLS','horizontalalignment','left','color','r','clipping','off','fontsize',16)

%SABER
plot([0.35,0.40],[1,1].*9.5,'k-','clipping','off','linewidth',2)
text(0.42,9.5,'SABER','horizontalalignment','left','color','k','clipping','off','fontsize',16)

%HIRDLS
plot([0.35,0.40],[1,1].*11,'b-','clipping','off','linewidth',2)
text(0.42,11,'HIRDLS','horizontalalignment','left','color','b','clipping','off','fontsize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%done!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% export_fig('global_effects','-pdf','-native')