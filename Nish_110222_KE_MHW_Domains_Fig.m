%initial figure including domain areas for the paper
lons = ncread('/project/pdpanalysis/nish/data/Copernicus_SST/2000.nc', 'lon');
lats = ncread('/projectnb/pdpanalysis/data/Copernicus_SST_Jan2021/East/2000_E.nc', 'lat');


load /project/pdpanalysis/nish/data/Cop_means_1981_2020.mat;
Cop_SST = reshape(Cop_means_1981_2020, 2400*1800*7, 40);
Cop_SST(isnan(Cop_SST)) = 0;
sst_anom = detrend(Cop_SST')'; 						%calculate the anomalies
data_anom = reshape(sst_anom, 2400, 1800, 7, 40);
data_anom_spr = data_anom(:, :, 4:5, 34);
Seasonal_sst = mean(data_anom_spr, 3);

%%
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 13.8 5.78]);

LatLim = [25 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [120 240];

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

h = worldmap(LatLim,LonLim);
	
setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 14, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
           'MLineLocation', 20,...
          'MLabelLocation',20,...
          'PLineLocation', 10,...
          'PLabelLocation',10);
framem('FLineWidth',1);
	
	
limit = 3;  %here you set the max/min value you want shaded
contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
cmap_mod = flipud(precip_cmap);
mycolormap = customcolormap(linspace(0,1,11), {...
    '#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
colormap(mycolormap);
%colormap(cmap_mod)
caxis([-limit limit]); % set the limits of the colorbar
	
	
temp5 = double(Seasonal_sst);   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 2400, 1800);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats), double(lons), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;


load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth', 1)
    
 tightmap
   gridm on 
   
   
hold on;

la = [30 30 42 42 30]; 
lo = [140 172 172 140 140];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','#000000', 'Linewidth', 2.5)

hold on;

%la = [20 20 65 65 20]; 
%lo = [170 240 240 170 170];

%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','k', 'Linewidth', 1.5)

hold on;

la = [30 30 60 60 30]; 
lo = [200 240 240 200 200];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','#238b45', 'Linewidth', 2.5)
   
   %subplot('position',[0.92 0.055 0.088 0.392])

   hold on
%subplot(5,2,10)
limit = 3; 
    caxis([-limit limit]); % set the limits of the colorbar
%h=colorbar('Location','southoutside'); 
h=colorbar('Location','eastoutside');
h.FontWeight = 'bold';
h.FontSize = 14;
h.LineWidth = 0.05;
h.Ticks = [-3 -2 -1 0 1 2 3];
%xlabel(h, '^{o}C')
%h.Label = 'mol m^{-3}';
%cbfreeze(h)
%x = get(cbhandle,'Position');
%x(2) = x(2)-0.125;
%x(4) = 1.2*x(4);
%set(cbhandle,'Position',x);
%set(cbhandle,'XTick',[-3:0.5:3]);   %here you set the max/min value you want labeled
%set(cbhandle,'XTickLabel',[-2:1:2]);  %here you set the max/min value you want labeled
%set(cbhandle,'YAxisLocation','right');
%h.XTicks = [-3:0.5:3];
%h.TickLabels = [-3:1:3];
%cmap_mod = flipud(precip_cmap);
%colormap(cmap_mod)
colormap(customcolormap);
axis 'off';
%axis('off')

print('Nish_110222_KE_MHW_Domains_Fig','-dpng', '-r500');