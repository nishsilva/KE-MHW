%%plotting observed atmospheric diagnostics during blob

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 7 9]);


%%
no3 = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'no3');
no3(isnan(no3))=0;

NO3 = reshape(no3, 481*221*12, 28);

NO3_anom = detrend(NO3')';
NO3_anom = reshape(NO3_anom, 481*221, 12, 28);
NO3_anom_spr = NO3_anom(:,2:3,:);
Seasonal_NO3 = mean(NO3_anom_spr, 2);
Seasonal_NO3 = reshape(Seasonal_NO3, 481, 221*1, 28);
Seasonal_NO3 = Seasonal_NO3(201:481, 1:181, :);
Seasonal_NO3 = reshape(Seasonal_NO3, 281*181, 28);

Blob_no3 = Seasonal_NO3(:,22:23);

lons = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'latitude');

lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);



%%
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];


load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = {'(a)' '(b)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 2014;
        hi = [1 3 5 7];
        for k = 1:2
                        
           ilag = ilag+1;
			o = ilag - 1;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        %subplot('position', [0.04+(k-1)*0.22 0.01 0.19 0.48])
        subplot('position', [0.05+(k-1)*0.44 0.67 0.38 0.28])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 4;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(Blob_no3(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['[NO_3^{2-}] Feb/Mar' clag],'FontSize',12,'FontWeight','normal');
   
    tightmap
   gridm on 

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
%subplot('position',[0.00 0.065 0.135 0.36])
%subplot('units','normalized','position',[0.00 0.065 0.135 0.36])
%subplot(5,2,9)
%subplot;
subplot('position',[0.88 0.679 0.18 0.264])

limit = 4; 
    caxis([-limit limit]); % set the limits of the colorbar
h=colorbar('Location', 'west'); 
%h = colorbar;
%h.Position = [0.46 0.07325 0.01 0.355];
cbfreeze(h)
%x = get(cbhandle,'Position');
%x(2) = x(2)-0.125;
%x(4) = 1.2*x(4);
%set(h,'Position',[0.1 0.001 0.48 0.14]);
%set(cbhandle,'YTick',[-5.0:2.0:5.0]);   %here you set the max/min value you want labeled
%set(cbhandle,'YTickLabel',[-5.0:2.0:5.0]);  %here you set the max/min value you want labeled
%set(cbhandle,'YAxisLocation','bottom');
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')
h.FontWeight = 'normal';
h.FontSize = 12;
h.LineWidth = 0.05;
%h.Ticks = [-1.5 -1.0 -0.5 0 0.5 1.0 1.5];


t=get(h,'Limits');
T=linspace(t(1),t(2),5);
set(h,'Ticks',T);
TL=arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
set(h,'TickLabels',TL);
%xlabel(h, 'Nm^{-2}s')

hold on
%%
po4 = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'po4');
po4(isnan(po4))=0;

po4 = reshape(po4, 481*221*12, 28);

po4_anom = detrend(po4')';
po4_anom = reshape(po4_anom, 481*221, 12, 28);
po4_anom_spr = po4_anom(:,2:3,:);
Seasonal_po4 = mean(po4_anom_spr, 2);
Seasonal_po4 = reshape(Seasonal_po4, 481, 221*1, 28);
Seasonal_po4 = Seasonal_po4(201:481, 1:181, :);
Seasonal_po4 = reshape(Seasonal_po4, 281*181, 28);

Blob_po4 = Seasonal_po4(:,22:23);

lons = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'latitude');

lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);

%%
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];


load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = {'(c)' '(d)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 2014;
        hi = [1 3 5 7];
        for k = 1:2
                        
           ilag = ilag+1;
			o = ilag - 1;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        subplot('position', [0.05+(k-1)*0.44 0.37 0.38 0.28])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 0.35;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(Blob_po4(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['[PO_4^{3-}] Feb/Mar' clag],'FontSize',12,'FontWeight','normal');
   
    tightmap
   gridm on 

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
subplot('position',[0.88 0.379 0.18 0.264])
%subplot('units','normalized','position',[0.00 0.065 0.135 0.36])
%subplot(5,2,9)
%subplot;

limit = 0.35; 
    caxis([-limit limit]); % set the limits of the colorbar
h=colorbar('Location', 'west'); 
%h = colorbar;
%h.Position = [0.46 0.07325 0.01 0.355];
cbfreeze(h)
%x = get(cbhandle,'Position');
%x(2) = x(2)-0.125;
%x(4) = 1.2*x(4);
%set(h,'Position',[0.1 0.001 0.48 0.14]);
%set(cbhandle,'YTick',[-5.0:2.0:5.0]);   %here you set the max/min value you want labeled
%set(cbhandle,'YTickLabel',[-5.0:2.0:5.0]);  %here you set the max/min value you want labeled
%set(cbhandle,'YAxisLocation','bottom');
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')
h.FontWeight = 'normal';
h.FontSize = 12;
h.Ticks = [-8 -4 0 4 8];
h.LineWidth = 0.05;

%xlabel(h, 'Nm^{-2}s')

t=get(h,'Limits');
T=linspace(t(1),t(2),5);
set(h,'Ticks',T);
TL=arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
set(h,'TickLabels',TL);

hold on

%%
phyc = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'phyc');
phyc(isnan(phyc))=0;

phyc = reshape(phyc, 481*221*12, 28);

phyc_anom = detrend(phyc')';
phyc_anom = reshape(phyc_anom, 481*221, 12, 28);
phyc_anom_spr = phyc_anom(:,4,:);
Seasonal_phyc = mean(phyc_anom_spr, 2);
Seasonal_phyc = reshape(Seasonal_phyc, 481, 221*1, 28);
Seasonal_phyc = Seasonal_phyc(201:481, 1:181, :);
Seasonal_phyc = reshape(Seasonal_phyc, 281*181, 28);

Blob_phyc = Seasonal_phyc(:,22:23);

lons = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'latitude');

lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);


%%
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];


load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = {'(e)' '(f)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 2014;
        hi = [1 3 5 7];
        for k = 1:2
                        
           ilag = ilag+1;
			o = ilag - 1;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        %subplot('position', [0.04+(k-1)*0.22 0.01 0.19 0.48])
      subplot('position', [0.05+(k-1)*0.44 0.07 0.38 0.28])

        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 2;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(Blob_phyc(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['Phy.Plnktn Apr Year ' clag],'FontSize',12,'FontWeight','normal');
   
    tightmap
   gridm on 

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
%subplot('position',[0.00 0.065 0.135 0.36])
%subplot('units','normalized','position',[0.00 0.065 0.135 0.36])
%subplot(5,2,9)
%subplot;
subplot('position',[0.88 0.079 0.18 0.264])

limit = 2; 
    caxis([-limit limit]); % set the limits of the colorbar
h=colorbar('Location', 'west'); 
%h = colorbar;
%h.Position = [0.46 0.07325 0.01 0.355];
cbfreeze(h)
%x = get(cbhandle,'Position');
%x(2) = x(2)-0.125;
%x(4) = 1.2*x(4);
%set(h,'Position',[0.1 0.001 0.48 0.14]);
%set(cbhandle,'YTick',[-5.0:2.0:5.0]);   %here you set the max/min value you want labeled
%set(cbhandle,'YTickLabel',[-5.0:2.0:5.0]);  %here you set the max/min value you want labeled
%set(cbhandle,'YAxisLocation','bottom');
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')
h.FontWeight = 'normal';
h.FontSize = 12;
h.Ticks = [-8 -4 0 4 8];
h.LineWidth = 0.05;
%xlabel(h, 'Nm^{-2}s')

hold on

t=get(h,'Limits');
T=linspace(t(1),t(2),5);
set(h,'Ticks',T);
TL=arrayfun(@(x) sprintf('%.2f',x),T,'un',0);
set(h,'TickLabels',TL);

%%
%before exporting the figure make sure to delete latitude labels manually
 %print('Nish_100622_KE_MHW_GOA_EOF_Fig','-dpng', '-r500');
% print('Nish_110322_KE_MHW_GOA_Productivity_obs_Fig','-dpng', '-r500');
print('Nish_021423_KE_MHW_GOA_Productivity_obs_Fig','-dpng', '-r500');