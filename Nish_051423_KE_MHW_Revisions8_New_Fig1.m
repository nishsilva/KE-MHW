lons = ncread('/project/pdpanalysis/nish/data/Copernicus_SST/2000.nc', 'lon');
lats = ncread('/projectnb/pdpanalysis/data/Copernicus_SST_Jan2021/East/2000_E.nc', 'lat');


load /project/pdpanalysis/nish/data/Cop_means_1981_2020.mat;
Cop_SST = reshape(Cop_means_1981_2020, 2400*1800*7, 40);
Cop_SST(isnan(Cop_SST)) = 0;
sst_anom = detrend(Cop_SST')'; 						%calculate the anomalies
data_anom = reshape(sst_anom, 2400, 1800, 7, 40);
data_anom_spr = data_anom(:, :, 4:5, :);
Seasonal_sst = mean(data_anom_spr, 3);

Seasonal_sst = Seasonal_sst(1001:end, 400:1301, :, :); %subset GOA
Seasonal_sst = reshape(Seasonal_sst, 1400*902, 1*40);

load HVD_Cop_LS_EOF_2_2020.mat

HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);

%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 40),1);
Var = Seasonal_sst;

rows = length(Var(1:end,:));
cols = 4;

lag_Cor_mat = zeros(rows, cols);
lag_Reg_mat = zeros(rows, cols);

for i = 1:4
	
	KI(:,end) = [];
	Var(:,1) = [];
	
	y = std(Var',1)'*ones(1,length(KI));
	x = mean(Var',1)'*ones(1,length(KI));

	T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
	A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

	T_norm_LS = detrend(T_norm_LS')';
	A_norm_LS = detrend(A_norm_LS')';
    
    A_corr_LS = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
	lag_Cor_mat(:,i) = A_corr_LS;
	A_regress = Var*T_norm_LS'/length(KI);
	lag_Reg_mat(:,i) = A_regress;
	
end

KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);
%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 40),1);
Var = Seasonal_sst;

rows = length(Var(1:end,:));
cols = 4;

for i = 1:4
	
	KI(:,1) = [];
	Var(:,end) = [];
	
	y = std(Var',1)'*ones(1,length(KI));
	x = mean(Var',1)'*ones(1,length(KI));

	T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
	A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

	T_norm_LS = detrend(T_norm_LS')';
	A_norm_LS = detrend(A_norm_LS')';
    
    A_corr_LS = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
	lead_Cor_mat(:,i) = A_corr_LS;
	A_regress = Var*T_norm_LS'/length(KI);
	lead_Reg_mat(:,i) = A_regress;
	
end

KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);
%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 40),1);
Var = Seasonal_sst;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';

Zero_lag_corr = A_norm_LS*T_norm_LS'/length(KI);
Zero_lag_reg = Var*T_norm_LS'/length(KI);


GOA_Blob_SST_Lead_Lag = [lead_Reg_mat(:,4) lead_Reg_mat(:,3) lead_Reg_mat(:,2) lead_Reg_mat(:,1)...
    Zero_lag_reg lag_Reg_mat];
	
GOA_Seasonal_sst = Seasonal_sst(:,31:35);
GOA_Blob_SST_Lead_Lag(isnan(GOA_Blob_SST_Lead_Lag)) = 0;
GOA_Blob_SST_regs = GOA_Blob_SST_Lead_Lag(:, 3:7);

for i = 1:4
	a = GOA_Seasonal_sst(:,i);
	b = GOA_Blob_SST_regs(:,i);
	R = corr2(a, b);
	
	BLOB_SST_cor_coeff(i) = R;
	
end

lons_GOA = lons(1001:end,:);
lats_GOA = lats(400:1301, :);

%% below chunk computes the significance at 95%
rng(1);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/40; % convert to power spectral density

Var = Seasonal_sst;
rand_t = randn(40, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

for i = 1:4
	
	rand_t(end, :) = [];
	Var(:, 1) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    lead_sig(:,i) = temp3;

end

Var = Seasonal_sst;
rand_t = randn(40, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

for i = 1:4
	
	rand_t(1, :) = [];
	Var(:, end) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    lag_sig(:,i) = temp3;

end

Var = Seasonal_sst;
rand_t = randn(40, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	nolag_sig = temp2(:, 950);
	
sig_matrix = [lead_sig(:,4) lead_sig(:,3) lead_sig(:,2) lead_sig(:,1)...
    nolag_sig lag_sig];
	
fig_sig_matrix = sig_matrix(:, 3:7);	
	
mask = (abs(GOA_Blob_SST_regs) > fig_sig_matrix);
mask = double(mask);
mask(mask == 0) = NaN;

lats_GOA_sparse = lats_GOA(1:30:902);
lons_GOA_sparse = lons_GOA(1:30:1400);

%LA = reshape(double(lats_GOA*ones([1 1400])), 902, 1400);
%LO = reshape(double(lons_GOA*ones([1 902])), 1400, 902);

LA = reshape(double(lats_GOA_sparse*ones([1 47])), 31, 47);
LO = reshape(double(lons_GOA_sparse*ones([1 31])), 47, 31);

%% Plot figure
%%
%tiledlayout(4,2);
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 13.3 5.78]);

LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = { '(a)' '(b)' '(c)' '(d)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 0;
        hi = [1 3 5 7];
        for k = 1:4
                        
           ilag = ilag+1;
			o = ilag - 3;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        %subplot('position', [0.11+(k-1)*0.2 0.49 0.19 0.48])
        subplot('position', [0.05+(k-1)*0.22 0.49 0.19 0.48])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
           'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 0.5;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(GOA_Blob_SST_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 1400, 902);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask(:,k));
plot_data1 = reshape(temp6, 1400, 902);
plot_data1_sparse = plot_data1(1:30:1400, 1:30:902);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 2.5);


load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(BLOB_SST_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',14)
    textm(63.5,171, tilelab(k),'FontWeight','Bold', 'FontSize',10)

    clag = int2str(o);
    title(['Year ' clag],'FontSize',16);
   
    tightmap
   gridm on 

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.01 0.58 0.1 0.3])
subplot('position',[0.92 0.535 0.088 0.392])
%subplot(5,2,9)
%subplot;

limit = 0.5; 
    caxis([-limit limit]); % set the limits of the colorbar
h=colorbar('Location', 'west'); 
%h = colorbar;
%h.Position = [0.1 0.001 0.48 0.14];
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
h.FontWeight = 'bold';
h.FontSize = 12;
h.LineWidth = 0.05;
h.Ticks = [-0.5 -0.25 0 0.25 0.5];
%xlabel(h, '^{o}C')

        
panellabel = { 'Feb-Mar 2012' 'Feb-Mar 2013' 'Feb-Mar 2014' 'Feb-Mar 2015'};

LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%nexttlie
%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };
tilelab = { '(e)' '(f)' '(g)' '(h)'};

%h = figure;
 %       set(h, 'Color', 'w');
		
		
		ilag = 0;
        li = [2 4 6 8];
        for k=1:4
			ilag = ilag+1;
			o = ilag - 2;
		
		%subplot(5, 2, li(k))
        %subplot('position',[0.51 .74-(k-1)*.25 0.48 .23])
		%subplot('position', [0.11+(k-1)*0.22 0.01 0.19 0.48])
        subplot('position', [0.05+(k-1)*0.22 0.01 0.19 0.48])
		
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

temp5 = double(GOA_Seasonal_sst(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 1400, 902);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;
%textm(60,-125, num2str(SST_cor_coeff(k)),'FontWeight','Bold')

load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth', 0.5)
    

    clag = int2str(o);
    title(panellabel(k),'FontSize',16);
    textm(63.5,171, tilelab(k),'FontWeight','Bold', 'FontSize',10)
   
    tightmap
   gridm on 

        end


%subplot('position',[0.51 -0.001 0.48 0.05])
%subplot('position',[0.01 0.1 0.088 0.3])
subplot('position',[0.92 0.055 0.088 0.392])

%subplot(5,2,10)
limit = 2; 
    caxis([-limit limit]); % set the limits of the colorbar
%h=colorbar('Location','southoutside'); 
h=colorbar('Location','west');
h.FontWeight = 'bold';
h.FontSize = 12;
h.LineWidth = 0.05;
h.Ticks = [-2 -1 0 1 2];
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
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')

%text = 'char(181) mol m^{-3}';

%%
%before exporting the figure make sure to delete latitude labels manually
 %print('Nish_100422_KE_MHW_MS_Fig1','-dpng', '-r500');
 %print('Nish_101222_KE_MHW_MS_Fig1','-dpng', '-r500');
%print('Nish_021423_KE_MHW_MS_Fig1','-dpng', '-r500');
print('Nish_051423_KE_MHW_MS_Fig1','-dpng', '-r500');
%% GOA EOF Variability

%% here we get the blob area and analyze the leading variations of blob area SSTs

load /project/pdpanalysis/nish/data/Cop_means_1981_2020.mat;
Cop_SST = reshape(Cop_means_1981_2020, 2400*1800*7, 40);
Cop_SST(isnan(Cop_SST)) = 0;
sst_anom = detrend(Cop_SST')'; 						%calculate the anomalies
data_anom = reshape(sst_anom, 2400, 1800, 7, 40);
data_anom_spr = data_anom(:, :, 4:5, :);
Seasonal_sst = mean(data_anom_spr, 3);

%Seasonal_sst = Season al_sst(1601:end, 600:1201, :, :); %subset blob region 30-60N 160-120W
%Seasonal_sst = reshape(Seasonal_sst, 800, 602, 1*40); 

Seasonal_sst = Seasonal_sst(1001:end, 400:1301, :, :); %subset GOA 20-65N & 170-60W
Seasonal_sst = reshape(Seasonal_sst, 1400, 902, 1*40);

[eofmap_large, pc_large, expvar_large] = eof(Seasonal_sst);

%% Explained variance

%bar(expvar_large(:,1:10));
%barvalues;

set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 6 4]);
Y = round(expvar_large(:, 1:10));
bar(expvar_large(:,1:10));
hold on 
plot(expvar_large(:, 1:10), 'LineWidth', 2);
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center');
%text(1:length(Y),Y,num2str(cumsum(Y')),'vert','top','horiz','center');
title('GOA EOF Explained Variance (Copernicus Data)');
%%
plot(pc_large(1,:));
plot(pc_large(2,:));

Blob_SST_PC1_TS = detrend((normalize(pc_large(1,:).*(-1)))')';
Blob_SST_PC2_TS = detrend((normalize(pc_large(2,:).*(1)))')';

save Blob_SST_PC1_TS.mat Blob_SST_PC1_TS; 
save Blob_SST_PC2_TS.mat Blob_SST_PC2_TS;


lons = ncread('/project/pdpanalysis/nish/data/Copernicus_SST/2000.nc', 'lon');
lats = ncread('/projectnb/pdpanalysis/data/Copernicus_SST_Jan2021/East/2000_E.nc', 'lat');

lons_GOA = lons(1001:end,:);
lats_GOA = lats(400:1301, :);


%% PC2 & GOA SST regression

load /project/pdpanalysis/nish/data/Cop_means_1981_2020.mat;
Cop_SST = reshape(Cop_means_1981_2020, 2400*1800*7, 40);
Cop_SST(isnan(Cop_SST)) = 0;
sst_anom = detrend(Cop_SST')'; 						%calculate the anomalies
data_anom = reshape(sst_anom, 2400, 1800, 7, 40);
data_anom_spr = data_anom(:, :, 4:5, :);
Seasonal_sst = mean(data_anom_spr, 3);

Seasonal_sst = Seasonal_sst(1001:end, 400:1301, :, :); %subset GOA
Seasonal_sst = reshape(Seasonal_sst, 1400*902, 1*40);

KI = Blob_SST_PC2_TS;
Var = Seasonal_sst;

rows = length(Var(1:end,:));
cols = 4;

lag_Cor_mat = zeros(rows, cols);
lag_Reg_mat = zeros(rows, cols);

for i = 1:4
	
	KI(:,end) = [];
	Var(:,1) = [];
	
	y = std(Var',1)'*ones(1,length(KI));
	x = mean(Var',1)'*ones(1,length(KI));

	T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
	A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

	T_norm_LS = detrend(T_norm_LS')';
	A_norm_LS = detrend(A_norm_LS')';
    
    A_corr_LS = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
	lag_Cor_mat(:,i) = A_corr_LS;
	A_regress = Var*T_norm_LS'/length(KI);
	lag_Reg_mat(:,i) = A_regress;
	
end


KI = Blob_SST_PC2_TS;
Var = Seasonal_sst;

rows = length(Var(1:end,:));
cols = 4;

lead_Cor_mat = zeros(rows, cols);
lead_Reg_mat = zeros(rows, cols);

for i = 1:4
	
	KI(:,1) = [];
	Var(:,end) = [];
	
	y = std(Var',1)'*ones(1,length(KI));
	x = mean(Var',1)'*ones(1,length(KI));

	T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
	A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

	T_norm_LS = detrend(T_norm_LS')';
	A_norm_LS = detrend(A_norm_LS')';
    
    A_corr_LS = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
	lead_Cor_mat(:,i) = A_corr_LS;
	A_regress = Var*T_norm_LS'/length(KI);
	lead_Reg_mat(:,i) = A_regress;
	
end

KI = Blob_SST_PC2_TS;
Var = Seasonal_sst;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';

Zero_lag_corr = A_norm_LS*T_norm_LS'/length(KI);
Zero_lag_reg = Var*T_norm_LS'/length(KI);

GOA_SST_Lead_Lag = [lead_Reg_mat(:,4) lead_Reg_mat(:,3) lead_Reg_mat(:,2) lead_Reg_mat(:,1)...
    Zero_lag_reg lag_Reg_mat];


GOA_Blob_SST_Lead_Lag = [lead_Reg_mat(:,4) lead_Reg_mat(:,3) lead_Reg_mat(:,2) lead_Reg_mat(:,1)...
    Zero_lag_reg lag_Reg_mat];
	
GOA_SST_regs = GOA_SST_Lead_Lag(:, 3:6);

lons = ncread('/project/pdpanalysis/nish/data/Copernicus_SST/2000.nc', 'lon');
lats = ncread('/projectnb/pdpanalysis/data/Copernicus_SST_Jan2021/East/2000_E.nc', 'lat');

lons_GOA = lons(1001:end,:);
lats_GOA = lats(400:1301, :);

%%

%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 12 6]);
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 13.3 5.78]);


LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = { '(a)' '(b)' '(c)' '(d)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 0;
        hi = [1 3 5 7];
        for k = 1:4
                        
           ilag = ilag+1;
			o = ilag - 3;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        %subplot('position', [0.11+(k-1)*0.2 0.49 0.19 0.48])
        subplot('position', [0.05+(k-1)*0.22 0.49 0.19 0.48])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 0.5;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(GOA_SST_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 1400, 902);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth', 0.5)
    %textm(61,-135, num2str(round(SST_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',14)
    textm(63.5,171, tilelab(k),'FontWeight','Bold', 'FontSize',10)

    clag = int2str(o);
    title(['Year ' clag],'FontSize',16);
   
    tightmap
   gridm off 
   
la = [30 30 60 60 30]; 
lo = [200 240 240 200 200];

%p = projcrs(53009, 'Authority','ESRI');
%[x,y] = mfwdtran (la,lo); 
%line(x,y,'color','k', 'Linewidth', 2)    

        end
        
hold on

%textm(63.5,230, 'GOA PC1','FontWeight','Bold', 'FontSize',10);
%annotation('textbox', [0.01, 0.9, 0.09, 0.06], 'string', 'GOA EOF 2',...
%    'FitBoxToText','on', 'FontWeight','Bold', 'FontSize',14);
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
subplot('position',[0.92 0.535 0.088 0.392])
%subplot(5,2,9)
%subplot;

limit = 0.5; 
    caxis([-limit limit]); % set the limits of the colorbar
h=colorbar('Location', 'west'); 
%h = colorbar;
%h.Position = [0.1 0.001 0.48 0.14];
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
h.FontWeight = 'bold';
h.FontSize = 16;
h.LinwWidth = 0.05;
%xlabel(h, '^{o}C')

%%


KI = Blob_SST_PC2_TS;
Var = Seasonal_sst;

rows = length(Var(1:end,:));
cols = 4;

lag_Cor_mat = zeros(rows, cols);
lag_Reg_mat = zeros(rows, cols);

for i = 1:4
	
	KI(:,end) = [];
	Var(:,1) = [];
	
	y = std(Var',1)'*ones(1,length(KI));
	x = mean(Var',1)'*ones(1,length(KI));

	T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
	A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

	T_norm_LS = detrend(T_norm_LS')';
	A_norm_LS = detrend(A_norm_LS')';
    
    A_corr_LS = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
	lag_Cor_mat(:,i) = A_corr_LS;
	A_regress = Var*T_norm_LS'/length(KI);
	lag_Reg_mat(:,i) = A_regress;
	
end

KI = Blob_SST_PC2_TS;
Var = Seasonal_sst;

rows = length(Var(1:end,:));
cols = 4;

lead_Cor_mat = zeros(rows, cols);
lead_Reg_mat = zeros(rows, cols);

for i = 1:4
	
	KI(:,1) = [];
	Var(:,end) = [];
	
	y = std(Var',1)'*ones(1,length(KI));
	x = mean(Var',1)'*ones(1,length(KI));

	T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
	A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

	T_norm_LS = detrend(T_norm_LS')';
	A_norm_LS = detrend(A_norm_LS')';
    
    A_corr_LS = A_norm_LS*T_norm_LS'/length(KI);  % This calculates the correlation values; the resulting vector is a map,  A_corr(x)
	lead_Cor_mat(:,i) = A_corr_LS;
	A_regress = Var*T_norm_LS'/length(KI);
	lead_Reg_mat(:,i) = A_regress;
	
end

KI = Blob_SST_PC2_TS;
Var = Seasonal_sst;

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';

Zero_lag_corr = A_norm_LS*T_norm_LS'/length(KI);
Zero_lag_reg = Var*T_norm_LS'/length(KI);

GOA_SST_Lead_Lag = [lead_Reg_mat(:,4) lead_Reg_mat(:,3) lead_Reg_mat(:,2) lead_Reg_mat(:,1)...
    Zero_lag_reg lag_Reg_mat];


GOA_Blob_SST_Lead_Lag = [lead_Reg_mat(:,4) lead_Reg_mat(:,3) lead_Reg_mat(:,2) lead_Reg_mat(:,1)...
    Zero_lag_reg lag_Reg_mat];
	
GOA_SST_regs = GOA_SST_Lead_Lag(:, 3:6);

panellabel = { '2012' '2013' '2014' '2015'};

LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];
%%
load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%nexttlie
%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };
tilelab = { '(e)' '(f)' '(g)' '(h)'};

%h = figure;
 %       set(h, 'Color', 'w');
		
		
		ilag = 0;
        li = [2 4 6 8];
        for k=1:4
			ilag = ilag+1;
			o = ilag - 3;
		
		%subplot(5, 2, li(k))
        %subplot('position',[0.51 .74-(k-1)*.25 0.48 .23])
		subplot('position', [0.11+(k-1)*0.2 0.01 0.19 0.48])
       
		
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'off','MeridianLabel', 'off', ...
          'FontSize', 8, ...
          'FontWeight', 'normal', ...
          'LabelFormat', 'none');
	framem('FLineWidth',0.5);
		
    limit = 0.5;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(GOA_SST_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 1400, 902);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;
%textm(60,-125, num2str(SST_cor_coeff(k)),'FontWeight','Bold')

load coast
    geoshow(lat, long, 'DisplayType','polygon','FaceColor','white', 'LineWidth', 0.5)
    

    clag = int2str(o);
    title(['Year ' clag],'FontSize',16);
    textm(63.5,171, tilelab(k),'FontWeight','Bold', 'FontSize',10)
   
    tightmap
   gridm off 
   
la = [30 30 60 60 30]; 
lo = [200 240 240 200 200];

%p = projcrs(53009, 'Authority','ESRI');
[x,y] = mfwdtran (la,lo); 
line(x,y,'color','k', 'Linewidth', 2)   

        end

textm(63.5,230, 'GOA PC2','FontWeight','Bold', 'FontSize',10)

%subplot('position',[0.51 -0.001 0.48 0.05])
subplot('position',[0.02 0.1 0.088 0.3])

%subplot(5,2,10)
limit = 0.5; 
    caxis([-limit limit]); % set the limits of the colorbar
%h=colorbar('Location','southoutside'); 
h=colorbar('Location','westoutside');
h.FontWeight = 'bold';
h.FontSize = 16;
h.LineWidth = 0.05;
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
cmap_mod = flipud(precip_cmap);
colormap(cmap_mod)
axis('off')

%text = 'char(181) mol m^{-3}';

%% Plot KI & Blob EOF

load Blob_SST_PC2_TS.mat;
a = nan;

GOA_PC2 = horzcat(a, Blob_SST_PC2_TS);

GOA_PC2 = GOA_PC2(:, 1:40);

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
HVD_Cop_LS_EOF_2 = HVD_Cop_LS_EOF_2(2:3, :);

subplot('position',[0.02 0.1 0.6 0.4])

KI = mean(reshape(HVD_Cop_LS_EOF_2, 2, 40),1);

[c,lags] = xcorr(KI, GOA_PC2, 5, 'normalized');
stem(lags,c, 'filled', 'LineWidth', 2);

Years = 1981:2020;

%tiledlayout(1, 4);
%nexttile([1,3])
x = plot(KI, '--', 'Color', '#0072BD', 'LineWidth',2.5);
yline(0,'--', 'LineWidth',1.5);
xlim([0 38]);
ylim([-3 3]);
hold on
a = plot(GOA_PC2, '--', 'Color', '#D95319', 'LineWidth',2.5);
%b = plot(Seasonal_KEE_EOF_3, 'Color', '#33a02c', 'LineWidth',1.5);
%c = plot(Seasonal_CP_EOF_2, 'Color', '#3182bd', 'LineWidth',1.5);
set(gca,'xtick',(1:5:40),'xticklabel',Years(1:5:40));
ax = gca;
ax.FontSize = 14; 
%ax.LineWidth = 2;
ax.FontWeight = 'bold';
%grid on
%legend([x, a, b, c],{'PDP','KE Index', 'CP EOF2', 'KE EOF 3'},...
%    'Location','NorthEast');
%legend([b, c],{'KE EOF 3', 'CP EOF2'},...
%    'Location','NorthEast');
legend([x, a],{'KI', 'GOA PC 2'},...
    'Location','NorthWest');
legend boxoff;
hold on;
xlabel('Year', 'FontSize', 14, 'FontWeight', 'bold');
%title('(a)', 'FontSize', 22, 'FontWeight', 'bold')	



subplot('position',[0.63 0.1 0.22 0.4])

[c,lags] = xcorr(KI(:,2:end), GOA_PC2(2:end), 5, 'normalized');
stem(lags,c, 'filled', 'LineWidth', 2);
ylim([-0.7 0.7]);
hold on
vline(0, '--', 'Color', '#bdbdbd');
text(-4, -0.55, 'KI leads', 'FontSize',15, 'FontWeight', 'bold')
text(1.5, 0.55, 'GOA PC2 leads', 'FontSize',15, 'FontWeight', 'bold')
set(gca,'YTickLabel',[]);
yyaxis right
ylim([-0.65 0.65]);
ax = gca;
ax.FontSize = 14; 
%ax.LineWidth = 2;
ax.FontWeight = 'bold';
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
xlabel('Lag Year', 'FontSize', 14, 'FontWeight', 'bold');
%title('(b)', 'FontSize', 22, 'FontWeight', 'bold')	