ncdisp('/project/pdpanalysis/nish/data/ERA5_Avg_SLP/ERA5_Avg_SLP.nc')

times = ncread('/project/pdpanalysis/nish/data/ERA5_Avg_SLP/ERA5_Avg_SLP.nc', 'time');
startday = '1900-01-01 00:00:00.0';
start = datetime(startday,'InputFormat','yyyy-MM-dd HH:mm:ss.s');
data_day = start + hours(times);

data_day = reshape(data_day, 12, 38);
data_day(2:11, :) = [];
data_day = reshape(data_day, 2*38, 1);
data_day = data_day(2:75, :);


lons = ncread('/project/pdpanalysis/nish/data/ERA5_Avg_SLP/ERA5_Avg_SLP.nc', 'longitude');
lats = ncread('/project/pdpanalysis/nish/data/ERA5_Avg_SLP/ERA5_Avg_SLP.nc', 'latitude');

slp = ncread('/project/pdpanalysis/nish/data/ERA5_Avg_SLP/ERA5_Avg_SLP.nc', 'msl');

%slp runs from 1981 Jan to 2018 Dec

slp = reshape(slp, 1440*261*12, 38);
slp_anom = detrend(slp')';
slp_anom = reshape(slp_anom, 1440*261, 12, 38);
slp_anom_Oct_Nov = slp_anom;
slp_anom_Oct_Nov(:, 2:11, :) = [];
slp_anom_Oct_Nov = reshape(slp_anom_Oct_Nov, 1440*261, 2*38);
slp_anom_Oct_Nov = slp_anom_Oct_Nov(:, 2:75);
slp_anom_Oct_Nov = reshape(slp_anom_Oct_Nov, 1440*261, 2, 37);
Seasonal_slp = mean(slp_anom_Oct_Nov, 2);
Seasonal_slp = reshape(Seasonal_slp, 1440*261*1, 37);
Seasonal_slp = detrend(Seasonal_slp')';

%when processed data runs from 1981-Dec to Jan-2018

load HVD_Cop_LS_EOF_2_Oct.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 7, 39);
KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);
%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 38),1);

KI = KI(:, 1:37);
Var = Seasonal_slp(:, 1:37);

rows = length(Var(1:end,:));
cols = 4;

lag_Cor_mat = zeros(rows, cols);

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

load HVD_Cop_LS_EOF_2_Oct.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 7, 39);
KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);
%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 38),1);

KI = KI(:, 1:37);
Var = Seasonal_slp(:, 1:37);

rows = length(Var(1:end,:));
cols = 4;

lead_Cor_mat = zeros(rows, cols);

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

load HVD_Cop_LS_EOF_2_Oct.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 7, 39);
KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);
%KI = mean(reshape(HVD_Cop_LS_EOF_2, 6, 38),1);

KI = KI(:, 1:37);
Var = Seasonal_slp(:, 1:37);

y = std(Var',1)'*ones(1,length(KI));
x = mean(Var',1)'*ones(1,length(KI));

T_norm_LS = (KI(:,1:end)-mean(KI(:,1:end)))/std(KI(:,1:end)')';
A_norm_LS = (Var-x)./y;   %these two lines calculate the standardized or normalized anomalies of T and A respectively

T_norm_LS = detrend(T_norm_LS')';
A_norm_LS = detrend(A_norm_LS')';

Zero_lag_corr = A_norm_LS*T_norm_LS'/length(KI);
Zero_lag_reg = Var*T_norm_LS'/length(KI);

GOA_MLD_Lead_Lag = [lead_Reg_mat(:,4) lead_Reg_mat(:,3) lead_Reg_mat(:,2) lead_Reg_mat(:,1)...
    Zero_lag_reg lag_Reg_mat];
	

Blob_Tau = Seasonal_slp(:,31:34); %selects 2011 Dec 2012Jan and next three years

GOA_MLD_regs = GOA_MLD_Lead_Lag(:, 3:6);

Blob_Tau(isnan(Blob_Tau)) = 0;
GOA_MLD_regs(isnan(GOA_MLD_regs)) = 0;

corvar1 = reshape(Blob_Tau, 1440, 261, 4);
corvar1 = corvar1(:, 1:181, :);
temp_1 = corvar1(1:242,:,:);
temp_2 = corvar1(1401:1440,:,:);
corvar1 = [temp_2; temp_1];
corvar1 = reshape(corvar1, 282*181, 4);

corvar2 = reshape(GOA_MLD_regs, 1440, 261, 4);
corvar2 = corvar2(:, 1:181, :);
temp_1 = corvar2(1:242,:,:);
temp_2 = corvar2(1401:1440,:,:);
corvar2 = [temp_2; temp_1];
corvar2 = reshape(corvar2, 282*181, 4);

for i = 1:4
	a = corvar1(:,i);
	b = corvar2(:,i);
	R = corr2(a, b);
	
	SST_cor_coeff(i) = R;
	
end

close

%% below chunk computes the significance at 95%
rng(1)

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/37; % convert to power spectral density

Var = Seasonal_slp(:, 1:37);
rand_t = randn(37, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

for i = 1:4
	
	rand_t(end, :) = [];
	Var(:, 1) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    lead_sig(:,i) = temp3;

end

Var = Seasonal_slp(:, 1:37);
rand_t = randn(37, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

for i = 1:4
	
	rand_t(1, :) = [];
	Var(:, end) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    lag_sig(:,i) = temp3;

end

Var = Seasonal_slp(:, 1:37);
rand_t = randn(37, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	nolag_sig = temp2(:, 950);
	
sig_matrix = [lead_sig(:,4) lead_sig(:,3) lead_sig(:,2) lead_sig(:,1)...
    nolag_sig lag_sig];
	
fig_sig_matrix = sig_matrix(:, 3:6);	
	
mask = (abs(GOA_MLD_regs) > fig_sig_matrix);
mask = double(mask);
mask(mask == 0) = NaN;

lons_GOA = lons;
lats_GOA = lats;

lats_GOA_sparse = lats_GOA(1:10:261);
lons_GOA_sparse = lons_GOA(1:10:1440);

%LA = reshape(double(lats*ones([1 1440])), 261, 1440);
%LO = reshape(double(lons*ones([1 261])), 1440, 261);

LA = reshape(double(lats_GOA_sparse*ones([1 144])), 27, 144);
LO = reshape(double(lons_GOA_sparse*ones([1 27])), 144, 27);

%% Plot figure
%tiledlayout(4,2);
%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 24 6.06]);
%set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 12 3.8]);
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 13.3 5.78]);


LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];

lons_GOA = lons;
lats_GOA = lats;

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
        subplot('position', [0.05+(k-1)*0.22 0.49 0.19 0.48])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 9, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 200;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(GOA_MLD_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 1440, 261);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask(:,k));
plot_data1 = reshape(temp6, 1440, 261);
plot_data1_sparse = plot_data1(1:10:1440, 1:10:261);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 2.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(62,-134, num2str(round(SST_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',12)
    textm(25,173, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['Year ' clag],'FontSize',14);
   
    tightmap
   gridm on 
   
   hold on;

% la = [30 30 60 60 30]; 
% lo = [180 235 235 180 180];
% 
% %p = projcrs(53009, 'Authority','ESRI');
% [x,y] = mfwdtran (la,lo); 
% line(x,y,'color','#636363', 'LineStyle', '--', 'Linewidth', 1.5)

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
%subplot('position',[0.85 0.572 0.088 0.32])
subplot('position',[0.92 0.535 0.088 0.392])
%subplot(5,2,9)
%subplot;

limit = 200; 
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
h.FontSize = 10;
h.Ticks = [-200 -100 0 100 200];
h.LineWidth = 0.05;
xlabel(h, 'Pa')

        
panellabel = { 'Dec-Jan 2011/2012' 'Dec-Jan 2012/2013' 'Dec-Jan 2013/2014' 'Dec-Jan 2014/2015'};

LatLim = [ 20 65];  %this is where you set the lat/lon limits you want to plot
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
		%subplot('position', [0.05+(k-1)*0.22 0.49 0.19 0.48])
        subplot('position', [0.05+(k-1)*0.22 0.01 0.19 0.48])
		
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 9, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 1000;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(Blob_Tau(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 1440, 261);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;
%textm(60,-125, num2str(SST_cor_coeff(k)),'FontWeight','Bold')

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    

    clag = int2str(o);
    title(panellabel(k),'FontSize',12);
    textm(25,173, tilelab(k),'FontWeight','Bold', 'FontSize',11)
   
    tightmap
   gridm on 

        end


%subplot('position',[0.51 -0.001 0.48 0.05])
%subplot('position',[0.02 0.1 0.088 0.3])
subplot('position',[0.92 0.055 0.088 0.392])


%subplot(5,2,10)
limit = 1000; 
    caxis([-limit limit]); % set the limits of the colorbar
%h=colorbar('Location','southoutside'); 
h=colorbar('Location','west');
h.FontWeight = 'bold';
h.FontSize = 10;
h.Ticks = [-1000 -500 0 500 1000];
h.LineWidth = 0.05;
xlabel(h, 'Pa')
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
 %print('Nish_100522_KE_MHW_MS_SLP_Fig','-dpng', '-r500');
%print('Nish_101322_KE_MHW_MS_SLP_Fig','-dpng', '-r500');
%print('Nish_102722_KE_MHW_MS_SLP_Fig','-dpng', '-r500');
%print('Nish_021423_KE_MHW_MS_SLP_Fig','-dpng', '-r500');
%print('Nish_051423_KE_MHW_MS_SLP_Fig','-dpng', '-r500');
%print('Nish_060223_KE_MHW_MS_SLP_Fig','-dpng', '-r500');
print('Nish_061423_KE_MHW_MS_SLP_Fig','-dpng', '-r500');