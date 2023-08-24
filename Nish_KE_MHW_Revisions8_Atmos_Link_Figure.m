load '/projectnb/pdpanalysis/data/ERA5_lsm_binary.mat'
lsm_binary = reshape(lsm_binary, 1440*261, 1);

%ncread('/projectnb/pdpanalysis/data/ERA5_land_sea_mask.nc')


set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [0 0 7 9]);

%% Windstress
lons = ncread('/project/pdpanalysis/nish/data/ERA5_Windstress/ERA5_Mon_Windstress_2021.nc', 'longitude');
lats = ncread('/project/pdpanalysis/nish/data/ERA5_Windstress/ERA5_Mon_Windstress_2021.nc', 'latitude');

lons1 = lons(1:241);
lons2 = lons(1401:1440);
lons = vertcat(lons2, lons1);
lats = lats(1:181);

turb_stress = ncread('/project/pdpanalysis/nish/data/ERA5_Windstress/ERA5_Mon_Windstress_2021.nc', 'magss');

turb_stress = turb_stress(:,:, 265:end);

turb_stress = reshape(turb_stress, 1440*261*12, 41);
turb_stress_anom = detrend(turb_stress')';
turb_stress_anom = reshape(turb_stress_anom, 1440*261, 12, 41);
turb_stress_anom_Oct_Jan = turb_stress_anom;
turb_stress_anom_Oct_Jan(:, 2:10, :) = [];
turb_stress_anom_Oct_Jan = reshape(turb_stress_anom_Oct_Jan, 1440*261, 3*41);
turb_stress_anom_Oct_Jan = turb_stress_anom_Oct_Jan(:, 2:121);
turb_stress_anom_Oct_Jan = reshape(turb_stress_anom_Oct_Jan, 1440*261, 3, 40);
Seasonal_turb_stress = mean(turb_stress_anom_Oct_Jan, 2);
Seasonal_turb_stress = reshape(Seasonal_turb_stress, 1440*261*1, 40);

load '/projectnb/pdpanalysis/data/ERA5_lsm_binary.mat'
lsm_binary = reshape(lsm_binary, 1440*261, 1);
Seasonal_turb_stress = Seasonal_turb_stress.*lsm_binary;

Seasonal_turb_stress = reshape(Seasonal_turb_stress, 1440, 261, 40);
Seasonal_turb_stress1 = Seasonal_turb_stress(1:241, :, :);
Seasonal_turb_stress2 = Seasonal_turb_stress(1401:1440, :, :);
Seasonal_turb_stress = vertcat(Seasonal_turb_stress2, Seasonal_turb_stress1);
Seasonal_turb_stress = Seasonal_turb_stress(:, 1:181, :);
Seasonal_turb_stress = reshape(Seasonal_turb_stress, 281*181, 40);

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);

KI = KI(:, 1:40);
Var = Seasonal_turb_stress(:, 1:40);

temp1 = NaN(1);
temp2 = NaN(size(Seasonal_turb_stress,1),1);

KI = [KI temp1];
Var = [temp2 Var];

rows = length(Var(1:end,:));
cols = 2;

lag_Cor_mat = zeros(rows, cols);
lag_Reg_mat = zeros(rows, cols);

for i = 1:2
	
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

windstress_regs = lag_Reg_mat;
Blob_Tau = Seasonal_turb_stress(:,33:34);

Blob_Tau(isnan(Blob_Tau)) = 0;
windstress_regs(isnan(windstress_regs)) = 0;

%load '/projectnb/pdpanalysis/data/ERA5_lsm_binary.mat'
%lsm_binary = reshape(lsm_binary, 1440*261, 1);
%windstress_regs = windstress_regs.*lsm_binary;
%Blob_Tau = Blob_Tau.*lsm_binary;

for i = 1:2
	a = Blob_Tau(:,i);
	b = windstress_regs(:,i);
	R = corr2(a, b);
	
	WS_cor_coeff(i) = R;
	
end

KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/40; % convert to power spectral density

Var = Seasonal_turb_stress(:, 1:40);
rand_t = randn(40, 1000);

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

temp1 = NaN(1, 1000);
rand_t = vertcat(rand_t, temp1);

temp2 = NaN(size(Seasonal_turb_stress,1),1);
Var = [temp2 Var];

clear reg_sig;
for i = 1:2
	
	rand_t(end, :) = [];
	Var(:, 1) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    reg_sig(:,i) = temp3;

end

mask = (abs(windstress_regs) > reg_sig);
mask = double(mask);
mask(mask == 0) = NaN;

%lats_sparse = lats(1:5:261);
%lons_sparse = lons(1:5:1440);

%LA = reshape(double(lats*ones([1 1440])), 261, 1440);
%LO = reshape(double(lons*ones([1 261])), 1440, 261);

%LA = reshape(double(lats_sparse*ones([1 288])), 53, 288);
%LO = reshape(double(lons_sparse*ones([1 53])), 288, 53);

lats_sparse = lats(1:5:181);
lons_sparse = lons(1:5:281);

LA = reshape(double(lats_sparse*ones([1 57])), 37, 57);
LO = reshape(double(lons_sparse*ones([1 37])), 57, 37);

%%
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];

lons_GOA = lons;
lats_GOA = lats;

load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = {'(a)' '(b)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 0;
        hi = [1 3 5 7];
        for k = 1:2
                        
           ilag = ilag+1;
			o = ilag - 1;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        subplot('position', [0.05+(k-1)*0.44 0.67 0.38 0.28])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'off', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 2000;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(windstress_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask(:,k));
plot_data1 = reshape(temp6, 281, 181);
plot_data1_sparse = plot_data1(1:5:281, 1:5:181);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 4.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(WS_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',12)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['WS Dec/Jan Year ' clag],'FontSize',14,'FontWeight','normal');
   
    tightmap
   gridm on 

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
subplot('position',[0.88 0.679 0.18 0.264])
%subplot(5,2,9)
%subplot;

limit = 2000; 
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
h.FontWeight = 'normal';
h.FontSize = 12;
h.LineWidth = 0.05;
%xlabel(h, 'Nm^{-2}s')

%% MLT
ncdisp('/projectnb/pdpanalysis/data/MLD_2020/MLD_2020.nc');

lons = ncread('/projectnb/pdpanalysis/data/MLD_2020/MLD_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/MLD_2020/MLD_2020.nc', 'latitude');
lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);

MLT = ncread('/projectnb/pdpanalysis/data/MLD_2020/MLD_2020.nc', 'mlotst');
MLT(isnan(MLT))=0;

MLD = reshape(MLT, 481*221*12, 28);

MLD = reshape(MLD, 481, 221, 12*28);
MLD_filtered = imboxfilt(MLD, 11); % 5x5 boxcar filter (21) 2.5x2.5 filter (11)
MLD_filtered = reshape(MLD_filtered, 481*221*12,28);
MLD_anom = detrend(MLD_filtered')';

MLD_anom = reshape(MLD_anom, 481*221, 12, 28);
MLD_anom_spr = MLD_anom(:,1:3,:);
Seasonal_MLD = mean(MLD_anom_spr, 2);
Seasonal_MLD = reshape(Seasonal_MLD, 481, 221*1, 28);
Seasonal_MLD = Seasonal_MLD(201:481, 1:181, :);
Seasonal_MLD = reshape(Seasonal_MLD, 281*181, 28);

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);

KI = KI(:, 12:39);
Var = Seasonal_MLD(:, 1:28);

temp1 = NaN(1);
temp2 = NaN(size(Seasonal_MLD,1),1);

KI = [KI temp1];
Var = [temp2 Var];

rows = length(Var(1:end,:));
cols = 2;

lag_Cor_mat = zeros(rows, cols);
lag_Reg_mat = zeros(rows, cols);

for i = 1:2
	
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

MLT_regs = lag_Reg_mat;
Blob_MLT = Seasonal_MLD(:,22:23);

Blob_MLT(isnan(Blob_MLT)) = 0;
MLT_regs(isnan(MLT_regs)) = 0;

for i = 1:2
	a = Blob_MLT(:,i);
	b = MLT_regs(:,i);
	R = corr2(a, b);
	
	MLT_cor_coeff(i) = R;
	
end

KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);
KI = KI(:, 12:39);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/28; % convert to power spectral density

Var = Seasonal_MLD(:, 1:28);
rand_t = randn(28, 1000);

clear temp;

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;


temp1 = NaN(1, 1000);
rand_t = vertcat(rand_t, temp1);

temp2 = NaN(size(Seasonal_MLD,1),1);
Var = [temp2 Var];


clear reg_sig;

for i = 1:2
	
	rand_t(end, :) = [];
	Var(:, 1) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    reg_sig(:,i) = temp3;

end

mask2 = (abs(MLT_regs) > reg_sig);
mask2 = double(mask2);
mask2(mask2 == 0) = NaN;

lats_sparse = lats_GOA(1:5:181);
lons_sparse = lons_GOA(1:5:281);

LA = reshape(double(lats_sparse*ones([1 57])), 37, 57);
LO = reshape(double(lons_sparse*ones([1 37])), 57, 37);


%%
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];


load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = {'(c)' '(d)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 0;
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
		
    limit = 8;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(MLT_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask2(:,k));
plot_data1 = reshape(temp6, 281, 181);
plot_data1_sparse = plot_data1(1:5:281, 1:5:181);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 4.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(MLT_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',12)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['MLT Jan-Mar Year ' clag],'FontSize',14,'FontWeight','normal');
   
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

limit = 8; 
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
%h.Ticks = [-24 -18 -12 -6 0 6 12 18 24];
%xlabel(h, 'Nm^{-2}s')

hold on
%% Tot Heat Flux
lhtfl = ncread('/project/pdpanalysis/nish/data/ERA_5_single_level_energy/ERA5_Single_Level_Sensible_Heat_Flux.nc', 'sshf');
shtfl = ncread('/project/pdpanalysis/nish/data/ERA_5_single_level_energy/ERA5_Single_Level_Latent_Heat_Flux.nc', 'slhf');


lons = ncread('/project/pdpanalysis/nish/data/ERA_5_single_level_energy/ERA5_Single_Level_Sensible_Heat_Flux.nc', 'longitude');
lats = ncread('/project/pdpanalysis/nish/data/ERA_5_single_level_energy/ERA5_Single_Level_Sensible_Heat_Flux.nc', 'latitude');

thflx = lhtfl + shtfl;

%pick Jan-Feb

%here use the mask to get rid of land values.
load '/projectnb/pdpanalysis/data/ERA5_lsm_binary.mat'
lsm_binary1 = lsm_binary(1:241, :);
lsm_binary2 = lsm_binary(1401:1440, :);
lsm_binary = vertcat(lsm_binary2, lsm_binary1);
lsm_binary = lsm_binary(:, 1:181);
lsm_binary = reshape(lsm_binary, 281*181, 1);


thflx_JFM = thflx(201:481,41:221,6:195);
thflx_JF = reshape(thflx_JFM, 281*181, 5, 38);
thflx_JF = thflx_JF(:, 1:2, :);
thflx_JF_mean = mean(thflx_JF, 2);
thflx_JF_mean = reshape(thflx_JF_mean, 281*181*1, 38);
thflx_JF_mean = detrend(thflx_JF_mean')';
thflx_JF_mean = thflx_JF_mean.*lsm_binary;

lons_GOA = lons(201:481,:);
lats_GOA = lats(41:221,:);

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);

KI = KI(:, 1:38);
Var = thflx_JF_mean(:, 1:38);

temp1 = NaN(1);
temp2 = NaN(size(thflx_JF_mean,1),1);

KI = [KI temp1];
Var = [temp2 Var];

rows = length(Var(1:end,:));
cols = 2;

lag_Cor_mat = zeros(rows, cols);
lag_Reg_mat = zeros(rows, cols);

for i = 1:2
	
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

thflx_regs = lag_Reg_mat;
Blob_thflx = thflx_JF_mean(:,34:35);

Blob_thflx(isnan(Blob_thflx)) = 0;
thflx_regs(isnan(thflx_regs)) = 0;



for i = 1:2
	a = Blob_thflx(:,i);
	b = thflx_regs(:,i);
	R = corr2(a, b);
	
	MLT_cor_coeff(i) = R;
	
end

KI = mean(HVD_Cop_LS_EOF_2(1:2, :), 1);
KI = KI(:, 1:38);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/38; % convert to power spectral density

Var = thflx_JF_mean(:, 1:38);
rand_t = randn(38, 1000);

clear temp;

for i = 1:1000
    y = ifft(sqrt(psd).*fft(rand_t(:,i)'));  % Generate a correlated time series
    norm = normalize(y,2);
    temp(:,i) = norm;
end

rand_t = temp;

temp1 = NaN(1, 1000);
rand_t = vertcat(rand_t, temp1);

temp2 = NaN(size(thflx_JF_mean,1),1);
Var = [temp2 Var];

clear reg_sig;

for i = 1:2
	
	rand_t(end, :) = [];
	Var(:, 1) = [];
	
	temp1 = Var*rand_t/size(rand_t,1);
	temp2 = sort(abs(temp1), 2);
	temp3 = temp2(:, 950);
    reg_sig(:,i) = temp3;

end

mask3 = (abs(thflx_regs) > reg_sig);
mask3 = double(mask3);
mask3(mask3 == 0) = NaN;

lats_sparse = lats_GOA(1:5:181);
lons_sparse = lons_GOA(1:5:281);

LA = reshape(double(lats_sparse*ones([1 57])), 37, 57);
LO = reshape(double(lons_sparse*ones([1 37])), 57, 37);

%%
LatLim = [20 65];  %this is where you set the lat/lon limits you want to plot
LonLim = [170 240];


load('PrecipColormaps','precip_cmap')
load('ClimColormaps','cmap_clim')

%panellabel = { '(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' };

tilelab = { '(e)' '(f)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 0;
        hi = [1 3 5 7];
        for k = 1:2
                        
           ilag = ilag+1;
			o = ilag - 1;
            		
		%subplot(5, 2, hi(k))
		%subplot('position',[0.01 .73-(k-1)*.19 0.48 .18])
        %subplot('position', [0.05+(k-1)*0.22 0.01 0.19 0.48])
        subplot('position', [0.05+(k-1)*0.44 0.07 0.38 0.28])
        
	h = worldmap(LatLim,LonLim);
	
	setm(h,'MapProjection','bsam','ParallelLabel', 'on','MeridianLabel', 'on', ...
          'FontSize', 10, ...
          'FontWeight', 'bold', ...
          'LabelFormat', 'none',...
          'MLineLocation', 30,...
          'MLabelLocation',30);
	framem('FLineWidth',0.5);
		
    limit = 1.2e+06;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(thflx_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask3(:,k));
plot_data1 = reshape(temp6, 281, 181);
plot_data1_sparse = plot_data1(1:5:281, 1:5:181);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 4.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(MLT_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',12)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['THF Jan/Feb Year ' clag],'FontSize',14, 'FontWeight','normal');
   
    tightmap
   gridm on 

        end
        
hold on
        
%subplot('position',[0.1 0.01 0.48 0.1])
%subplot('position',[0.02 0.58 0.08 0.3])
%subplot('position',[0.465 0.055 0.088 0.392])
subplot('position',[0.88 0.079 0.18 0.264])
%subplot(5,2,9)
%subplot;

limit = 1.2e+06; 
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
h.FontWeight = 'normal';
h.FontSize = 12;
h.LineWidth = 0.05;
%xlabel(h, 'Nm^{-2}s')

%%
%before exporting the figure make sure to delete latitude labels manually
 %print('Nish_100622_KE_MHW_GOA_EOF_Fig','-dpng', '-r500');
% print('Nish_111722_KE_MHW_GOA_Atmoslink_Fig_v2','-dpng', '-r500');
%print('Nish_021423_KE_MHW_GOA_Atmoslink_Fig_v2','-dpng', '-r500');
print('Nish_051523_KE_MHW_GOA_Atmoslink_Fig_v2','-dpng', '-r500');