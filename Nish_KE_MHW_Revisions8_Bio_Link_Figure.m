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

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);

KI = KI(:, 12:39);
Var = Seasonal_NO3(:, 1:28);

temp1 = NaN(1);
temp2 = NaN(size(Seasonal_NO3,1),1);

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

no3_regs = lag_Reg_mat;
Blob_no3 = Seasonal_NO3(:,22:23);

Blob_no3(isnan(Blob_no3)) = 0;
no3_regs(isnan(no3_regs)) = 0;

Blob_no3 = reshape(Blob_no3, 281, 181, 2);
Blob_no3 = Blob_no3(121:281, 41:161, :);
Blob_no3 = reshape(Blob_no3, 161*121, 2);

no3_regs = reshape(no3_regs, 281, 181, 2);
no3_regs = no3_regs(121:281, 41:161, :);
no3_regs = reshape(no3_regs, 161*121, 2);

for i = 1:2
	a = Blob_no3(:,i);
	b = no3_regs(:,i);
	R = corr2(a, b);
	
	no3_cor_coeff(i) = R;
	
end

KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);
KI = KI(:, 12:39);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/28; % convert to power spectral density

Var = Seasonal_NO3(:, 1:28);
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

temp2 = NaN(size(Seasonal_NO3,1),1);
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

no3_regs = lag_Reg_mat;

mask2 = (abs(no3_regs) > reg_sig);
mask2 = double(mask2);
mask2(mask2 == 0) = NaN;

lons = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'latitude');

lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);

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
		
    limit = 1.5;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(no3_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask2(:,k));
plot_data1 = reshape(temp6, 281, 181);
plot_data1_sparse = plot_data1(1:5:281, 1:5:181);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 4.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(no3_cor_coeff(k), 1)),'FontWeight','Bold', 'FontSize',12)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['[NO_3^{2-}] Feb/Mar Year ' clag],'FontSize',12,'FontWeight','normal');
   
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

limit = 1.5; 
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
h.Ticks = [-1.5 -1.0 -0.5 0 0.5 1.0 1.5];
h.LineWidth = 0.05;


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

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);

KI = KI(:, 12:39);
Var = Seasonal_po4(:, 1:28);

temp1 = NaN(1);
temp2 = NaN(size(Seasonal_po4,1),1);

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

po4_regs = lag_Reg_mat;
Blob_po4 = Seasonal_po4(:,22:23);

Blob_po4(isnan(Blob_po4)) = 0;
po4_regs(isnan(po4_regs)) = 0;

Blob_po4 = reshape(Blob_po4, 281, 181, 2);
Blob_po4 = Blob_po4(121:281, 41:161, :);
Blob_po4 = reshape(Blob_po4, 161*121, 2);

po4_regs = reshape(po4_regs, 281, 181, 2);
po4_regs = po4_regs(121:281, 41:161, :);
po4_regs = reshape(po4_regs, 161*121, 2);

for i = 1:2
	a = Blob_po4(:,i);
	b = po4_regs(:,i);
	R = corr2(a, b);
	
	po4_cor_coeff(i) = R;
	
end

KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);
KI = KI(:, 12:39);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/28; % convert to power spectral density

Var = Seasonal_po4(:, 1:28);
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

temp2 = NaN(size(Seasonal_po4,1),1);
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

po4_regs = lag_Reg_mat;

mask2 = (abs(po4_regs) > reg_sig);
mask2 = double(mask2);
mask2(mask2 == 0) = NaN;

lons = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'latitude');

lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);

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
		
    limit = 0.1;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(po4_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask2(:,k));
plot_data1 = reshape(temp6, 281, 181);
plot_data1_sparse = plot_data1(1:5:281, 1:5:181);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 4.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(po4_cor_coeff(k), 2)),'FontWeight','Bold', 'FontSize',12)
    textm(24,172, tilelab(k),'FontWeight','Bold', 'FontSize',11)

    clag = int2str(o);
    title(['[PO_4^{3-}] Feb/Mar Year ' clag],'FontSize',12,'FontWeight','normal');
   
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

limit = 0.1; 
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

load HVD_Cop_LS_EOF_2_2020.mat;
HVD_Cop_LS_EOF_2 = reshape(HVD_Cop_LS_EOF_2, 6, 40);
KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);

KI = KI(:, 12:39);
Var = Seasonal_phyc(:, 1:28);

temp1 = NaN(1);
temp2 = NaN(size(Seasonal_phyc,1),1);

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

phyc_regs = lag_Reg_mat;
Blob_phyc = Seasonal_phyc(:,22:23);

Blob_phyc(isnan(Blob_phyc)) = 0;
phyc_regs(isnan(phyc_regs)) = 0;

Blob_phyc = reshape(Blob_phyc, 281, 181, 2);
Blob_phyc = Blob_phyc(121:281, 41:161, :);
Blob_phyc = reshape(Blob_phyc, 161*121, 2);

phyc_regs = reshape(phyc_regs, 281, 181, 2);
phyc_regs = phyc_regs(121:281, 41:161, :);
phyc_regs = reshape(phyc_regs, 161*121, 2);

for i = 1:2
	a = Blob_phyc(:,i);
	b = phyc_regs(:,i);
	R = corr2(a, b);
	
	phyc_cor_coeff(i) = R;
	
end

KI = mean(HVD_Cop_LS_EOF_2(2:3, :), 1);
KI = KI(:, 12:39);

[maxlags,~,~] = size(KI');
[r_KE,lags] = autocorr(KI,maxlags-1); %calculate the autocorrelation of KI

% Compute the power spectral density function
psd = fft(r_KE);
psd = real(psd.*conj(psd))/28; % convert to power spectral density

Var = Seasonal_phyc(:, 1:28);
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

temp2 = NaN(size(Seasonal_phyc,1),1);
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

phyc_regs = lag_Reg_mat;

mask2 = (abs(phyc_regs) > reg_sig);
mask2 = double(mask2);
mask2(mask2 == 0) = NaN;

lons = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'longitude');
lats = ncread('/projectnb/pdpanalysis/data/Nutrients_2020/Nutrients_2020.nc', 'latitude');

lons_GOA = lons(201:481,:);
lats_GOA = lats(1:181,:);

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

tilelab = {'(e)' '(f)'};

%h = figure;
%        set(h, 'Color', 'w');
		
		
		ilag = 0;
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
		
    limit = 0.5;  %here you set the max/min value you want shaded
	contourcmap([-limit:limit*2/255:limit],'gray','colorbar','off')
	cmap_mod = flipud(precip_cmap);
	colormap(cmap_mod)
	caxis([-limit limit]); % set the limits of the colorbar

temp5 = double(phyc_regs(:,k));   %Here you assign the lag-correlation field
plot_data = ...
 reshape(temp5, 281, 181);    %set these for the y/x dimentsion, e.g. 193 is #grids in lat-dir; 256 is #grids in lon-dir

surfm(double(lats_GOA), double(lons_GOA), plot_data','Facecolor', 'interp');  % rlat and rlon needs to be lat and lon values in (ydim,xdim) form;

temp6 = double(mask2(:,k));
plot_data1 = reshape(temp6, 281, 181);
plot_data1_sparse = plot_data1(1:5:281, 1:5:181);

plotm(double(LA), double(LO'), plot_data1_sparse', 'k.', 'MarkerSize', 4.5);

load coast
    geoshow(lat, long, 'DisplayType','line','Color','black', 'LineWidth', 0.5)
    textm(61,-135, num2str(round(phyc_cor_coeff(k), 1)),'FontWeight','Bold', 'FontSize',12)
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

limit = 0.5; 
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
% print('Nish_111722_KE_MHW_GOA_Productivity_Figv2','-dpng', '-r500');
%print('Nish_021423_KE_MHW_GOA_Productivity_Figv2','-dpng', '-r500');

print('Nish_051523_KE_MHW_GOA_Productivity_Figv2','-dpng', '-r500');