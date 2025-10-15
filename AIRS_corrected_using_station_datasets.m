clc; clear; close all

% Part of Sarkar et al. (2025) Section 5 (Usage Notes)

% AIRS Surface δD Bias Correction using Good et al. (2015)-style Regression
% This script performs a bias correction of AIRS δD retrievals
% using collocated in-situ water vapor isotope data. The in-situ water
% vapor isotope datasets are obtained from previous works of Bonne et al. [2014], Steen-Larsen
% et al. [2015], Jacob and Sonntag [1991], Bastrikov et al. [2014], Tremoy et al. [2012], Angert
% et al. [2008], Berkelhammer et al. [2016], Samuels-Crow et al. [2014], Galewsky [2018], Wen
% et al. [2010], Lee et al. [2006], Salamalikis et al. [2015], Gonz´alez et al. [2016], Christensen
% and Knezek [2017], Wei et al. [2019], Dai et al. [2021], Tian et al. [2020], Kurita et al. [2015,
% 2011], Bailey et al. [2015], Leroy-Dos Santos et al. [2020], Vimeux and Risi [2021], Vimeux
% et al. [2024], Wu et al. [2025], Chen et al. [2024], Landais et al. [2024], Bagheri Dastgerdi et al.
% [2020], Leroy-Dos Santos et al. [2023], Bonne et al. [2020], Griffis et al. [2016], Laskar et al.
% [2014]

%  Steps:
%    1.Load AIRS monthly retrieval-space data
%    2.Load TROPESS apriori+DOF data
%    3.Apply DOF,error,and minimum count filters
%    4.Collocate stations in 3D (lon-lat-alt) and aggregate
%    5.Perform Good et al.(2015)-style multivariate regression
%    6.Evaluate via raw,corrected,jackknife,and Monte Carlo stats
%    7.Plot validation diagnostics

min_dof=1;         % DOF screen threshold
nMC=1000;      % Monte Carlo ensemble size
ilevel_mode='z';     % reference mode for level matching


%  Load AIRS monthly data 
fprintf('\n Loading AIRS retrieval-space monthly fields \n');
file='/home/ms377/remote_rdf/Michelle/wisonet/Products/AIRS_RS_monthly_5deg.nc';

ymd=ncread(file,'YYYYMMDD'); ymd=round(ymd./100);  % Convert YYYYMMDD → YYYYMM
lon=ncread(file,'lon');
lat=ncread(file,'lat');
dD=ncread(file,'dD');                     % (lon,lat,level,time)
dD_nobs=ncread(file,'dD_nobs');            % # of observations per grid/month
dd_col_p_error=ncread(file,'dd_col_p_error'); % column precision error
Pressure=ncread(file,'Pressure');          % (lon,lat,level,time)
H2O_vmr=ncread(file,'H2O_vmr');           % (lon,lat,level,time)
level=squeeze(nanmean(Pressure,[1,2,4]));  % representative level structure

fprintf('AIRS data loaded: %d×%d grid,%d levels,%d months.\n',...
    length(lon),length(lat),length(level),size(dD,4));

%  Load TROPESS (Apriori and DOF fields available from 2021-2024) 
fprintf('\n Loading TROPESS apriori+DOF fields \n');
trop_file='/home/ms377/remote_rdf/Michelle/wisonet/Products/AIRS_Monthly/AIRS_TROPESS_all_monthly_5deg.nc';

XA=ncread(trop_file,'XA');                 % (level,lat,lon,time) a-priori isotope ratio
XA=nanmean(XA,4);                          % Time-mean apriori
Rvsmow=3.1152e-4;
XA_dD=((XA./Rvsmow)-1)*1000;         % Convert ratio to δD [permil]


DOF=ncread(trop_file,'DOF');               % (lat,lon,time)
DOF=nanmean(DOF,3)';                       % Time-mean (lon,lat)

fprintf('TROPESS apriori+DOF fields loaded.\n');

%  Quality screening
fprintf('\n Applying uncertainty,count,and DOF masks \n');
[ni,nj,nlev,nt]=size(dD);

% (i) Mask/Remove top 10% error pixels
error90=prctile(dd_col_p_error,90,3);
mask_high_error=dd_col_p_error>error90;
mask_high_error_4d=repmat(mask_high_error,[1 1 nlev 1]);
dD(mask_high_error_4d)=NaN;

% (ii) Mask low-observation months
mask_count=dD_nobs<2;
dD(mask_count)=NaN;

% (iii) Mask low-DOF areas
bad_dof=isnan(DOF)|(DOF<min_dof);
for i=1:ni
    for j=1:nj
        if bad_dof(i,j)
            dD(i,j,:,:)=NaN;
            H2O_vmr(i,j,:,:)=NaN;
        end
    end
end
fprintf('Masking complete.\n');

%  Build pseudo-altitude for level matching 
p0=level(1);     % reference pressure,hPa
T0=288.15;       % surface T,K
g=9.80665;      % gravity,m/s²
R=287.05;       % gas constant,J/kg/K
z=(T0/(g/R))*(1-(level./p0).^(R/g/287.05));  % approximate z based on std atm,m
deg2m=111.32e3;  % approximate
[LON,LAT,Z]=ndgrid(lon,lat,z);

%  Load station dataset 
fprintf('\n Loading station dataset \n');
load('/home/ms377/remote_rdf/Mampi/Data/vapor_station_data.mat'); % loads table T

Tmat=table2array(T(:,2:end));  % exclude name column
T_dD=Tmat(:,4)+8.*(Tmat(:,3));  % δD=8×δ18O+dexcess
T_alt=Tmat(:,5); T_lon=Tmat(:,2); T_lat=Tmat(:,1);
st=Tmat(:,6); et=Tmat(:,7);

% Sort by altitude (for diagnostics)
[~,isort]=sort(T_alt);
T_dD=T_dD(isort); T_alt=T_alt(isort);
T_lon=T_lon(isort); T_lat=T_lat(isort);
st=st(isort); et=et(isort);

fprintf('Station dataset ready: %d total sites.\n',numel(T_lat));

%  3D collocation (lon-lat-alt)+annual aggregation 
fprintf('\n Collocating stations with AIRS grid \n');
AIRS_mean=[]; AIRS_std=[]; H2O_coll=[]; APR_coll=[];
STA_dD=[]; DOF_coll=[]; ALT_coll=[]; lon_coll=[]; lat_coll=[];
lev_idx_keep=[];

all_ym=unique(ymd(:));
for ii=1:numel(T_lat)
    stlat=T_lat(ii); stlon=T_lon(ii); stalt=T_alt(ii);

    % Time range for station record
    stind=find(all_ym==st(ii)); etind=find(all_ym==et(ii));
    if isempty(stind)||isempty(etind)
        continue; 
    end
    ym_span=all_ym(stind:etind);

    % Nearest grid in (lon,lat,z)
    dlat=(LAT-stlat)*deg2m;
    dlon=(LON-stlon)*deg2m.*cosd(stlat);
    dalt=(Z-stalt);
    D=sqrt(dlat.^2+dlon.^2+dalt.^2);
    [~,idx3d]=min(D(:));
    [ilon,ilatg,ilev_z]=ind2sub(size(D),idx3d);

    % Extract monthly subset within station window
    tidx=ismember(ymd,ym_span);
    dD_series=squeeze(dD(ilon,ilatg,ilev_z,tidx));
    h2o_series=squeeze(H2O_vmr(ilon,ilatg,ilev_z,tidx));

    if all(isnan(dD_series))
        continue; 
    end

    % Treat the full station span as one annualized window: get the std dev
    % and mean over the year.
    dD_mean=nanmean(dD_series);
    dD_std=nanstd(dD_series);
    h2o_mean=nanmean(h2o_series);

    apriori_here=XA_dD(ilev_z,ilatg,ilon);
    dof_here=DOF(ilon,ilatg);

    if any(isnan([dD_mean,h2o_mean,apriori_here,dof_here]))
        continue; 
    end

    AIRS_mean(end+1,1)=dD_mean;
    AIRS_std(end+1,1)=dD_std;
    H2O_coll(end+1,1)=h2o_mean;
    APR_coll(end+1,1)=apriori_here;
    STA_dD(end+1,1)=T_dD(ii);
    DOF_coll(end+1,1)=dof_here;
    ALT_coll(end+1,1)=stalt;
    lon_coll(end+1,1)=stlon;
    lat_coll(end+1,1)=stlat;
    lev_idx_keep(end+1,1)=ilev_z;
end
fprintf('Collocation complete: %d matched pairs.\n',numel(STA_dD));

%  Regression setup (Good et al. (2015)-style) 
X=[AIRS_mean,APR_coll,H2O_coll];
y=STA_dD;
valid=all(~isnan([X y]),2);
X=X(valid,:); y=y(valid);
AIRS_std=AIRS_std(valid);
DOF_coll=DOF_coll(valid);
lon_coll=lon_coll(valid); 
lat_coll=lat_coll(valid);
lev_idx_keep=lev_idx_keep(valid);

fprintf('\nMatched pairs after screening: %d\n',numel(y));

%  Raw AIRS vs Station diagnostics 
bias_raw=mean(X(:,1)-y);
rmse_raw=sqrt(mean((X(:,1)-y).^2));
r_raw=corr(X(:,1),y,'rows','complete');

%  Good et al. (2015)-style multivariate regression 
mdl=fitlm(X,y);   % y ~ 1+AIRS_dD+Apriori_dD+H2O
y_hat=predict(mdl,X);

bias_corr=mean(y_hat-y);
rmse_corr=sqrt(mean((y_hat-y).^2));
r_corr=corr(y_hat,y,'rows','complete');

fprintf('\n SUMMARY \n');
fprintf('Raw  bias=%.1f ‰,RMSE=%.1f ‰,r=%.2f\n',bias_raw,rmse_raw,r_raw);
fprintf('Corr bias=%.1f ‰,RMSE=%.1f ‰,r=%.2f\n',bias_corr,rmse_corr,r_corr);

%  Jackknife leave-one-out 
n=length(y);
betas=nan(n,size(X,2)+1);
y_loocv=nan(n,1);
for i=1:n
    idx=setdiff(1:n,i);
    mdli=fitlm(X(idx,:),y(idx));
    betas(i,:)=mdli.Coefficients.Estimate';
    y_loocv(i)=predict(mdli,X(i,:));
end
beta_mean=mean(betas,1);
beta_stdB=std(betas,0,1);

bias_jk=mean(y_loocv-y);
rmse_jk=sqrt(mean((y_loocv-y).^2));
r_jk=corr(y_loocv,y,'rows','complete');

fprintf('\n Jackknife (coeffs & LOO predictions) \n');
fprintf('β₀=%.2f ± %.2f\n',beta_mean(1),beta_stdB(1));
fprintf('β_AIRS=%.2f ± %.2f\n',beta_mean(2),beta_stdB(2));
fprintf('β_aprior=%.2f ± %.2f\n',beta_mean(3),beta_stdB(3));
fprintf('β_H₂O=%.2f ± %.2f\n',beta_mean(4),beta_stdB(4));
fprintf('LOO Bias=%.2f ‰,RMSE=%.2f ‰,r=%.2f\n',bias_jk,rmse_jk,r_jk);

%  Monte Carlo ensemble (around β mean +/- std) 
y_mc=nan(n,nMC);
for k=1:nMC
    b=beta_mean+beta_stdB.*randn(size(beta_mean));
    y_mc(:,k)=b(1)+X*b(2:end)';
end
y_mc_mean=mean(y_mc,2);

bias_mc=mean(y_mc_mean-y);
rmse_mc=sqrt(mean((y_mc_mean-y).^2));
r_mc=corr(y_mc_mean,y,'rows','complete');

fprintf('\n Monte Carlo (n=%d) \n',nMC);
fprintf('Bias=%.2f ‰,RMSE=%.2f ‰,r=%.2f\n',bias_mc,rmse_mc,r_mc);

%  Validation Plots
resid_corr=y-y_hat;

% Diverging red–white–blue colormap
ncolors=256;
base=[0 0 1; 1 1 1; 1 0 0];
x=linspace(0,1,ncolors);
cmap_rwb=zeros(ncolors,3);
for i=1:3
    cmap_rwb(:,i)=interp1([0 0.5 1],base(:,i),x);
end


% (a) Residual map
figure('Color','w'); hold on; box on; grid on;
scatter(lon_coll,lat_coll,90,resid_corr,'filled','MarkerEdgeColor','k');
load coastlines.mat; plot(coastlon,coastlat,'k-','LineWidth',0.8);
cmax=max(abs(resid_corr)); caxis([-cmax cmax]);
colormap(cmap_rwb);
cb=colorbar; ylabel(cb,'Station – Corrected AIRS δD (‰)');
xlabel('Longitude'); ylabel('Latitude');
set(gca,'FontSize',13,'FontWeight','bold');
title('Residuals (Station − Corrected AIRS)');

% (b) Corrected scatter vs Station (color=residual)
figure('Color','w'); hold on; box on; grid on;
for i=1:numel(y)
    errorbar(y_hat(i),y(i),0,0,AIRS_std(i),AIRS_std(i),...
        'o','Color',[0 0 0 0.3],'CapSize',6,'LineWidth',1);
end
scatter(y_hat,y,80,resid_corr,'filled','MarkerEdgeColor','k');
plot([-250 0],[-250 0],'k--','LineWidth',1.2);
cmax2=max(abs(resid_corr)); caxis([-cmax2 cmax2]);
colormap(cmap_rwb);
cb=colorbar; ylabel(cb,'Residual (‰)');
xlabel('Bias-corrected AIRS δD (‰)'); ylabel('Station δD (‰)');
xlim([-250 0]); ylim([-250 0]);
set(gca,'FontSize',13,'FontWeight','bold');
title(sprintf('Bias-corrected vs Station (r=%.2f)',r_corr));

% (c) Corrected scatter vs Station (color=DOF)
figure('Color','w'); hold on; box on; grid on;
for i=1:numel(y)
    errorbar(y_hat(i),y(i),0,0,AIRS_std(i),AIRS_std(i),...
        'o','Color',[0 0 0 0.3],'CapSize',6,'LineWidth',1);
end
scatter(y_hat,y,80,DOF_coll,'filled','MarkerEdgeColor','k');
plot([-250 0],[-250 0],'k--','LineWidth',1.2);
colormap(turbo);
cb=colorbar; ylabel(cb,'DOF (trace)');
xlabel('Bias-corrected AIRS δD (‰)'); ylabel('Station δD (‰)');
xlim([-250 0]); ylim([-250 0]);
set(gca,'FontSize',13,'FontWeight','bold');
title('Bias-corrected vs Station (color=DOF)');
