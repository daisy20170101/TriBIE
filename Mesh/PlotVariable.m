% read from accuacc.data to plot snapshot of nucleation
% var = load('./lk22_h80_N25_T1/');
% var = load('./lk22_h80_N25_T1/vardep-lk22_h80_N25_T1.dat');
% var = load('./lk20_h90_N25_T2/vardep-lk20_h90_N25_T2.dat');
% var = load('./lk 20_h90_N25_T2_slip/vardep-lk20_h90_N25_T2.dat');
var = load('h90_N25_T2/vardep-h90_N25_T2.dat');
% var = load('eq_h90_N25_T2/NG/vardep-lk20_h90_N25_T2.dat');

depth = var(:,1);
eff = var(:,2);
dc = var(:,3);
pab = var(:,4);
pa = var(:,5);

len=length(pa);

%% calculate accumulative slip on the fault
% num = load('../Patch_num4.txt');
data = load('400km_1km_smooth.gts'); % load meshing - read geometry
nvex = data(1,1); ncell = data(1,3);
vertex = data(2:nvex+1,:);
cell =  data(2+nvex:1+nvex+ncell,:);
 
%%
load('/Users/duoli/Documents/Mexico/2014SSE/Contour/guerrero.mat');
[x,y] = km2lonlat_mex(vertex(:,1)/1000,vertex(:,2)/1000,-55); % rotation back to longitude-latitude

figure;
set(gcf,'position',[100 100 850 750]);
subplot(2,2,1);
hold on;
box on;
set(gca,'linewidth',1.2,'fontsize',12);
% scatter(x(1:1:end,1),y(1:1:end,1),8,100*accslip(1:1:end),'filled');
set(gca,'xlim',[-102.5, -99],'ylim',[16.5,19],'fontsize',12);
% title('Mapview of a-b');
xlabel('longitude');
ylabel('latitude');
colormap jet;

cl = colorbar ;
% cl.Limits = [0 50] ;
caxis([-0.0035 0.0035 ]) ;
cl.Label.String = 'a-b' ;

tr = triangulation(cell(1:len,:),x,y,-vertex(:,3));
trisurf(tr,pab,'edgecolor','none','facealpha',0.9); 

% earthquake
eq = load('/Users/duoli/Documents/Mexico//2014SSE/Contour/2014Eq_USGS.txt');
plot(eq(1,1),eq(1,2),'pr','markerfacecolor','r','markersize',16); % 2014 Papanoa epicenter from USGS
af = load('/Users/duoli/Documents/Mexico//2014SSE/Contour/2014Eq_aftershock.txt');
% plot(af(:,1),af(:,2),'pk','markerfacecolor','y','markersize',16); % 2014 Papanoa epicenter from USGS
% city
city = load('/Users/duoli/Documents/Mexico//2014SSE/Contour/city2.txt');
plot(city(:,1),city(:,2),'sk','markerfacecolor','k','markersize',6); % 2014 Papanoa epicenter from USGS

% set(gca,'Yscale','log','YLim',[1.91e20 1.98e20]);

% trench, coastline and depth contour
trch = load('/Users/duoli/Documents/Mexico//2014SSE/Contour/MAT_trench.txt');
[tr_xr,tr_yr] = lonlat2km_mex(trch(:,1),trch(:,2),0);
plot(trch(:,1),trch(:,2),'-k','linewidth',2.0);
% plot(trch(1:12:end,1),trch(1:12:end,2)+0.015,'<k','markerface','k');
[ncst_xr,ncst_yr] = lonlat2km_mex(ncst(:,1),ncst(:,2),0); %% sphere plattern + rotate 70 degree.
mapshow(ncst(:,1),ncst(:,2),'DisplayType','line','color','k','linewidth',1.2);

depcon = load('/Users/duoli/Documents/Mexico//2014SSE/Contour/contour_latlon2.txt');
mapshow(depcon(:,1),depcon(:,2),'displaytype','line','color','w','linestyle','-.','linewidth',1.5);


% legend('data','2014 Mw7.2','Aftershocks','GPS stations','trench','coastline','isodepth');

%
subplot(2,2,2);
hold on;
box on;
set(gca,'linewidth',1.2,'fontsize',12);
% scatter(x(1:1:end,1),y(1:1:end,1),8,100*accslip(1:1:end),'filled');
set(gca,'xlim',[-102.5, -99],'ylim',[16.5,19],'fontsize',12);
% title('Mapview of effective normal stress');
xlabel('longitude');
ylabel('latitude');
colormap jet;

cl = colorbar ;
% cl.Limits = [0 50] ;
caxis([0 50 ]) ;
cl.Label.String = 'effective normal stress (MPa)' ;

trisurf(tr,eff/10,'edgecolor','none','facealpha',0.9); 

% earthquake
plot(eq(1,1),eq(1,2),'pr','markerfacecolor','r','markersize',16); % 2014 Papanoa epicenter from USGS
% plot(af(:,1),af(:,2),'pk','markerfacecolor','y','markersize',16); % 2014 Papanoa epicenter from USGS
% city
plot(city(:,1),city(:,2),'sk','markerfacecolor','k','markersize',6); % 2014 Papanoa epicenter from USGS

% set(gca,'Yscale','log','YLim',[1.91e20 1.98e20]);

% trench, coastline and depth contour
plot(trch(:,1),trch(:,2),'-k','linewidth',2.0);
% plot(trch(1:12:end,1),trch(1:12:end,2)+0.015,'<k','markerface','k');
mapshow(ncst(:,1),ncst(:,2),'DisplayType','line','color','k','linewidth',1.2);
mapshow(depcon(:,1),depcon(:,2),'displaytype','line','color','w','linestyle','-.','linewidth',1.5);

% legend('data','2014 Mw7.2','Aftershocks','GPS stations','trench','coastline','isodepth');

subplot(2,2,3);
hold on;
box on;
set(gca,'linewidth',1.2,'fontsize',12);
% scatter(x(1:1:end,1),y(1:1:end,1),8,100*accslip(1:1:end),'filled');
set(gca,'xlim',[-102.5, -99],'ylim',[16.5,19],'fontsize',12);
xlabel('longitude');
ylabel('latitude');
colormap jet;

cl = colorbar ;
% cl.Limits = [0 50] ;
caxis([0 200 ]) ;
cl.Label.String = 'dc (mm)' ;

trisurf(tr,dc,'edgecolor','none','facealpha',0.9); 

% plot(af(:,1),af(:,2),'pk','markerfacecolor','y','markersize',16); % 2014 Papanoa epicenter from USGS
% city
plot(city(:,1),city(:,2),'sk','markerfacecolor','k','markersize',6); % 2014 Papanoa epicenter from USGS

% set(gca,'Yscale','log','YLim',[1.91e20 1.98e20]);

% trench, coastline and depth contour
plot(trch(:,1),trch(:,2),'-k','linewidth',2.0);
% plot(trch(1:12:end,1),trch(1:12:end,2)+0.015,'<k','markerface','k');
mapshow(ncst(:,1),ncst(:,2),'DisplayType','line','color','k','linewidth',1.2);
mapshow(depcon(:,1),depcon(:,2),'displaytype','line','color','w','linestyle','-.','linewidth',1.5);

% legend('data','2014 Mw7.2','Aftershocks','GPS stations','trench','coastline','isodepth');

subplot(2,2,4);
hold on;
box on;
set(gca,'linewidth',1.2,'fontsize',12);
% scatter(x(1:1:end,1),y(1:1:end,1),8,100*accslip(1:1:end),'filled');
set(gca,'xlim',[-102.5, -99],'ylim',[16.5,19],'fontsize',12);
xlabel('longitude');
ylabel('latitude');
colormap jet;

cl = colorbar ;
% cl.Limits = [0 50] ;
caxis([0.01 0.05]) ;
cl.Label.String = 'a' ;

trisurf(tr,pa,'edgecolor','none','facealpha',0.9); 
% plot(af(:,1),af(:,2),'pk','markerfacecolor','y','markersize',16); % 2014 Papanoa epicenter from USGS
plot(city(:,1),city(:,2),'sk','markerfacecolor','k','markersize',6); % 2014 Papanoa epicenter from USGS

% set(gca,'Yscale','log','YLim',[1.91e20 1.98e20]);
% trench, coastline and depth contour
plot(trch(:,1),trch(:,2),'-k','linewidth',2.0);
% plot(trch(1:12:end,1),trch(1:12:end,2)+0.015,'<k','markerface','k');
mapshow(ncst(:,1),ncst(:,2),'DisplayType','line','color','k','linewidth',1.2);
mapshow(depcon(:,1),depcon(:,2),'displaytype','line','color','w','linestyle','-.','linewidth',1.5);

% legend('data','2014 Mw7.2','Aftershocks','GPS stations','trench','coastline','isodepth');
