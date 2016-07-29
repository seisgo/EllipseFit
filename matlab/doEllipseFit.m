% Perform ellipse fitting by calling customized functions.
clear; clc; close all;
RBrownOK = 1;   % In some case, Richard Browm method doesn't work
%% Input Raw Data
%
% Fit-failed data from Prestack Analysis in EPoffice
load hyperbolaData;
RBrownOK = 0;
%{
% Fit-success data from Prestack Analysis in EPoffice
load hyperEllipData;
RBrownOK = 0;
%
% Simple Test Data
load ellipseData;
%
% Random Noisy Ellipse Data
load noisyEllipData
%}

figure;
set(gcf,'Units','centimeters','Position',[10 10 18 13]);
% plot raw data
plot(x(:,1),x(:,2),'b.','Linewidth',1.5);
hold on;

%% Do RANSAC filter of input raw data;
filter=1;   % Flag of whether do RANSAC filtering
if filter
    tmp=ellipseDataFilter_RANSAC(x);
    x=tmp;
    % plot filtered data
    plot(x(:,1),x(:,2),'yo','Linewidth',1.5);
end

%% Do Ellipse Fitting by calling funcEllipseFit_nlinfit
[F,p]=funcEllipseFit_nlinfit(x);
% plot fitted curve
xmin=min(x(:,1));
xmax=max(x(:,1));
ymin=min(x(:,2));
ymax=max(x(:,2));
shift=0.4*(xmax-xmin);
xlmt1=xmin-shift;
xlmt2=xmax+shift;
shift=0.2*(ymax-ymin);
ylmt1=ymin-shift;
ylmt2=ymax+shift;
axis equal;
ezplot(@(x,y)F(p,[x,y]),[xlmt1,xlmt2,ylmt1,ylmt2]);
title('Ellipse fitting by different methods','FontWeight','bold','fontsize',14);
ylabel('Y','FontWeight','bold','fontsize',12)
xlabel('X','FontWeight','bold','fontsize',12)
set(gca,'fontweight','bold','fontsize',10);

%% Do Ellipse Fitting by calling funcEllipseFit_OGal
ellipse_t = funcEllipseFit_OGal(x(:,1), x(:,2));
if isempty(ellipse_t.status)==0
    warning('Ellipse fitting failed!');
    plot(x(:,1),x(:,2),'r.');
else
% Default azimuth should be in degree unit in function 'ellipse1'
phi=-ellipse_t.phi*180/pi;
[lat,lon] = ellipse1(ellipse_t.X0_in,ellipse_t.Y0_in,[ellipse_t.long_axis/2,...,
    axes2ecc(ellipse_t.long_axis/2,ellipse_t.short_axis/2)],phi);
plot(lat,lon,'c');
end

%% Do Ellipse Fitting by calling funcEllipseFit_BFisher, a direct method

ellipse_1 = funcEllipseFit_BFisher(x(:,1), x(:,2));
maxR=max(ellipse_1(3),ellipse_1(4));
minR=min(ellipse_1(3),ellipse_1(4));
ecc=axes2ecc(maxR,minR);
% Default azimuth should be in degree unit in function 'ellipse1'
if ellipse_1(3)>ellipse_1(4)
    phi=ellipse_1(5)*180/pi;
else
    phi=90+ellipse_1(5)*180/pi;
end
[lat,lon] = ellipse1(ellipse_1(1),ellipse_1(2),[maxR,ecc],phi);
plot(lat,lon,'r*');

%% Do Ellipse Fitting by calling funcEllipseFit_direct, a direct method
para = funcEllipseFit_direct(x);
ezplot(@(x,y)F(para,[x,y]),[xlmt1,xlmt2,ylmt1,ylmt2]);
title('Ellipse fitting by different methods','FontWeight','bold','fontsize',14);

%% Do Ellipse Fitting by calling funcEllipseFit_RBrown
if RBrownOK
    [zg, ag, bg, alphag] = funcEllipseFit_RBrown(x');
    [zb, ab, bb, alphab] = funcEllipseFit_RBrown(x','linear');
    [zt, at, bt, alphat] = funcEllipseFit_RBrown(x','linear', 'constraint', 'trace');
    plotellipse(zg, ag, bg, alphag, 'k')
    plotellipse(zb, ab, bb, alphab, 'm--')
    plotellipse(zt, at, bt, alphat, 'g--')
end
hold off;

if filter
    legend('RawData','FilterData','Nlinfit','OhadGal','BFisher','Direct',...,
        'Geometric','Bookstein','Trace');
else
    legend('RawData','Nlinfit','OhadGal','BFisher','Direct','Geometric',...,
        'Bookstein','Trace');
end
