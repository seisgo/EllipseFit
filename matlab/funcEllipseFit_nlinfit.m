% Ellipse Fitting using nlinfit function and adopting conic section as
% nonlinear regression model function.
% This function can be extended to other conic section fitting, like
% circle, parabola and hyperbola.
function [modelF,fitCoef]=funcEllipseFit_nlinfit(inMtrx)
% Ellipse Function
modelF=@(fitCoef,inMtrx)fitCoef(1)*inMtrx(:,1).^2 ...,
    +fitCoef(2)*inMtrx(:,1).*inMtrx(:,2) ...,
    +fitCoef(3)*inMtrx(:,2).^2 ...,
    +fitCoef(4)*inMtrx(:,1) ...,
    +fitCoef(5)*inMtrx(:,2) ...,
    +fitCoef(6);
% Initial Coefficients
coef=[1 1 1 1 1 1];
% Do Nonlinear Regression
fitCoef=nlinfit(inMtrx,zeros(size(inMtrx,1),1),modelF,coef);
