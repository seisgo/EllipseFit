% Prepare scatter data for ellipse fitting
%
clear; clc; close all;
%% Fit-failed data from Prestack Analysis in EPoffice
x=[3177.48 -1330.11
   1607.04 -459.156
    371.41  136.436
   1583.09  434.007
   3067.27  115.600
   2370.14 -86.1865
   3059.32 -1873.05
   5077.39 -3665.81
   3824.05  2657.03
   2147.30  1171.25
         0        0
         0        0
  -3177.48  1330.11
  -1607.04  459.156
   -371.41 -136.436
  -1583.09 -434.007
  -3067.27 -115.600
  -2370.14  86.1865
  -3059.32  1873.05
  -5077.39  3665.81
  -3824.05 -2657.03
  -2147.30 -1171.25
         0        0
         0        0
];
x=x/max(max(x));
save hyperbolaData
%
%% Fit-success data from Prestack Analysis in EPoffice
x=[
    3971.31 -1662.41
    758.616 -216.747
    88.4086  32.4766
    1664.47  456.316
    4625.38  174.323
    4261.29 -154.955
    1635.06 -1001.05
    3935.46 -2841.35
    4174.71  2900.68
          0        0
          0        0
          0        0
   -3971.31  1662.41
   -758.616  216.747
   -88.4086 -32.4766
   -1664.47 -456.316
   -4625.38 -174.323
   -4261.29  154.955
   -1635.06  1001.05
   -3935.46  2841.35
   -4174.71 -2900.68
          0        0
          0        0
          0        0
    ];
x=x/max(max(x));
save hyperEllipData
%
%% Ideal Simple Test Data
x=[1.7729 1.9228
   1.7338 1.9072
   2.0539 1.6137
   2.0656 1.6412
   1.8611 1.48765
   1.9005 1.4971
   2.0732 1.6546
   1.8338 1.9405
   1.9375 1.5104
   1.6878 1.5177
   1.7031 1.5097
   1.9577 1.5201
   1.9872 1.5437
   2.0341 1.5805
   2.0723 1.6546
   2.0681 1.8284
   2.0557 1.8483
   2.0491 1.5651
];
x=x/max(max(x));
save ellipseData

%%  Random Noisy Ellipse Data
% Ref:  http://www.cnblogs.com/yingying0907/archive/2012/10/13/2722149.html
%       Extract noisy ellipse data generation in ellipsefit.m and transpose
%       2*n matrix to n*2 matrix
% 参数初始化
g_NumOfPoints = 500;   % 点数
g_ErrPointPart = 0.4;  % 噪声
g_NormDistrVar = 3;    % 标准偏差
a=10;b=20;             %长轴短轴
angle=60;              %倾斜角
%% 椭圆生成
beta = angle * (pi / 180);
alpha = linspace(0, 360, g_NumOfPoints) .* (pi / 180); 
X = (a * cos(alpha) * cos(beta)- b * sin(alpha) * sin(beta) )+wgn(1,length(alpha),g_NormDistrVar^2,'linear');    
Y = (a * cos(alpha) * sin(beta)+ b * sin(alpha) * cos(beta) )+wgn(1,length(alpha),g_NormDistrVar^2,'linear');
x=[X' Y'];
% x=x/max(max(x));
clearvars -except x;
save noisyEllipData