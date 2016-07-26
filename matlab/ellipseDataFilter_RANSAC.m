function filterData=ellipseDataFilter_RANSAC(data)
% Do ellipse scatter data filtering for ellipse fitting by RANSAC method.
% Author: Zhenyu Yuan
% Date: 2016/7/26
% Ref:  http://www.cnblogs.com/yingying0907/archive/2012/10/13/2722149.html
%       Extract RANSAC filtering in ellipsefit.m and make some modification
%%  参数初始化
nSampLen = 3;               %设定模型所依据的点数
nDataLen = size(data, 1);   %数据长度
nIter = 50;                 %最大循环次数
dThreshold = 2;             %残差阈值
nMaxInlyerCount=-1;         %点数下限
A=zeros([2 1]);
B=zeros([2 1]);
P=zeros([2 1]);
%%  主循环
for i = 1:nIter 
    ind = ceil(nDataLen .* rand(1, nSampLen)); %抽样，选取nSampLen个不同的点
    %%  建立模型，存储建模需要的坐标点,焦点和过椭圆的一个点
    %椭圆定义方程：到两定点之间距离和为常数
    A(:,1)=data(ind(1),:);    %焦点
    B(:,1)=data(ind(2),:);    %焦点
    P(:,1)=data(ind(3),:);    %椭圆上一点
    DIST=sqrt((P(1,1)-A(1,1)).^2+(P(2,1)-A(2,1)).^2)+sqrt((P(1,1)-B(1,1)).^2+(P(2,1)-B(2,1)).^2);
    xx=[];
    nCurInlyerCount=0;        %初始化点数为0个
    %%  是否符合模型？
    for k=1:nDataLen
        %         CurModel=[A(1,1)   A(2,1)  B(1,1)  B(2,1)  DIST ];
        pdist=sqrt((data(k,1)-A(1,1)).^2+(data(k,2)-A(2,1)).^2)+sqrt((data(k,1)-B(1,1)).^2+(data(k,2)-B(2,1)).^2);
        CurMask =(abs(DIST-pdist)< dThreshold);     %到直线距离小于阈值的点符合模型,标记为1
        nCurInlyerCount =nCurInlyerCount+CurMask;             %计算符合椭圆模型的点的个数
        if(CurMask==1)
            xx =[xx;data(k,:)];
        end
    end
    %% 选取最佳模型
    if nCurInlyerCount > nMaxInlyerCount   %符合模型的点数最多的模型即为最佳模型
        nMaxInlyerCount = nCurInlyerCount;
        %             Ellipse_mask = CurMask;
        %              Ellipse_model = CurModel;
        %              Ellipse_points = [A B P];
        filterData =xx;
    end
end