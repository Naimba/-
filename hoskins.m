% clear;clc;close all
% 复现Hoskins的结果

%% 第一次积分，启动！
nexttile
m_proj('Equidistant Cylindrical','lon',[0 360],'lat',[-90 90]);
m_grid('box','on','tickdir','in','xtick',0:30:360,'ytick',-90:30:90);
m_coast('line','color','k','linewidth',0.2,'linestyle','-');hold on
for ii = 1:5
a = 6400000;
d = 86400;h = 0.1*d;ks = 1;
u_M = 6.4e6/30.875*7.292e-5;k=ii/a;

alpha = acosd(a*k*sqrt(1/63.75));
% l = sqrt((beta_M-u_M*k^2)/u_M);
% vg =2.*k.*l.*u_M^2./beta_M/a;
% ug = (2*u_M^2*k^2/beta_M)/(a*cosd(alpha));

n = 1230;% n为积分次数

ug = @ugroup;vg = @vgroup;
phi = 0;lambda = 0;lat = zeros(1,n+1);lon = zeros(1,n+1);
lat(1) = phi;lon(1)=lambda;

i = 1;
% tic
while abs(phi)<abs(alpha)
    [Z,Y,T] = Runge_Kutta_4_2eq(ug,vg,0,lambda,phi,h,h,k,ks);
    phi = Z(2);lambda = Y(2);
    % disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
    lat(i+1) = phi;lon(i+1)=lambda;
    i = i+1;
end

lat = lat(1:i-1);lon = lon(1:i-1);

m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on

% set(gcf,'Position',[286.6,237,854,498.4])
% toc

hold on
%% 转向，启动！
ks = -ks;
n = 2440;% n为积分次数
lambda = interp1([lat(end-1),lat(end)],[lon(end-1),lon(end)],alpha,"linear",'extrap');

phi = Z(1);
% lambda = Y(1);
lambda = 2*lambda-Y(1);
m_text(90,phi+3,num2str(ii,'%d'))
m_plot([Y(1),lambda],[Z(1) Z(1)],'color','r','linewidth',1,'marker','none')
lat = zeros(1,n+1);lon = zeros(1,n+1);
lat(1) = phi;lon(1)=lambda;

i = 1;
% tic
while abs(phi)<abs(alpha)
    [Z,Y,T] = Runge_Kutta_4_2eq(ug,vg,0,lambda,phi,h,h,k,ks);
    phi = Z(2);lambda = Y(2);
    % disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
    lat(i+1) = phi;lon(i+1)=lambda;
    i = i+1;
end
lat = lat(1:i-1);lon = lon(1:i-1);

m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on
% toc

%% 二次转向，启动！
ks = -ks;
n = 1000;% n为积分次数
% disp(['共需积分 ',num2str(n*h/d,'%4.2f'),'d'])
lambda = interp1([lat(end-1),lat(end)],[lon(end-1),lon(end)],-alpha,"linear",'extrap');
lambda = 2*lambda-Y(1);
m_text(270,phi-3,num2str(ii,'%d'))
m_plot([Y(1),lambda],[Z(1) Z(1)],'color','r','linewidth',1,'marker','none')
phi = Z(1);
% lambda = Y(1);

lat = zeros(1,n+1);lon = zeros(1,n+1);
lat(1) = phi;lon(1)=lambda;
i = 1;
% tic
while abs(phi)<abs(alpha)
    [Z,Y,T] = Runge_Kutta_4_2eq(ug,vg,0,lambda,phi,h,h,k,ks);
    phi = Z(2);lambda = Y(2);
    % disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
    lat(i+1) = phi;lon(i+1)=lambda;
    i = i+1;
end
lat = lat(1:i-1);lon = lon(1:i-1);

h1 = m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on
% toc
%% 理论上的波射线
hold on
lambda0 = 0;
lambda = 0:360;
phi = atand(tand(alpha).*sind(lambda-lambda0));
% m_proj('Equidistant Cylindrical','lon',[0 360],'lat',[-90 90]);
h2 = m_plot(lambda,phi,'color','b','linewidth',1,'linestyle',':');
% m_grid('box','on','tickdir','in','xtick',0:30:360,'ytick',-90:30:90);
% m_coast('line','color','k','linewidth',0.2,'linestyle','-');
% set(gcf,'Position',[286.6,237,854,498.4])
end
m_text(345,80,'(a)')
% legend([h1 h2],{'计算解','理论解'})
% print(gcf,['F:\学习\毕业论文\复现李艳杰\hoskins1-5波'],'-dpng','-r400');

%% 经向波数
function l = l_num(phi,k)
% phi为纬度，单位为角度
% k为波数

beta_M = beta_M_(phi);a = 6.4e6;
% u_M0 = beta_M/k/k-(6.4e6)*(7.292e-5);%u_M0为转向临界风速
u_M = 6.4e6/30.875*7.292e-5;

l2 = (beta_M-u_M*k^2)/u_M;
if l2 < 0
    l = nan;
    error('需要转向了')
else
    l = sqrt(l2);
end


end
%% 群速度
function ug = ugroup(~,~,phi,k,ks)
% ugroup(t,lambda,phi,k,l)
% phi为纬度，单位为角度
% k为波数

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

if ks ~= -1 && ks ~= 1
    error('ks的取值只有-1或1')
end

l = l_num(phi,k);
if ks == -1
    l = -l;
end

beta_M = beta_M_(phi);a = 6.4e6;u_M = 6.4e6/30.875*7.292e-5;

ug = (2.*u_M.^2.*k.^2/beta_M)/(a.*cosd(phi)).*cosd(phi);

end
%-------------------------------------------------
function vg = vgroup(~,~,phi,k,ks)
% phi为纬度，单位为角度
% k为波数

beta_M = beta_M_(phi);

if ks ~= -1 && ks ~= 1
    error('ks的取值只有-1或1')
end

l = l_num(phi,k);
if ks == -1
    l = -l;
end
u_M = 6.4e6/30.875*7.292e-5;a = 6.4e6;

vg = (2.*k.*l.*u_M.^2)./beta_M./(a).*cosd(phi);


end

%% beta_M  
function beta_M = beta_M_(phi)

a = 6400000;
beta_M = 2*cosd(phi)^2/a*(7.292e-5+1/30.875*7.292e-5);
% beta_M = 2e-11;
end

%% 有两个方程的四阶Runge-Kutta方法
function [Z,Y,T] = Runge_Kutta_4_2eq(f,g,t0,y0,z0,h,t,k,ks)
% 输入参数说明：
% y'=f(t,y,z)；z'=g(t,y,z)；y = y(t)；z = z(t)；但是在这里ugroup(t,lambda,phi,k,l)
% t0为初始t值
% y0为初始t值对于的y值，y0=y(y0)；z0为初始t值对于的z值，z0=z(z0)
% h为积分步长
% t为积分取值
% k为波数
% ks为l取正负的选项，1为正，-1为负
% --------------------------------------------------------------------------------------------
% 输出参数说明：
% X为x0到x区间上等间距h的数组
% Y(i)=y(T(i))

if ks ~= -1 && ks ~= 1
    error('ks的取值只有-1或1')
end

n = floor((t-t0)/h);%求步数
T = zeros(1,n+1);
Y = zeros(1,n+1);
Z = zeros(1,n+1);
T(1) = t0;Y(1) = y0;Z(1) = z0;

for i = 1:n
    k1 = h*f(T(i),Y(i),Z(i),k,ks);
    l1 = h*g(T(i),Y(i),Z(i),k,ks);
    k2 = h*f(T(i)+h/2,Y(i)+k1/2,Z(i)+l1/2,k,ks);
    l2 = h*g(T(i)+h/2,Y(i)+k1/2,Z(i)+l1/2,k,ks);
    k3 = h*f(T(i)+h/2,Y(i)+k2/2,Z(i)+l2/2,k,ks);
    l3 = h*g(T(i)+h/2,Y(i)+k2/2,Z(i)+l2/2,k,ks);
    k4 = h*f(T(i)+h,Y(i)+k3,Z(i)+l3,k,ks);
    l4 = h*g(T(i)+h,Y(i)+k3,Z(i)+l3,k,ks);
    Y(i+1) = Y(i)+(k1+2*k2+2*k3+k4)/6*180/pi;
    Z(i+1) = Z(i)+(l1+2*l2+2*l3+l4)/6*180/pi;
    T(i+1) = T(i)+h;
    %按照龙格库塔方法进行数值求解
end

end