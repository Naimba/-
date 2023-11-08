% 水平非均匀基流中波动传播理论——波射线理论2
clear;clc;close all


%% 求出转向纬度解析解
% for ii = 1:1
ii = 1;
    syms phi u_M beta_M  a  Omega k
    k = ii/a;
    % u_M = a*omega;
    u_M = (18.*sin(3.*pi./2.*(1+sin(phi)))+14.*(1-sin(phi).^2))/cos(phi);
    beta_M = 2*Omega*cos(phi)^2/a-cos(phi)/a*diff(1/a/cos(phi)*diff(u_M*cos(phi)^2,phi),phi);
    eq = beta_M./k^2 == u_M;
    eq = subs(eq,{a,Omega},{6.4e6,7.292e-5});
    % eq = subs(eq,{a,Omega,omega},{6.4e6,7.292e-5,(7.292e-5)/30.875});
    alpha = vpasolve(eq,phi,[0 pi/2]);
    alpha = double(alpha*180/pi);
    clear u_M phi eq beta_M Omega k a
    %% 求出截陷纬度解析解
    syms phi u_M
    u_M = (18.*sin(3.*pi./2.*(1+sin(phi)))+14.*(1-sin(phi).^2))/cos(phi);
    phi0 = vpasolve(u_M == 0,phi,[0 pi/2]);
    phi0 = double(phi0*180/pi);
    clear phi u_M
    %% 非均匀基流的第1次计算
    a = 6.4e6;d = 86400;Omega = 7.292e-5;n = 10000;
    k=ii/a;%纬向波数
    h = 0.25*d;%积分步长
    ks = 1;

    ug = @ugroup;vg = @vgroup;
    phi = 30;lambda = 0;lat = zeros(1,n+1);lon = zeros(1,n+1);
    lat(1) = phi;lon(1)=lambda;

    i = 1;
    tic
    while abs(phi)<abs(alpha)
        if abs(phi-phi0)<1e-2
            break
            disp('注意！进入截陷纬度')
        end
        [Z,Y,T] = Runge_Kutta_4_2eq(ug,vg,0,lambda,phi,h,h,k,ks);
        phi = Z(2);lambda = Y(2);
        disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
        lat(i+1) = phi;lon(i+1)=lambda;
        i = i+1;
    end

    lat = lat(1:i-1);lon = lon(1:i-1);

    m_proj('Equidistant Cylindrical','lon',[-90 270],'lat',[-90 90]);
    m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on
    m_grid('box','on','tickdir','in','xtick',-90:30:270,'ytick',-90:30:90);
    m_coast('line','color','k','linewidth',0.2,'linestyle','-');
    set(gcf,'Position',[286.6,237,854,498.4])
    toc
    %% 只考虑纬向基流的第2次计算
    ks = -ks;
    n = 10000;% n为积分次数
    phi = Z(1);
    lambda = interp1([lat(end-1),lat(end)],[lon(end-1),lon(end)],alpha,"linear",'extrap');
    % lambda = Y(1);
    m_plot([Y(1),lambda],[Z(1) Z(1)],'color','r','linewidth',1,'marker','none')
    lat = zeros(1,n+1);lon = zeros(1,n+1);
    lat(1) = phi;lon(1)=lambda;
    i = 1;tic
    while abs(phi)<abs(alpha)
        if abs(phi-phi0)<1e-2
            break
            disp('注意！进入截陷纬度')
        end
        [Z,Y,T] = Runge_Kutta_4_2eq(ug,vg,0,lambda,phi,h,h,k,ks);
        phi = Z(2);lambda = Y(2);
        disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
        lat(i+1) = phi;lon(i+1)=lambda;
        i = i+1;

    end
    lat = lat(1:i);lon = lon(1:i);

    m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on
    toc
% end
%% beta_M
function beta_M = beta_M_(phi)

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

if abs(phi) == 90
    beta_M = 0;
else
    phi0 =phi*pi/180;
    beta_M= (2690273155709801*cos(phi0)^2)/118059162071741130342400000 +...
        (cos(phi0)*((cos(phi0)*(28*cos(phi0)^2 - 28*sin(phi0)^2 + ...
        (81*pi^2*sin((3*pi*(sin(phi0) + 1))/2)*cos(phi0)^2)/2 + 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*sin(phi0))...
        - 2*sin(phi0)*(28*cos(phi0)*sin(phi0) - 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*cos(phi0)) + ...
        cos(phi0)*(18*sin((3*pi*(sin(phi0) + 1))/2) - 14*sin(phi0)^2 + 14))/(6400000*cos(phi0)) + ...
        (sin(phi0)*(sin(phi0)*(18*sin((3*pi*(sin(phi0) + 1))/2) - 14*sin(phi0)^2 + 14) + ...
        cos(phi0)*(28*cos(phi0)*sin(phi0) - 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*cos(phi0))))...
        /(6400000*cos(phi0)^2)))/6400000;
end

end
%% 经向波数
function l = l_num(phi,k)
% phi为纬度，单位为角度
% k为波数

beta_M = beta_M_(phi);
u_M = u_M_fun(phi);
v_M = v_M_fun(phi);

if v_M == 0
    l2 = (beta_M-u_M_fun(phi)*k^2)/u_M_fun(phi);
    if l2 < 0
        l = nan;
        error('需要转向了')
    else
        l = sqrt(l2);
    end
else
    syms aa
    % [~,l,~] = solve3(v_M,u_M*k,v_M*k^2,u_M*k^3-beta_M*k);
    l = vpasolve(v_M*aa^3+u_M*k*aa^2+v_M*k^2*aa+u_M*k^3-beta_M*k==0,aa);% 经向波数
    % l = l(2);
end


end
%% 群速度ug vg
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

beta_M = beta_M_(phi);a = 6.4e6;
l = l_num(phi,k);
if ks == -1
    l = -l;
end


ug = (v_M_fun(phi).*l./k+2*beta_M*k^2/((k.^2+l.^2).^2))/(a*cosd(phi)).*cosd(phi);


end
%-------------------------------------------------
function vg = vgroup(~,~,phi,k,ks)
% phi为纬度，单位为角度
% k为波数

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

beta_M = beta_M_(phi);a = 6.4e6;
l = l_num(phi,k);
if ks == -1
    l = -l;
end

vg = (v_M_fun(phi)+2.*k.*l.*beta_M./((k.^2+l.^2).^2))/(a).*cosd(phi);

% vg = (k*l/(k^2+l^2)*u_M_fun(phi)+(1+(l^2)/(k^2+l^2))*v_M_fun(phi)+ ...
% 2.*k.*l.*beta_M./((k.^2+l.^2).^2))/(a);

end

%% 基流
function u_M = u_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

u_M = (18.*sind(3.*180./2.*(1+sind(phi)))+14.*(1-sind(phi).^2))./cosd(phi);
end

function v_M = v_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

y = sind(phi);y0 = 0.3;
% v_M = 0;
v_M = zeros(1,length(phi));
for i = 1:length(phi)
    if y(i)<=y0 && y(i)>=(-0.5)
        v_M(i) = 3.2.*sind(180.*(y(i)-y0)./(y0+0.5));
    end
    if y(i)<=(0.5) && y(i)>y0
        v_M(i) = 0.8.*sind(180.*(y(i)-y0)./(0.5-y0));
    end
end
v_M = v_M./cosd(phi);
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