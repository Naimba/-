% clear;clc;close all
% 使用dl/dt的方法计算非均匀基流
%% dl的构造
% syms y phi k a Omega l
% y0 = 0.3;
% y = sin(phi);
% u_M = (18*sin(3*pi/2*(1+y))+14*(1-y.^2))/cos(phi);
% beta_M = 2*Omega*cos(phi)^2/a-cos(phi)/a*diff(1/a/cos(phi)*diff(u_M*cos(phi)^2,phi),phi);
% % -0.5<=y<=y0
% v_M = 3.2.*sin(pi.*(y-y0)./(0.5+y0))/cos(phi);
% dl = -k*cos(phi)/a*diff(u_M,phi)-l*cos(phi)/a*diff(v_M,phi)+cos(phi)/a*diff(beta_M,phi)*k/(k^2+l^2);
% dl = simplify(dl);
% % y0<=y<=0.5
% v_M = 0.8.*sin(pi.*(y-y0)./(0.5-y0))/cos(phi);
% dl = -k*cos(phi)/a*diff(u_M,phi)-l*cos(phi)/a*diff(v_M,phi)+cos(phi)/a*diff(beta_M,phi)*k/(k^2+l^2);
% dl = simplify(dl);
% abs(y)>30
% v_M = 0;
% dl = -k*cos(phi)/a*diff(u_M,phi)-l*cos(phi)/a*diff(v_M,phi)+cos(phi)/a*diff(beta_M,phi)*k/(k^2+l^2);
% dl = simplify(dl);
%% 非均匀基流的计算
a = 6.4e6;d = 86400;Omega = 7.292e-5;n = 500;
h = 0.1*d;%积分步长
ug = @ugroup;vg = @vgroup;dl = @dl_fun;

% tic
% tiledlayout(1,1)
% nexttile
m_proj('Equidistant Cylindrical','lon',[-90 270],'lat',[-90 90]);
m_grid('box','on','tickdir','in','xtick',-90:30:270,'ytick',-90:30:90);
m_coast('line','color','k','linewidth',0.2,'linestyle','-');hold on
for ks = [1 -1]
    for ii = 1:5
        k =ii/a;

        phi = -30;lambda = 0;lat = zeros(1,n+1);lon = zeros(1,n+1);
        lat(1) = phi;lon(1)=lambda;
        beta_M = beta_M_(phi);
        l = ks*sqrt((beta_M-u_M_fun(phi)*k^2)/u_M_fun(phi));

        i = 1;
        while i<=n
            [Lambda,Phi,L,T] = Runge_Kutta_4_2eq(ug,vg,dl,0,lambda,phi,l,h,h,k);
            phi = Phi(2);
            lambda = Lambda(2);l = L(2);
            % disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
            lat(i+1) = phi;lon(i+1)=lambda;
            i = i+1;

        end

        lat = lat(1:i-1);lon = lon(1:i-1);

        m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on
        m_scatter(lon(1:2*d/h:end),lat(1:2*d/h:end),8,'k','filled')

        % set(gcf,'Position',[286.6,237,854,498.4])
    end
end
% toc
%% 图题等等
% title('非均匀基流的定常波射线及传播距离1','fontsize',12,'fontname','Yahei','Interpreter','tex','FontAngle','italic')
% m_text(35,78,'1');m_text(45,70,'2');m_text(60,62,'3');m_text(55,50,'4');m_text(35,40,'5');
% m_text(-8,27,'1');m_text(15,27,'5');
m_text(35,-78,'1');m_text(45,-70,'2');m_text(60,-62,'3');m_text(55,-50,'4');m_text(35,-40,'5');
m_text(-8,-27,'1');m_text(15,-27,'5');
% print(gcf,'F:\学习\毕业论文\复现李艳杰\非均匀基流的定常波射线及传播距离1','-dpng','-r400');
% close
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
% v_M =0;
end

%% 群速度ug vg
function ug = ugroup(~,~,phi,l,k)
% ugroup(t,lambda,phi,l,k)
% phi为纬度，单位为角度
% k为波数

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

beta_M = beta_M_(phi);a = 6.4e6;



ug = (-v_M_fun(phi).*l./k+2*beta_M*k^2/((k.^2+l.^2).^2))/(a*cosd(phi)).*cosd(phi);


end
%-------------------------------------------------
function vg = vgroup(~,~,phi,l,k)
% phi为纬度，单位为角度
% k为波数

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

beta_M = beta_M_(phi);a = 6.4e6;


vg = (v_M_fun(phi)+2.*k.*l.*beta_M./((k.^2+l.^2).^2))/(a).*cosd(phi);

end
%% dl/dt
function dl = dl_fun(~,~,phi,l,k)
if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

y0 = 0.3;y = sind(phi);a = 6.4e6;Omega = 7.292e-5;
phi = phi*pi/180;
if y<=y0 && y>=(-0.5)
    dl = -(360*k*sin(phi)*cos((3*pi*sin(phi))/2) - 840*k*cos(phi)^2*sin(phi) + 5040*k*cos(phi)^4*sin(phi) - 280*a^2*k^3*cos(phi)^2*sin(phi) + 540*k*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 - 2160*k*pi*sin((3*pi*sin(phi))/2)*cos(phi)^4 - 1215*k*pi^3*sin((3*pi*sin(phi))/2)*cos(phi)^6 - 360*a^2*k^3*sin(phi)*cos((3*pi*sin(phi))/2) - 360*a^2*k*l^2*sin(phi)*cos((3*pi*sin(phi))/2) + 32*a^2*l^3*sin((5*pi*sin(phi))/4)*sin(phi)*(2 - 2^(1/2))^(1/2) - 280*a^2*k*l^2*cos(phi)^2*sin(phi) + 540*a^2*k^3*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 - 4050*k*pi^2*cos(phi)^4*sin(phi)*cos((3*pi*sin(phi))/2) + 80*Omega*a*k*cos(phi)^3*sin(phi) - 32*a^2*l^3*sin(phi)*cos((5*pi*sin(phi))/4)*(2^(1/2) + 2)^(1/2) - 32*a^2*k^2*l*sin(phi)*cos((5*pi*sin(phi))/4)*(2^(1/2) + 2)^(1/2) + 40*a^2*l^3*pi*sin((5*pi*sin(phi))/4)*cos(phi)^2*(2^(1/2) + 2)^(1/2) + 32*a^2*k^2*l*sin((5*pi*sin(phi))/4)*sin(phi)*(2 - 2^(1/2))^(1/2) + 540*a^2*k*l^2*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 + 40*a^2*l^3*pi*cos(phi)^2*cos((5*pi*sin(phi))/4)*(2 - 2^(1/2))^(1/2) + 40*a^2*k^2*l*pi*cos(phi)^2*cos((5*pi*sin(phi))/4)*(2 - 2^(1/2))^(1/2) + 40*a^2*k^2*l*pi*sin((5*pi*sin(phi))/4)*cos(phi)^2*(2^(1/2) + 2)^(1/2))/(20*a^3*cos(phi)*(k^2 + l^2));
elseif y<=(0.5) && y>y0
    dl = (840*k*cos(phi)^2*sin(phi) - 360*k*sin(phi)*cos((3*pi*sin(phi))/2) - 5040*k*cos(phi)^4*sin(phi) + 280*a^2*k^3*cos(phi)^2*sin(phi) - 540*k*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 + 2160*k*pi*sin((3*pi*sin(phi))/2)*cos(phi)^4 + 1215*k*pi^3*sin((3*pi*sin(phi))/2)*cos(phi)^6 + 360*a^2*k^3*sin(phi)*cos((3*pi*sin(phi))/2) - 16*a^2*l^3*sin(phi)*cos(5*pi*sin(phi)) + 360*a^2*k*l^2*sin(phi)*cos((3*pi*sin(phi))/2) - 16*a^2*k^2*l*sin(phi)*cos(5*pi*sin(phi)) + 280*a^2*k*l^2*cos(phi)^2*sin(phi) - 540*a^2*k^3*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 + 80*a^2*l^3*pi*sin(5*pi*sin(phi))*cos(phi)^2 + 4050*k*pi^2*cos(phi)^4*sin(phi)*cos((3*pi*sin(phi))/2) - 80*Omega*a*k*cos(phi)^3*sin(phi) - 540*a^2*k*l^2*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 + 80*a^2*k^2*l*pi*sin(5*pi*sin(phi))*cos(phi)^2)/(20*a^3*cos(phi)*(k^2 + l^2));
else
    dl = (168*k*cos(phi)^2*sin(phi) - 72*k*sin(phi)*cos((3*pi*sin(phi))/2) - 1008*k*cos(phi)^4*sin(phi) +...
        56*a^2*k^3*cos(phi)^2*sin(phi) - 108*k*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 + ...
        432*k*pi*sin((3*pi*sin(phi))/2)*cos(phi)^4 + 243*k*pi^3*sin((3*pi*sin(phi))/2)*cos(phi)^6 +...
        72*a^2*k^3*sin(phi)*cos((3*pi*sin(phi))/2) + 72*a^2*k*l^2*sin(phi)*cos((3*pi*sin(phi))/2) +...
        56*a^2*k*l^2*cos(phi)^2*sin(phi) - 108*a^2*k^3*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2 +...
        810*k*pi^2*cos(phi)^4*sin(phi)*cos((3*pi*sin(phi))/2) - 16*Omega*a*k*cos(phi)^3*sin(phi) -...
        108*a^2*k*l^2*pi*sin((3*pi*sin(phi))/2)*cos(phi)^2)/(4*a^3*cos(phi)*(k^2 + l^2));
end
end
%% 有三个方程的四阶Runge-Kutta方法
function [X,Y,Z,T] = Runge_Kutta_4_2eq(f,g,p,t0,x0,y0,z0,h,t,k)
% 输入参数说明：
% x'=f(t,x,y,z)；y'=g(t,x,y,z)；z'=p(t,x,y,z)；x=x(y)；y=y(t)；z=z(t)；
% t0为初始t值
% y0为初始t值对于的y值，y0=y(y0)；z0为初始t值对于的z值，z0=z(z0)
% h为积分步长
% t为积分取值
% k为波数
% --------------------------------------------------------------------------------------------
% 输出参数说明：
% X为x0到x区间上等间距h的数组
% Y(i)=y(T(i))

n = floor((t-t0)/h);%求步数
T = zeros(1,n+1);
X = zeros(1,n+1);
Y = zeros(1,n+1);
Z = zeros(1,n+1);
T(1) = t0;X(1) = x0;Y(1) = y0;Z(1) = z0;

for i = 1:n
    k1 = h*f(T(i),X(i),Y(i),Z(i),k);
    l1 = h*g(T(i),X(i),Y(i),Z(i),k);
    m1 = h*p(T(i),X(i),Y(i),Z(i),k);

    k2 = h*f(T(i)+h/2,X(i)+k1/2,Y(i)+l1/2,Z(i)+m1/2,k);
    l2 = h*g(T(i)+h/2,X(i)+k1/2,Y(i)+l1/2,Z(i)+m1/2,k);
    m2 = h*p(T(i)+h/2,X(i)+k1/2,Y(i)+l1/2,Z(i)+m1/2,k);

    k3 = h*f(T(i)+h/2,X(i)+k2/2,Y(i)+l2/2,Z(i)+m2/2,k);
    l3 = h*g(T(i)+h/2,X(i)+k2/2,Y(i)+l2/2,Z(i)+m2/2,k);
    m3 = h*p(T(i)+h/2,X(i)+k2/2,Y(i)+l2/2,Z(i)+m2/2,k);

    k4 = h*f(T(i)+h,X(i)+k3,Y(i)+l3,Z(i)+m3,k);
    l4 = h*g(T(i)+h,X(i)+k3,Y(i)+l3,Z(i)+m3,k);
    m4 = h*p(T(i)+h,X(i)+k3,Y(i)+l3,Z(i)+m3,k);

    X(i+1) = X(i)+(k1+2*k2+2*k3+k4)/6*180/pi;
    Y(i+1) = Y(i)+(l1+2*l2+2*l3+l4)/6*180/pi;
    Z(i+1) = Z(i)+(m1+2*m2+2*m3+m4)/6;
    T(i+1) = T(i)+h;
    %按照龙格库塔方法进行数值求解
end

end