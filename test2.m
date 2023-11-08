% 水平非均匀基流中波动传播理论——波射线理论
clear;clc;close all

%% 构造理想基本气流函数
phi = -90:0.5:90;phi0 = asind(0.3);
u_M = u_M_fun(phi);v_M = v_M_fun(phi);

% 以下是画图代码
TL = tiledlayout(1,2);
set(gcf,'Position',[326.2,195.4,862.8,482])
nexttile
plot(phi,u_M,'LineWidth',1,'Color','r');
xlim([-90 90]);xticks(-90:30:90);ylim([-10 30]);yticks(-5:5:30)
xticklabels({'90\circ S','60\circ S','30\circ S','EQ','30\circ N','60\circ N','90\circ N'})
set(gca,'XTickLabelRotation',0)
set(gca, 'GridLineStyle', ':','GridAlpha', 0.2,'MinorGridAlpha',0.2,...
    'XMinorGrid','on','YMinorGrid','on','LineWidth',0.8);
h = ylabel('m\cdot s^{-1}');
h.Rotation = 0;
h.Position = [-103.75,30.47,-1];
% title('u','fontweight','bold','FontAngle',"italic")
text(76,28,'(a)')

nexttile
plot(phi,v_M,'LineWidth',1,'Color','r');
xlim([-90 90]);xticks(-90:30:90);ylim([-4 1.5]);yticks(-4:0.5:1.5)
xticklabels({'90\circ S','60\circ S','30\circ S','EQ','30\circ N','60\circ N','90\circ N'})
set(gca,'XTickLabelRotation',0)
set(gca, 'GridLineStyle', ':','GridAlpha', 0.2,'MinorGridAlpha',0.2,...
    'XMinorGrid','on','YMinorGrid','on','LineWidth',0.8);
h = ylabel('m\cdot s^{-1}');
h.Rotation = 0;
h.Position = [-103.75,1.6,-1];
% title('v','fontweight','bold','FontAngle',"italic")
text(76,1.3,'(b)')

print(gcf,'F:\学习\毕业论文\复现李艳杰\构造理想基本气流','-dpng','-r400');
close
%% 基流
function u_M = u_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

u_M = (18.*sind(3.*180./2.*(1+sind(phi)))+14.*(1-sind(phi).^2));
% u_M = u_M./cosd(phi);
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
% v_M = v_M./cosd(phi);
% v_M =0;
end
%% 四阶Runge-Kutta方法
% function [Y,X] = Runge_Kutta_4(f,x0,y0,x,h,k,l)
% % 输入参数说明：
% % fun为函数形式y'=fun(x,y),y=y(x)
% % x0为初始x值
% % y0为初始x值对于的y值，y0=y(x0)
% % x为积分取值
% % h为积分步长
% % --------------------------------------------------------------------------------------------
% % 输出参数说明：
% % X为x0到x区间上等间距h的数组
% % Y(i)=y(X(i))
%
%
%
% n = floor((x-x0)/h);%求步数
% X = zeros(1,n+1);
% Y = zeros(1,n+1);
% X(1) = x0;Y(1) = y0;
%
% for i = 1:n
%     k1 = f(X(i),Y(i),k,l);
%     k2 = f(X(i)+h/2,Y(i)+h*k1/2,k,l);
%     k3 = f(X(i)+h/2,Y(i)+h*k2/2,k,l);
%     k4 = f(X(i)+h,Y(i)+h*k3,k,l);
%     Y(i+1) = Y(i)+h*(k1+2*k2+2*k3+k4)/6*180/pi;
%     X(i+1) = X(i)+h;
%     %按照龙格库塔方法进行数值求解
% end
%
% end
%% 有两个方程的四阶Runge-Kutta方法
function [Z,Y,T] = Runge_Kutta_4_2eq(f,g,t0,y0,z0,h,t,k)
% 输入参数说明：
% y'=f(t,y,z)；z'=g(t,y,z)；y = y(t)；z = z(t)；但是在这里ugroup(t,lambda,phi,k,l)
% t0为初始t值
% y0为初始t值对于的y值，y0=y(y0)；z0为初始t值对于的z值，z0=z(z0)
% h为积分步长
% t为积分取值
% k，l为波数
% --------------------------------------------------------------------------------------------
% 输出参数说明：
% X为x0到x区间上等间距h的数组
% Y(i)=y(T(i))



n = floor((t-t0)/h);%求步数
T = zeros(1,n+1);
Y = zeros(1,n+1);
Z = zeros(1,n+1);
T(1) = t0;Y(1) = y0;Z(1) = z0;

for i = 1:n
    k1 = h*f(T(i),Y(i),Z(i),k);l1 = h*g(T(i),Y(i),Z(i),k);
    k2 = h*f(T(i)+h/2,Y(i)+k1/2,Z(i)+l1/2,k);
    l2 = h*g(T(i)+h/2,Y(i)+k1/2,Z(i)+l1/2,k);
    k3 = h*f(T(i)+h/2,Y(i)+k2/2,Z(i)+l2/2,k);
    l3 = h*g(T(i)+h/2,Y(i)+k2/2,Z(i)+l2/2,k);
    k4 = h*f(T(i)+h,Y(i)+k3,Z(i)+l3,k);
    l4 = h*g(T(i)+h,Y(i)+k3,Z(i)+l3,k);
    Y(i+1) = Y(i)+(k1+2*k2+2*k3+k4)/6*180/pi;
    Z(i+1) = Z(i)+(l1+2*l2+2*l3+l4)/6*180/pi;
    T(i+1) = T(i)+h;
    %按照龙格库塔方法进行数值求解
end

end