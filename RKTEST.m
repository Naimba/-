clear;clc;close all
% 有两个方程的四阶Runge-Kutta方法的验证

f = @(t,x,y) x+3.*y;
g = @(t,x,y) -x-y;
t0 = 0;x0 = 2;y0=3;t = sqrt(2)*pi;h = 0.001*pi;

[Y,X,T] = Runge_Kutta_4_2eq(f,g,t0,x0,y0,h,t);
% disp(['R-K计算出的x = ',num2str(X(end))])
% disp(['R-K计算出的y = ',num2str(Y(end))])
x = @(t) x0.*cos(sqrt(2).*t)+(x0+3*y0)./sqrt(2).*sin(sqrt(2).*t);
y = @(t) y0.*cos(sqrt(2).*t)-(x0+y0)./sqrt(2).*sin(sqrt(2).*t);


X0 = x(T);Y0 = y(T);
plot(X,Y,'color','b','linewidth',1,'LineStyle','-.');hold on
% figure
plot(X0,Y0,'color','r','linewidth',0.8,'LineStyle',':')
legend("R-K值",'实际值')
% disp(' ')
% disp(['实际上的x = ',num2str(x)])
% disp(['实际上的y = ',num2str(y)])
%% 有两个方程的四阶Runge-Kutta方法
function [Z,Y,T] = Runge_Kutta_4_2eq(f,g,t0,y0,z0,h,t)
% 输入参数说明：
% y'=f(t,y,z)；z'=g(t,y,z)；y = y(t)；z = z(t)；
% t0为初始t值
% y0为初始t值对于的y值，y0=y(y0)；z0为初始t值对于的z值，z0=z(z0)
% h为积分步长
% t为积分取值
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
    k1 = h*f(T(i),Y(i),Z(i));
    l1 = h*g(T(i),Y(i),Z(i));
    k2 = h*f(T(i)+h/2,Y(i)+k1/2,Z(i)+l1/2);
    l2 = h*g(T(i)+h/2,Y(i)+k1/2,Z(i)+l1/2);
    k3 = h*f(T(i)+h/2,Y(i)+k2/2,Z(i)+l2/2);
    l3 = h*g(T(i)+h/2,Y(i)+k2/2,Z(i)+l2/2);
    k4 = h*f(T(i)+h,Y(i)+k3,Z(i)+l3);
    l4 = h*g(T(i)+h,Y(i)+k3,Z(i)+l3);
    Y(i+1) = Y(i)+(k1+2*k2+2*k3+k4)/6;
    Z(i+1) = Z(i)+(l1+2*l2+2*l3+l4)/6;
    T(i+1) = T(i)+h;
    %按照龙格库塔方法进行数值求解
end

end