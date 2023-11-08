clear;clc;close all
% 找出转向纬度-hoskins

%% beta_M的构造
a = 6.4e6;
k = 1/a;
phi = 0:0.1:90;

syms phi0 u_M beta_M  a  Omega
% u_M = a*omega;
u_M = (18.*sin(3.*pi./2.*(1+sin(phi0)))+14.*(1-sin(phi0).^2))/cos(phi0);
beta_M = 2*Omega*cos(phi0)^2/a-cos(phi0)/a*diff(1/a/cos(phi0)*diff(u_M*cos(phi0)^2,phi0),phi0);
beta_M = subs(beta_M,{a,Omega},{6.4e6,7.292e-5});

beta_M0 = zeros(1,length(phi));
for i = 1:length(phi)-1
    beta_M0(i) = double(subs(beta_M,{phi0},{phi(i)*pi/180}));
end
beta_M0(end) = 0;

%% uM
u_M = (18.*sind(3.*180./2.*(1+sind(phi)))+14.*(1-sind(phi).^2))./cosd(phi);
u_M0 = beta_M0./k^2;
plot(phi,u_M,'LineWidth',1,'Color','r');hold on
plot(phi,u_M0,'LineWidth',0.5,'Color','b','LineStyle','-.')
legend('实际速度','临界速度')
xlabel('$\varphi$','Interpreter','latex')
ylabel([' $u_M$ (m/s)'],'Interpreter','latex','Rotation',0)
xlim([0 90])
ylim([-10 40])
yticks(-10:5:40)

%% 求出转向纬度解析解
syms phi u_M beta_M  a  Omega k
k = 1/a;
u_M = (18.*sin(3.*pi./2.*(1+sin(phi)))+14.*(1-sin(phi).^2))/cos(phi);
beta_M = 2*Omega*cos(phi)^2/a-cos(phi)/a*diff(1/a/cos(phi)*diff(u_M*cos(phi)^2,phi),phi);
eq = beta_M./k^2 == u_M;
eq = subs(eq,{a,Omega},{6.4e6,7.292e-5});
phi0 = vpasolve(eq,phi,[0 pi/2]);
phi0 = double(phi0*180/pi);
%% 求出截陷纬度解析解
syms phi u_M
u_M = (18.*sin(3.*pi./2.*(1+sin(phi)))+14.*(1-sin(phi).^2))/cos(phi);
phi0 = vpasolve(u_M == 0,phi,[0 pi/2]);
phi0 = double(phi0*180/pi);