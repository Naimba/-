clear;clc;close all
% 找出转向纬度-考虑v_M
ii = 1;phi00 = asind(0.3);
%% 传播解的个数
% syms phi u_M beta_M  a  Omega k Delta v_M
% u_M = (18.*sin(3.*pi./2.*(1+sin(phi)))+14.*(1-sin(phi).^2))/cos(phi);
% phi0 =phi*pi/180; k = ii/a;
% beta_M= (2690273155709801*cos(phi0)^2)/118059162071741130342400000 +...
%     (cos(phi0)*((cos(phi0)*(28*cos(phi0)^2 - 28*sin(phi0)^2 + ...
%     (81*pi^2*sin((3*pi*(sin(phi0) + 1))/2)*cos(phi0)^2)/2 + 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*sin(phi0))...
%     - 2*sin(phi0)*(28*cos(phi0)*sin(phi0) - 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*cos(phi0)) + ...
%     cos(phi0)*(18*sin((3*pi*(sin(phi0) + 1))/2) - 14*sin(phi0)^2 + 14))/(6400000*cos(phi0)) + ...
%     (sin(phi0)*(sin(phi0)*(18*sin((3*pi*(sin(phi0) + 1))/2) - 14*sin(phi0)^2 + 14) + ...
%     cos(phi0)*(28*cos(phi0)*sin(phi0) - 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*cos(phi0))))...
%     /(6400000*cos(phi0)^2)))/6400000;
% % % abs(phi)>30 v_M=0
% disp('abs(phi)>30 v_M=0')
% v_M = 0;
% disp('退化为一元二次方程')% 无解，全是三个不等的实根
% % % phi00<phi<30 v_M~=0
% disp('phi00<phi<30 v_M~=0')
%  y00 = 0.3;y = sin(phi);
% v_M = 0.8.*sin(pi.*(phi-y00)./(0.5-y00));
% Delta = (27*v_M^2*(- u_M*k^3 + beta_M*k) - 2*k^3*u_M^3 ...
%     + 9*k^3*u_M*v_M^2)^2/(2916*v_M^6) - (k^2*u_M^2 - 3*k^2*v_M^2)^3/(729*v_M^6);
% Delta = subs(Delta,{a},{6.4e6});
% vpasolve(Delta==0,phi,[asin(0.3) pi/6])*180/pi;
% disp('无解')% 无解，全是三个不等的实根
% % % -30<phi<phi00 v_M~=0
% disp('-30<phi<phi00 v_M~=0')
% v_M = 3.2.*sin(pi.*(phi-y00)./(0.5+y00));
% Delta = (27*v_M^2*(- u_M*k^3 + beta_M*k) - 2*k^3*u_M^3 ...
%     + 9*k^3*u_M*v_M^2)^2/(2916*v_M^6) - (k^2*u_M^2 - 3*k^2*v_M^2)^3/(729*v_M^6);
% Delta = subs(Delta,{a},{6.4e6});
% double(vpasolve(Delta==0,phi,[-pi/6 0])*180/pi)% 
% double(vpasolve(Delta==0,phi,[0 asin(0.3)])*180/pi)% 

%% Delta作图
phi = -29.9:0.1:29.9;Delta = zeros(1,length(phi));a = 6.4e6;k = ii/a;
for i = 1:length(phi)
    v_M = v_M_fun(phi(i));u_M = u_M_fun(phi(i));
    phi0 =phi(i)*pi/180; 
    beta_M= (2690273155709801*cos(phi0)^2)/118059162071741130342400000 +...
        (cos(phi0)*((cos(phi0)*(28*cos(phi0)^2 - 28*sin(phi0)^2 + ...
        (81*pi^2*sin((3*pi*(sin(phi0) + 1))/2)*cos(phi0)^2)/2 + 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*sin(phi0))...
        - 2*sin(phi0)*(28*cos(phi0)*sin(phi0) - 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*cos(phi0)) + ...
        cos(phi0)*(18*sin((3*pi*(sin(phi0) + 1))/2) - 14*sin(phi0)^2 + 14))/(6400000*cos(phi0)) + ...
        (sin(phi0)*(sin(phi0)*(18*sin((3*pi*(sin(phi0) + 1))/2) - 14*sin(phi0)^2 + 14) + ...
        cos(phi0)*(28*cos(phi0)*sin(phi0) - 27*pi*cos((3*pi*(sin(phi0) + 1))/2)*cos(phi0))))...
        /(6400000*cos(phi0)^2)))/6400000;
    Delta(i) = (27*v_M^2*(- u_M*k^3 + beta_M*k) - 2*k^3*u_M^3 ...
        + 9*k^3*u_M*v_M^2)^2/(2916*v_M^6) - (k^2*u_M^2 - 3*k^2*v_M^2)^3/(729*v_M^6);
end
plot(phi,Delta,'LineWidth',1,'Color','r','LineStyle','-')
xlim([-30 30]);xticks(-30:10:30)
xticklabels({'30\circ S','20\circ S','10\circ S','EQ','10\circ N','20\circ N','30\circ N'})
ylim([-5e-36 5e-36]);hold on
plot([-30 30],[0 0],'Color','k','LineStyle','--')

%% l^2作图
figure
 a = 6.4e6;k = ii/a;
phi = -50:0.1:50;l = zeros(3,length(phi));
for i = 1:length(phi)
    l(:,i) = l_num(phi(i),k);
end
% % l2 = l.^2;
plot(phi,l(2,:),'LineWidth',1,'Color','r','LineStyle','-')
xlabel('$\varphi$','Interpreter','latex')
ylabel(['$l$'],'Interpreter','latex','Rotation',0)
xlim([-90 90]);xticks(-90:30:90)
xticklabels({'90\circ S','60\circ S','30\circ S','EQ','30\circ N','60\circ N','90\circ N'})
% ylim([-10 40])
% yticks(-10:5:40)
hold off
set(gca, 'GridLineStyle', ':','GridAlpha', 0.2,'MinorGridAlpha',0.2,...
    'XMinorGrid','on','YMinorGrid','on','LineWidth',0.8);
l(4,:)= phi;

%% 转向纬度
% % |phi|>30 v_M=0
% syms phi u_M beta_M  a  Omega k
% k = ii/a;
% % u_M = a*omega;
% u_M = (18.*sin(3.*pi./2.*(1+sin(phi)))+14.*(1-sin(phi).^2))/cos(phi);
% beta_M = 2*Omega*cos(phi)^2/a-cos(phi)/a*diff(1/a/cos(phi)*diff(u_M*cos(phi)^2,phi),phi);
% eq = beta_M./k^2 == u_M;
% eq = subs(eq,{a,Omega},{6.4e6,7.292e-5});
% alpha = vpasolve(eq,phi,[pi/6 pi/2]);
% alpha = double(alpha*180/pi);
% disp(['|phi|>30°时，转向纬度alpha = ',num2str(alpha)])
% 
% % phi00<phi<30 v_M~=0
% y00 = 0.3;y = sin(phi);
% v_M = 0.8.*sin(pi.*(phi-y00)./(0.5-y00));
% l = ((27*v_M^2*(- u_M*k^3 + beta_M*k) - 2*k^3*u_M^3 +...
%     9*k^3*u_M*v_M^2)/(54*v_M^3) - ((27*v_M^2*(- u_M*k^3 ...
%     + beta_M*k) - 2*k^3*u_M^3 + 9*k^3*u_M*v_M^2)^2/(2916*v_M^6) ...
%     - (k^2*u_M^2 - 3*k^2*v_M^2)^3/(729*v_M^6))^(1/2))^(1/3) +...
%     ((27*v_M^2*(- u_M*k^3 + beta_M*k) - 2*k^3*u_M^3 + ...
%     9*k^3*u_M*v_M^2)/(54*v_M^3) + ((27*v_M^2*(- u_M*k^3 + beta_M*k) ...
%     - 2*k^3*u_M^3 + 9*k^3*u_M*v_M^2)^2/(2916*v_M^6) - (k^2*u_M^2 - 3*k^2*v_M^2)^3/...
%     (729*v_M^6))^(1/2))^(1/3) - (k*u_M)/(3*v_M);
% vg = v_M+2*beta_M*k*l/((k^2+l^2)^2);
% vg = subs(vg,{a,Omega},{6.4e6,7.292e-5});
% alpha = vpasolve(vg==0,phi,[phi00*pi/180 pi/6]);
% alpha = double(alpha*180/pi);
% disp([num2str(phi00),'<phi<30°时，',newline,'转向纬度alpha = ',num2str(alpha)])
%% 经向波数
function l = l_num(phi,k)
% phi为纬度，单位为角度
% k为波数

beta_M = beta_M_(phi);
u_M = u_M_fun(phi);
v_M = v_M_fun(phi);

if v_M == 0
    l2 = (beta_M-u_M_fun(phi).*k.^2)/u_M_fun(phi);
    if l2 < 0
        l = nan;
        % error('需要转向了')
    else
        l = sqrt(l2);
    end
else
    disp(['此时phi = ',num2str(phi)])
    l = solve3(v_M,u_M*k,v_M*k^2,u_M*k^3-beta_M*k);
    
    % l = vpasolve(v_M*aa^3+u_M*k*aa^2+v_M*k^2*aa+u_M*k^3-beta_M*k==0,aa,[-1e-5 1e5]);% 经向波数
    % l = l(2);
end
end
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
% v_M = 0;
end

%% 群速度ug vg
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