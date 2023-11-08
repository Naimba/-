clear;clc;close all
a = 6400000;beta_M = 2e-11;d = 86400;h = 2*d;
for j = 1:5
    phi = -90:0.5:90;k = j/a;l = zeros(1,length(phi));
    for i = 1:length(phi)
        l(i) = l_num(phi(i),k);
    end
    l = l.^2;
    plot(phi,l,'LineWidth',1,'Color','r');
    xlim([-90 90]);xticks(-90:30:90);
    xticklabels({'90\circ S','60\circ S','30\circ S','EQ','30\circ N','60\circ N','90\circ N'})
    set(gca,'XTickLabelRotation',0)
    set(gca, 'GridLineStyle', ':','GridAlpha', 0.2,'MinorGridAlpha',0.2,...
        'XMinorGrid','on','YMinorGrid','on','LineWidth',0.8);
    title(['纬向',num2str(j),'波的l^2分布'])
    print(gcf,['F:\学习\毕业论文\复现李艳杰\纬向',num2str(j),'波的l2分布.png'],'-dpng','-r400');
    close
end
%% beta_M
function beta_M = beta_M_(phi)

syms phi0 u_M beta_M beta_M2
u_M = (18.*sind(3.*180./2.*(1+sind(phi0)))+14.*(1-sind(phi0).^2))/cosd(phi0);
beta_M = 2*7.292e-5*cosd(phi0)^2/6.4e6;
beta_M2 = cosd(phi0)/6.4e6*diff(1/6.4e6/cosd(phi0)*diff(u_M*cosd(phi0)^2,phi0),phi0);
beta_M = beta_M-beta_M2;
beta_M = double(subs(beta_M,{phi0},{phi*pi/180}));
% beta_M = 2e-11;
end
%% 经向波数
function l = l_num(phi,k)
% phi为纬度，单位为角度
% k为波数

beta_M = 2e-11;a = 6.4e6;
u_M0 = beta_M/k/k;%u_M0为转向临界风速
u_M = u_M_fun(phi);
v_M = v_M_fun(phi);
syms aa

if v_M == 0
    l = sqrt((beta_M-u_M_fun(phi)*k^2)/u_M_fun(phi));
else
    l = solve(v_M*aa^3+u_M*k*aa^2+v_M*k^2*aa+u_M*k^3-beta_M*k==0,aa,'Real',true);% 经向波数
end

if u_M>u_M0
    l = -l;
end
if ~isreal(l)
    l =nan;
end
end
%% 群速度ug vg
function ug = ugroup(~,~,phi,k)
% ugroup(t,lambda,phi,k,l)
% phi为纬度，单位为角度
% k为波数

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

beta_M = beta_M_(phi);a = 6.4e6;
l = l_num(phi,k);

ug = (v_M_fun(phi)*l/k+2*beta_M*k^2/((k^2+l^2)^2))/(a*cosd(phi));

% ug = ((1+(k^2)/(k^2+l^2))*u_M_fun(phi)+k*l/(k^2+l^2)*v_M_fun(phi)+ ...
%     2*beta_M*k*k/((k^2+l^2)^2))/(a*cosd(phi));
end
%-------------------------------------------------
function vg = vgroup(~,~,phi,k)
% phi为纬度，单位为角度
% k为波数

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

beta_M = beta_M_(phi);a = 6.4e6;
l = l_num(phi,k);


vg = (v_M_fun(phi)+2.*k.*l.*beta_M./((k.^2+l.^2).^2))/(a);

% vg = (k*l/(k^2+l^2)*u_M_fun(phi)+(1+(l^2)/(k^2+l^2))*v_M_fun(phi)+ ...
% 2.*k.*l.*beta_M./((k.^2+l.^2).^2))/(a);

end

%% 基流
function u_M = u_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

% u_M = (6.4e6)*30.875*(7.292e-5).*cosd(phi)/1000;
u_M = (18.*sind(3.*180./2.*(1+sind(phi)))+14.*(1-sind(phi).^2));
u_M = u_M./cosd(phi);
end

function v_M = v_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

v_M = 0;
% v_M = zeros(1,length(phi));phi0 = asind(0.3);
% for i = 1:length(phi)
%     if phi(i)<=phi0 && phi(i)>=asind(-0.5)
%         v_M(i) = 3.2.*sind(180.*(sind(phi(i))-sind(phi0))/(sind(phi0)+0.5));
%     end
%     if phi(i)<=asind(0.5) && phi(i)>phi0
%         v_M(i) = 0.8*sind(180*(sind(phi(i))-sind(phi0))/(0.5-sind(phi0)));
%     end
% end
% v_M = v_M./cosd(phi);
end