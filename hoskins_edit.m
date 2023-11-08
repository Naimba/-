% clear;clc;close all
% 使用dl/dt的方法计算Hoskins均匀基流
%% dl的构造
% syms y phi k a Omega l omega
% u_M = a*omega;
% beta_M = 2*Omega*cos(phi)^2/a-cos(phi)/a*diff(1/a/cos(phi)*diff(u_M*cos(phi)^2,phi),phi);
% v_M = 0;
% dl = -k*cos(phi)/a*diff(u_M,phi)-l*cos(phi)/a*diff(v_M,phi)+cos(phi)/a*diff(beta_M,phi)*k/(k^2+l^2);
% dl = simplify(dl);
%% Hoskins的计算
% tic
% tiledlayout(1,1)
nexttile
m_proj('Equidistant Cylindrical','lon',[0 360],'lat',[-90 90]);
m_grid('box','on','tickdir','in','xtick',0:30:360,'ytick',-90:30:90);
m_coast('line','color','k','linewidth',0.2,'linestyle','-');hold on
for ks = [1 -1]
    for ii = 1:5
        a = 6.4e6;k =ii/a;
        d = 86400;Omega = 7.292e-5;
        n = 1200;
        h = 0.1*d;%积分步长


        ug = @ugroup;vg = @vgroup;dl = @dl_fun;
        phi = 0;lambda = 0;
        lat = zeros(1,n+1);lon = zeros(1,n+1);
        lat(1) = phi;lon(1)=lambda;
        beta_M = beta_M_(phi);
        l = ks*sqrt((beta_M-u_M_fun(phi)*k^2)/u_M_fun(phi));

        i = 1;

        while i<=3000
            [Lambda,Phi,L,T] = Runge_Kutta_4_2eq(ug,vg,dl,0,lambda,phi,l,h,h,k);
            phi = Phi(2);
            lambda = Lambda(2);l = L(2);
            % disp(['第',num2str(i,'%2d'),'次积分完成，当前积分第',num2str(h/d*i,'%5.2f'),'天'])
            lat(i+1) = phi;lon(i+1)=lambda;
            i = i+1;
        end

        lat = lat(1:i-1);lon = lon(1:i-1);

        h1 = m_plot(lon(1:end),lat(1:end),'color','r','linewidth',1,'marker','none');hold on
        m_scatter(lon(1:2*d/h:end),lat(1:2*d/h:end),8,'k','filled')

        % set(gcf,'Position',[286.6,237,854,498.4])
        if ks == 1
            if ii>=3
                [~,ind]=max(lat(1:200));
                m_text(90,lat(ind)+4,num2str(ii));m_text(270,lat(ind)+4,num2str(ii));
            else
                [~,ind]=max(lat(1:400));
                m_text(90,lat(ind)+4,num2str(ii));m_text(270,lat(ind)+4,num2str(ii))
            end
        else
            if ii>=3
                [~,ind]=min(lat(1:200));
                m_text(90,lat(ind)-4,num2str(ii));m_text(270,lat(ind)-4,num2str(ii));
            else
                [~,ind]=min(lat(1:400));
                m_text(90,lat(ind)-4,num2str(ii));m_text(270,lat(ind)-4,num2str(ii));
            end
        end
    end
end
% toc
%% 理论上的波射线
for ks = [1 -1]
    for ii = 1:5
        k = ii/a;
        alpha = acosd(a*k*sqrt(1/63.75));
        lambda0 = 0;
        lambda = 0:360;
        phi = ks*atand(tand(alpha).*sind(lambda-lambda0));
        h2 = m_plot(lambda,phi,'color','b','linewidth',1,'linestyle',':');
    end
end
m_text(345,80,'(b)')
%% 图题等等
% title('Hoskins论文中的定常波射线及传播距离','fontsize',12,'fontname','Yahei','Interpreter','tex','FontAngle','italic')
% h = legend([h1 h2],{'计算解','理论解'});
% h.Position = [0.74,0.8587,0.204,0.04];h.Orientation = "horizontal";
% print(gcf,'F:\学习\毕业论文\复现李艳杰\Hoskins论文中的定常波射线及传播距离','-dpng','-r400');
% close
%% beta_M
function beta_M = beta_M_(phi)

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

a = 6400000;
beta_M = 2*cosd(phi)^2/a*(7.292e-5+1/30.875*7.292e-5);

end
%% 基流
function u_M = u_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

u_M = 6.4e6/30.875*7.292e-5;
end

function v_M = v_M_fun(phi)
% phi为纬度，单位为角度

if abs(phi)>90
    error('latitude should be in the range [-90 90]');
end

% v_M = 0;
% v_M = zeros(1,length(phi));
% for i = 1:length(phi)
%     if y(i)<=y0 && y(i)>=(-0.5)
%         v_M(i) = 3.2.*sind(180.*(y(i)-y0)./(y0+0.5));
%     end
%     if y(i)<=(0.5) && y(i)>y0
%         v_M(i) = 0.8.*sind(180.*(y(i)-y0)./(0.5-y0));
%     end
% end
% v_M = v_M./cosd(phi);
v_M =0;
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



ug = (v_M_fun(phi).*l./k+2*beta_M*k^2/((k.^2+l.^2).^2))/(a*cosd(phi)).*cosd(phi);


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

a = 6.4e6;Omega = 7.292e-5;omega = 1/30.875*Omega;
phi = phi*pi/180;
dl = -(4*k*(Omega + omega)*(sin(phi) - sin(phi)^3))/(a^2*(k^2 + l^2));
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