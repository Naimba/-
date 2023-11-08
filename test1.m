clear;clc;close all;

%%
clear;clc;close all;
T = [inf,30*86400,-30*86400];a = 6.4e6;beta_M = 2e-11;
k1 = 1:15;k = k1 /a;% k是数学意义上的波数，k=2pi/omega；k1是纬向波数,k1=ka；
TL = tiledlayout(1,3);Text = {'(a)','(b)','(c)'};
for i = 1:3
    u_up = 2*pi./(T(i)*k);
    u_down = 2*pi./(T(i)*k)+beta_M./(k.^2);
    nexttile
    plot(k1,u_down,'LineWidth',1,'Color','r')
    hold on
    plot(k1,u_up,'LineStyle','--','LineWidth',1,'Color','b')
    xlim([1 15]);ylim([-20 100]);xticks(1:2:15);yticks(-20:10:100);
    fill([k1,fliplr(k1)],[u_down,fliplr(u_up)],'k','FaceAlpha',0.2,'EdgeColor','none')
    set(gca,'XTickLabelRotation',0)
    text(12,95,Text{i})
end
% title(TL,'不同纬向波数波动传播的纬向基流范围')
xlabel(TL,'Zonal Wave Number')
ylabel(TL,'Zonal Basic flow/(m s^{-1})')
print(gcf,'F:\学习\毕业论文\复现李艳杰\1','-dpng','-r400');
close

%%
clear;clc;close all;
u_M = [0 10 -10];a = 6.4e6;beta_M = 2e-11;
k1 = 1:15;k = k1 /a;% k是数学意义上的波数，k=2pi/omega；k1是纬向波数,k1=ka；
TL = tiledlayout(1,3);Text = {'(a)','(b)','(c)'};
for i = 1:3
    c_up = repelem(u_M(i),length(k1));
    c_down = u_M(i)-beta_M./(k.^2);
    nexttile
    plot(k1,c_down,'LineWidth',1,'Color','r')
    hold on
    plot(k1,c_up,'LineStyle','--','LineWidth',1,'Color','b')
    xlim([1 15]);ylim([-100 20]);xticks(1:2:15);yticks(-100:10:20);
    fill([k1,fliplr(k1)],[c_down,fliplr(c_up)],'k','FaceAlpha',0.2,'EdgeColor','none')
    set(gca,'XTickLabelRotation',0)
    text(12,15,Text{i})
end
% title(TL,'不同纬向波数波动传播的相速度范围')
xlabel(TL,'Zonal Wave Number')
ylabel(TL,'c/(m s^{-1})')

print(gcf,'F:\学习\毕业论文\复现李艳杰\2','-dpng','-r400');
close
%%
clear;clc;close all;
u_M = [0 10 -10];a = 6.4e6;beta_M = 2e-11;
k1 = 1:15;k = k1 /a;% k是数学意义上的波数，k=2pi/omega；k1是纬向波数,k1=ka；
TL = tiledlayout(1,3);

nexttile% 东西风临界线
T_down = 2*pi./(beta_M./k)/86400;
plot(k1,T_down,'LineWidth',1,'Color','r')
hold on
fill([k1,fliplr(k1)],[T_down,fliplr(repelem(9,length(k1)))],'k','FaceAlpha',0.2,'EdgeColor','none')
set(gca,'XTickLabelRotation',0)
xlim([1 15]);ylim([0 9]);xticks(1:2:15);yticks(1:9);
text(13,8.75,'(a)')
%%%%%%%%%%%%%%%%%%

nexttile
% 纬向波数较大
T_down = 2.*pi./(u_M(2).*k)./86400;
T_up = 2.*pi./(u_M(2).*k-beta_M./k)./86400;
plot(k1,T_down,'LineWidth',1,'Color','r')
hold on
plot(k1(T_up>0),T_up(T_up>0),'LineStyle','--','LineWidth',1,'Color','b')

% 纬向波数较小 c>0
% T_down = 2.*pi./(beta_M./k)./86400;
% plot(k1,T_down,'LineWidth',1,'Color','g')

% % 纬向波数较小 c<0
T_down1 = -2.*pi./(u_M(2).*k-beta_M./k)./86400;
plot(k1(1:8),T_down1(1:8),'LineWidth',1,'Color','k')
% fill([k1,fliplr(k1)],[T_down1(1:6),T_down(7:end),T_up(end:-1:9),repelem(50,7)],'k','FaceAlpha',...
    % 0.2,'EdgeColor','none')
fill([k1(1:6),k1(6:-1:1)],[T_down1(1:6),repelem(50,6)],'k','FaceAlpha',0.2,'EdgeColor','none')
fill([k1(6:end),k1(end:-1:6)],[T_down(6:end),T_up(end:-1:10),repelem(50,4)],'k','FaceAlpha',0.2,'EdgeColor','none')
fill([6,6.3404,6],[T_down1(6),7.3799,T_down(6)],'k','FaceAlpha',0.2,'EdgeColor','none')

set(gca,'XTickLabelRotation',0)
xlim([1 15]);ylim([0 50]);xticks(1:2:15);yticks(0:5:50);
text(13,48.5,'(b)')
text(2.5,8,'B');text(6,15,'A');text(9,10,'C');

%%%%%%%%%%%%%%%%%%%%%
nexttile
T_up = -2.*pi./(u_M(3).*k)./86400;
T_down = -2.*pi./(u_M(3).*k-beta_M./k)./86400;
plot(k1,T_down,'LineWidth',1,'Color','r')
hold on
plot(k1,T_up,'LineStyle','--','LineWidth',1,'Color','b')
fill([k1,fliplr(k1)],[T_down,fliplr(T_up)],'k','FaceAlpha',0.2,'EdgeColor','none')
set(gca,'XTickLabelRotation',0)
xlim([1 15]);ylim([0 50]);xticks(1:2:15);yticks(0:5:50);
text(13,48.5,'(c)')

set(gcf,'Position',[139,225,906.4,420])

xlabel(TL,'Zonal Wave Number')
ylabel(TL,'T/d')
print(gcf,'F:\学习\毕业论文\复现李艳杰\3','-dpng','-r400');
close


