clear;clc;close all
% ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01]);
%  for ii = 1:6 
%      axes(ha(ii)); 
%      plot(randn(10,ii)); 
%      hold on
%      plot([1 1],[0 1])
%  end
%  set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

tic
ha = tight_subplot(3,1,[.05 .03],[.05 .01],[.01 .02]);
axes(ha(1));
test3_edit
m_text(225,80,'(a)')

axes(ha(2));
test4_edit
m_text(225,80,'(b)')

axes(ha(3));
test5
m_text(225,80,'(c)')

set(gcf,"Position",[514.2,65.8,533.8,707.2])
toc

print(gcf,'F:\学习\毕业论文\复现李艳杰\定常波射线及传播距离','-dpng','-r400');
close