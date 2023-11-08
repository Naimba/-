function [x1,x2,x3] = solve3(a,b,c,d)
% ax^3+bx^2+cx+d = 0
%一元三次方程的求根公式———精确解

p = (3*a*c-b^2)/(3*a^2);
q = (27*(a^2)*d-9*a*b*c+2*(b^3))/(27*a^3);
omega =  complex(-1,sqrt(3))/2;

Delta = q^2/4+p^3/27;

r1 = -b/(3*a);
r2 = (-q/2+sqrt(Delta))^(1/3);
r3 = (-q/2-sqrt(Delta))^(1/3);


x1 = r1+r2+r3;
x2 = r1+omega*r2+omega^2*r3;
x3 = r1+omega^2*r2+omega*r3;

if Delta > 0
    % x3 = real(x3);
    disp('有一个实根和两个共轭虚根')
else
    x1 = real(x1);x2 = real(x2);x3 = real(x3);
    if Delta == 0
         if q==0 && p==0
            disp('三个实根都相等')
         else
             disp('三个实根中有两个相等')
         end
    else
        disp('有三个不相等的实根')
    end
end
        
X = [x1 x2 x3];
end