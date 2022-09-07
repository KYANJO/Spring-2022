
close all 
clear

m = [-0.4326 -1.6656 0.1253 0.2877 -1.1465 1.1909 1.1892 ...
    -0.0376 0.3273 0.1746]

mbar = sum(m)/length(m)

mean(m)

s = std(m)

s1 = 0;

for i = 1:length(m)
    s1 = s1 + (m(i) - mbar)^2;
end
s2 = sqrt(s1/(length(m)-1))

N = 1000;
n5 = 5; n50 =50;

[m5] = generate(N,n5);
[m50] = generate(N,n50);

% figure(1)
% hist(m5)
% 
% figure(2)
% qqplot(m5)
% 
% 
% figure(3)
% hist(m50)
% 
% figure(4)
% qqplot(m50)

t = [3.4935 4.2853 5.1374 5.8181 6.8632 8.1841]';

x = [6 10.1333 14.2667 18.4000 22.5333 26.6667]';

m = length(t)

G = [ones(m,1) x]; % matrix G
M_L2 = inv(G'*G)*G'*t % least square solution

sig = 0.1;
W = (1/sig)*eye(m);

a = inv(G'*(W^2)*G);
b = G'*(W^2)*t;
a*b
M_mle = inv(G'*(W^2)*G)*G'*(W^2)*t

C = (sig^2)*inv(G'*G)
Delta = chi2inv(0.95,4)

plot_ellipse(Delta,C,M_L2)
title('Error Ellipsoid')
xlabel('t_o(s)'); ylabel('s_2 (s/km)')
%ylim([0.1,0.45])

function plot_ellipse(DELTA2,C,m)
n=5000;
%construct a vector of n equally-spaced angles from (0,2*pi)
theta=linspace(0,2*pi,n)';
%corresponding unit vector
xhat=[cos(theta),sin(theta)];
Cinv=inv(C);
%preallocate output array
r=zeros(n,2);
for i=1:n
%store each (x,y) pair on the confidence ellipse
%in the corresponding row of r
r(i,:)=sqrt(DELTA2/(xhat(i,:)*Cinv*xhat(i,:)'))*xhat(i,:);
end

% Plot the ellipse and set the axes.
plot(m(1)+r(:,1), m(2)+r(:,2));
axis equal
end


function [ave] = generate(N,n)
    ave = [];
    for i=1:N
        m = mean(exprnd(10,n,1)); 
        ave = [ave m];
    end
end
