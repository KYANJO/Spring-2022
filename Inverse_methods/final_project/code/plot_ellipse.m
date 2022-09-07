function plot_ellipse(DELTA2,C,m,first,second)
n=5000;
%first - first parameter
%second - second parameter
%construct a vector of n equally-spaced angles from (0,2*pi)
theta=linspace(0,2*pi,n)';

%corresponding unit vector
xhat=[cos(theta),sin(theta)];
Cinv=inv(C);

%preallocate output array
r=zeros(n,2);
for i=1:n
%store each (x,y) pair on the confidence ellipse in the corresponding row of r
r(i,:)=sqrt(DELTA2/(xhat(i,:)*Cinv*xhat(i,:)'))*xhat(i,:);
end

% Plot the ellipse and set the axes.
for i = 1:3
    a=['r' 'b' 'k'];
    plot(m(i,first)+r(:,1), m(i,second)+r(:,2),'color',a(i)); hold on
    plot(m(i,first),m(i,second),'.','color',a(i),'MarkerSize',17);grid on
    axis equal
end
end