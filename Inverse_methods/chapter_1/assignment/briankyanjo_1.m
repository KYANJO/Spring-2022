
close all
clc

mtrue = [0.18 16.21 9.81]';

m = 20; 
n = 3;

tj = linspace(0,3,N);

fprintf("N0.2 a) \n\n")

fprintf("(i). Form G \n ")

y = @(t) m1 + m2*t - 0.5*m3*t^2

G = zeros(m,n);

for i = 1:n
    for j = 1:m
       G(j,1) = 1;
       G(j,2) = tj(j);
       G(j,3) = -0.5*tj(j)^2;
    end
end 
fprintf
disp(G)

dtrue =  G*mtrue

%noise
noise = 2*randn(m,1);

%noise data
d = dtrue + noise

%plot
figure(1)
plot(dtrue,'b.','MarkerSize',20)
hold on
plot(d,'r.','MarkerSize',20)
title('A graph of d_{true} and d')
legend('dtrue','d')
xlabel('Number of Observations'); 
ylabel('d_{true} & d')


