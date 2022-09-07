%%====================================================================
% Boise State University
% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
% inverse model script
%%=====================================================================

clc
close all
global g; 
warning('off','all')

% problem
hl = 2;                       % left depth
hr = 1;                       % right depth
ul = 0;                       % left velocity
ur = 0;                       % right velocity

% Spatial domain
ax = -5;
bx = 5;
ay = -2;
by = 4;
meqn = 2;                       % Number of equations in the system


g = 1;                          % Gravity 

to = 0;                         % initial time
Tfinal = 1;                     % final time

ql = [hl; hl*ul];            % left conservation variable
qr = [hr; hr*ur];           % right consrvation variable

N = 10;                        % Number of spartial steps

dx = (bx - ax)/N;               % spartial step size

cfl = 0.9;                      % cfl number

a = 1.5;                        % maximum velocity

dt_est = cfl*dx/a;

M = (floor(Tfinal/dt_est) + 1); % number of time steps

dt = Tfinal/(M);                % temporal step size

t = linspace(to,Tfinal,M);      % temporal domain

xe = linspace(ax,bx,N+1);       % edge locations
x = xe(1:end-1) + dx/2;         % Cell-center locations

mq = 4;                         % 1-height, 
                                % 2-velocity, 
                                % 3-momentum, 
                                % 4-all fields
        
% %Initial conditions
h = zeros(N,1);
hu = zeros(N,1);
hn = zeros(N,1);
un = zeros(N,1);

for i=1:N-1
    if (x(i)<0) 
        h(i) = ql(1);
        hu(i) = ql(2);
    else
        h(i) = qr(1);
        hu(i) = qr(2);
    end
end

% writing a video
% v = VideoWriter('dam.avi');
% open(v);
% 
% 
% close(v);


% Inverse Riemann solver

for n=1:M
    for i=1:N
       if i == N
           q1l = [h(i); hu(i)];
           q1r = [h(i); hu(i)];
       else
           q1l = [h(i); hu(i)];
           q1r = [h(i+1); hu(i+1)];
       end

       %at the inteface
       [hs1] = Newton(q1l(1),q1r(1),q1l(2)./q1l(1),q1r(2)./q1r(1));
       us1 = q1l(2)./q1l(1) - phi(hs1,q1l(1));
       q1e = [hs1; hs1*us1]; 
       f1 = flux(q1e); %f_{i+1/2}

       if i == 1
           q2l = [h(i); hu(i)];
           q2r = [h(i); hu(i)];
       else
           q2l = [h(i-1); hu(i-1)];
           q2r = [h(i); hu(i)];
       end
        
       %at the inteface
       [hs2] = Newton(q2l(1),q2r(1),q2l(2)./q2l(1),q2r(2)./q2r(1));
       us2 = q2l(2)./q2l(1) - phi(hs2,q2l(1));
       q2e = [hs2; hs2*us2]; 
       f2 = flux(q2e); %f_{i-1/2}

       %fluctuations
        q = [h(i); hu(i)];
        ampdq = flux(q) - f2;
        apdq = f1 - flux(q);
        
        hn(i) = h(i) - (dt/dx)*(apdq(1) + ampdq(1))
        un(i) = hu(i) - (dt/dx)*(apdq(2) + ampdq(2));

    end
    %plot h
    plot(hn(1:end-1))
    %overwrite the solution
    h = hn;
    hu = un
end

% fulx
function [f] = flux(q)
    global g;
    q1 = q(1); q2 =q(2);
    f = zeros(2,1);
    f(1) = q2;
    f(2) = (((q2)^2)/q1) + (0.5*g*(q1)^2);
end

%boundary conditions
function [qext] = bc(Q)
    %extend Q with extrapolation boundary conditions
    qext = [Q([2 1],:); Q; Q([end-1 end-2],:)];
end

% Newton solver
function [hs] = Newton(hl,hr,ul,ur)
    global g;
    hs = ((sqrt(hl) + sqrt(hr) - (ur-ul)/2/sqrt(g))^2)/4;
    tol = 1e-12;
    max_iter = 100;

    for i=1:max_iter
        gk = func(hs,hl,hr,ul,ur);
        res = abs(gk);

        if (res<tol)
            break
        else
            continue
        end

        dg = dfun(hs,hl,hr,ul,ur);
        dh = -gk/dg;
        delta = 1;

        for i=1:500
            if (abs(func(hs+dh*delta,hl,hr,ul,ur)) >= res)
                delta = 0.5*delta;
            else
                break
            end
        end
        hs = hs + delta*dh;
    end
end

% phi function
function [h] = phi(hs,hlr)
    global g; 
    if (hs>hlr)
        h = sqrt(0.5*g*(hs + hlr)/(hs*hlr))*(hs - hlr);
    else
        h = 2*sqrt(g)*(sqrt(hs) - sqrt(hlr));     
    end
end

% function f
function [f] = func(hs,hl,hr,ul,ur)
    global g; 
    f = phi(hs,hl) + phi(hs,hr) + ur - ul;
end

% Jacobian of f
function [df] = dfun(hs,hl,hr,ul,ur)
    global g; 
    eps = 1e-7;
    df = (func(hs+eps,hl,hr,ul,ur) - func(hs-eps,hl,hr,ul,ur))/(2*eps);
end

