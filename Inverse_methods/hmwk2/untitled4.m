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

ql = [hl; ul];            % left conservation variable
qr = [hr; ur];           % right consrvation variable

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
Q = zeros(N,2);
for i=1:N
    if (x(i)<=0) 
        Q(i,:) = ql;
    else
        Q(i,:) = qr;
    end
end

q0 = Q;

% writing a video
v = VideoWriter('dam.avi');
open(v);

%inverse solver
[Qinv,xc,tvec] = inverse_solve(ax,bx,N,Tfinal,M,ql,qr,q0,meqn);

%intialize the inverse solver
q0 = Qinv(1,1,:);

for j= 2:M
    inver = Qinv(j,1,:);
    %plot(inver)
    for i=1:N
        xi = x(i)/t(j);
        %Q(:,i) = forward_solver(ql,qr,xi);
    end
    
%     if mq == 1
%         plot(x,Q(mq,:),'k-'); 
%         title('Height field')
%         ylim([hr-0.25,hl+0.25]);
%         ylabel('height');xlabel('x')
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     elseif mq == 2
%         plot(x,Q(mq,:),'k-'); 
%         title('velocity field')
%         ylabel('velocity');xlabel('x')
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     elseif mq == 3
%         plot(x,Q(mq,:),'k-'); 
%         title('momentum field')
%         ylabel('Momentum');xlabel('x')
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     else
%         subplot(2,2,1)
%         plot(x,Q(1,:),'k-'); title('Height field')
%         ylim([hr-0.25,hl+0.25])
%         ylabel('height');xlabel('x')
%         subplot(2,2,2);
%         plot(x,Q(2,:),'r-'); title('Velocity field')
%         ylabel('velocity');xlabel('x')
%         subplot(2,2,[3,4]);
%         plot(x,Q(3,:),'b-'); title('Momentum field')
%         ylabel('Momentum');xlabel('x')
%         %legend('Height','velocity','momentum')
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     end
%     
end
close(v);


% Inverse_solver
function [Q,xc,tvec] = inverse_solve(ax,bx,mx,Tfinal,nout,ql,qr,q0,meqn)
    
    %intial height fields(left and right)
    hl = ql(1); ul = ql(2);
    hr = qr(1); ur = qr(2);

    dx = (bx-ax)/mx;
    xe = linspace(ax,bx,mx+1); % edge locations
    xc = xe(:,end-1) + dx/2; % cell-center locations

    %temporal mesh
    t0 = 0;
    tvec = linspace(t0,Tfinal,nout+1);
    dt = Tfinal/nout;

    %store time solutions
    %Q = zeros(mx,meqn,nout+1);
    Q = zeros(nout+1,2,mx);
    %q0 = qinit(xc,meqn,ql,qr)
   
    Q(1,:,:) = q0';

    q = q0;

    dtdx = dt/dx;
    t = t0;

    for n=1:nout-1
        t = tvec(n);
        
        % Add 2 ghost cells at each end of the domain
        qext = bc(q);

        [waves,speeds,amdq,apdq] = rp(qext,meqn,xc,dx);

        %first order update
        q = q - dtdx*(apdq(2:end-2,:) + amdq(3:end-1,:));

        Q(n+1,:,:) = q';
    end

end

% Inverse Riemann solver
%function [fwaves,speeds,ampdq,apdq] = rp(Q,meqn,x,dx)
function [Q,xc,tvec] = rp(Q,meqn,x,dx)
    global g;

    %jump in Q at each interface
    delta = Q(2:end,:) - Q(1:end-1,:); 

    %d0 = delta(:,(1)); d1 = delta(:,(1));
    h = Q(:,1);
    u = Q(:,2);

    %[m,n] = size(delta);

    [mx,my] = size(Q);

    % Array of wave 1 and 2
    z1 = zeros(mx,my); 
    z2 = zeros(mx,my);
    
    % Array of speed 1 and 2
    s1 = zeros(mx,1);
    s2 = zeros(mx,1);

    for i=2:mx-1
        % Roe-averaged values
        uhat = (sqrt(h(i-1))*u(i-1) + sqrt(h(i))*u(i))/(sqrt(h(i-1))+sqrt(h(i)));
        hbar = (1/2)*(h(i-1) + h(i));
        chat = sqrt(g*hbar);

        % Eigen values
        l1 = uhat - chat; l2 = uhat + chat;

        % Eigen vectors
        r1 = [1 l1]; r2 = [1 l2];

        % Matrix of eigenvalues
        R = [r1;r2]';

        %vector of eigenvalues
        evals = [l1 l2];

        %left and right states based on backward differencing
        ql = [h(i-1); h(i-1)*u(i-1)];
        qr = [h(i); h(i)*u(i)];

        %flux at each boundary
        fl = reshape(flux(ql),[2,1]);
        fr = reshape(flux(qr),[2,1]);
%%%%%%%%%%%%%%
        q = [h(i); u(i)];
        ampdq = flux(q) - fl;
        apdq = fr - flux(q);

        
%%%%%%%%%%%%%
%         %flux difference
%         delta_flux = fr - fl - dx*psi(x,hbar);
%         d0 = delta_flux(1,:);
%         d1 = delta_flux(2,:);
% 
%         %beta_i+1/2
%         b1 = ((uhat + chat)*d0 - d1)/(2*chat);
%         b2 = (-(uhat - chat)*d0 + d1)/(2*chat);
% 
%         %f-wave and speed 1
%         z1(i-1,:) = b1.*R(:,(1))'; %%
%         s1(i-1) = evals(1);
% 
%         %f-wave and speed 2
%         z2(i-1,:) = b2.*R(:,(2))';
%         s2(i-1) = evals(2);
% 
%         %pth wave at each interface
%         fwaves = [z1;z2];
%         speeds = [s1;s2]; 
% 
%         %fluctuations
%         [m,n] = size(delta);
%         ampdq = zeros(m,n);
%         apdq = zeros(m,n);
%         for mw = 1:meqn
%             if speeds(mw) <= 0
%                 sm = 1;
%             else
%                 sm = 0;
%             end
%             ampdq = ampdq + sm.*fwaves;
% 
%             if speeds(mw) > 0
%                 sp = 1;
%                 %apdq = apdq + sp*fwaves(mw,:);
%             else
%                 sp = 0;
%                 %apdq = apdq + sp*fwaves(mw,:);
%             end
%             apdq = apdq + sp.*fwaves(:,mw);
%         end

    end

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

function [dB] = dBdx(x)
    if abs(x-0.5) < 0.1
        dB = (0.25*pi/0.1)*(sin((pi*(x - 0.5)/0.1) + 1));
    else
        dB = 0;
    end
end

function [arr] = psi(x,hm)
    global g;
    for i=1:length(x)
        a = g*hm*dBdx(x(i));
        arr = [0 a];
    end
end

function [q0] = h_init(x,hl,hr)
    if x<0
        q0 = hl
    else
        q0 = hr;
    end
end

function [q0] = hu_init(x,hl,ul,hr,ur)
    if x<0
        q0 = hl*ul;
    else
        q0 = hr*ur;
    end
end

function [q] = qinit(x,meqn,ql,qr)
    hl = ql(1);
    ul = ql(2);
    hr = qr(1); ur = qr(2);
    
    q = zeros(length(x),meqn);
    q(:,1) = h_init(x,hl,hr)
    q(:,2) = hu_init(x,hl,ul,hr,ur);
end