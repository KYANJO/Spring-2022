%%====================================================================
% Boise State University
% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
% main script
%%=====================================================================

clc
clear
close all

global g; global N;
global sig; global h;
global t; global x;
global M;
warning('off','all')

% problem
hl = 2;                          % left depth
hr = 1;                          % right depth
ul = -0.5;                       % left velocity
ur = -0.5;                       % right velocity

% Spatial domain
ax = -5;
bx = 5;
ay = -2;
by = 4;
meqn = 2;                       % Number of equations in the system


g = 1;                          % Gravity 

to = 0;                         % initial time
Tfinal = 1;                     % final time

ql = [hl; ul; hl*ul];           % left conservation variable
qr = [hr; ur; hr*ur];           % right consrvation variable
qm = [(hl+hr)/2 (ul+ur)/2 ...   % intermediate initial state
    (hl*ul+hr*ur)/2];

N = 64;                         % Number of spartial steps

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
h = 0.9;

sig = 1e-7;                     % standard deviation

L0 = eye(3);          % generate roughening matrices
L1 = get_l_rough(3,1);
L2 = get_l_rough(3,2); 


qll = [hl ul hl*ul];            % left conservation variable
qrr = [hr ur hr*ur];            % right conservation variable
xi = x(1)/t(2);                 % initial speed
[Q(:,1),hs,us] = forestclaw(ql,qr,xi); % initial conservation variable
qmm = [hs us hs*us];                   % intermediate state
%mo = [qll;qm;qrr];
mtrue = [qll; qmm; qrr];               % true parameters

mo = 1.5*mtrue;
main_function(mtrue,mo,h,L0,x,ql,qr,M,t,mq) %calling the main function script


