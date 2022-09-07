% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
%  [Q] = model(m)
%
% INPUT
%   m - a guess at the model
%
% OUTPUT
%   Q   - a model matrix G(ql,qr,x/t) for the forward problem
%%====================================================================

function [Q] = model(m)
    global N; global M; global h;
    global t; global x;

    for j= 2:M % time loop
        for i=1:N % spartial loop
            xi = x(i)/t(j);
            [Q(:,i),hs,us] = forestclaw(m(1,:),m(3,:),xi);
        end
    end
    Q = Q'; 