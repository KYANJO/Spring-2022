%%====================================================================
% Boise State University
% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
%  [Q,J] = model_jac(m)
%
% INPUT
%   m - a guess at the model
%
% OUTPUT
%   Q   - a model matrix G(ql,qr,x/t) for the forward problem
%   J   - its corresponding Jacobian
%%====================================================================

function [J] = Jacobian(m)
    global N; global M; global h;
    global t; global x;

    for j= 2:M % time loop
        J = [];
        for i=1:N % spartial loop
            xi = x(i)/t(j);
            [Q(:,i),hs,us] = forestclaw(m(1,:),m(3,:),xi);
            
            % Formulation of Jacobian Matrix
            ej = ones(3,1); 
            [Qmin(:,i),hs,us] = forestclaw(m(1,:)-h*ej,m(3,:)-h*ej,xi);
            [Qmax(:,i),hs,us] = forestclaw(m(1,:)+h*ej,m(3,:)+h*ej,xi);
            J = [J (Qmax(:,i) - Qmin(:,i))./(2*h)]; % Jacobian
        end
        Q = Q'; J = J';
    end