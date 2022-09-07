%%====================================================================
% Boise State University
% Author: Brian KYANJO
% supervised by: Prof. Jodi Mead
% class: Inverse Methods
% Date: March 17th, 2022
% final project
%
% Description:
% -----------
% Uses the Initial Riemann problem to find an intemediate state (qm) which either 
% the intial left or right state connects to it via any combination of shocks and
% rarefactions in the two families.
% 
% For Riemann problems with an initial dry state on one side, the exact Riemann 
% solution contains only a single rarefaction connecting the wet to dry state.
% The evolving wet dry interface is therefore simply one edge of the rarefaction. 
% The propagation speed of this interface can be exactly determined using the 
% Riemann invariants of the corresponding characteristics field.
% 
% Input:
% -----
% x  - array of spacial points 
% t  - array of temporal points
% mq - specifies the output (0 and 1 corresponds to h and hu respectively) 
% ql - left initial state
% qr - right initial state
% Returns:
% h  - array of hieght field values
% hu - array of momentum field values
%%=====================================================================

% Exact solver
function [Q,hs,us] = forestclaw(ql,qr,xi)
    global g; 
    hl = ql(1); hr = qr(1);
    ul = ql(2); ur = qr(2);

    hs = Newton(hl,hr,ul,ur);  %calling the newton solver
    us = ul - phi(hs,hl);
    
    if xi <= us
        if hs > hl
            s = ul - sqrt(0.5*g*hs/hl*(hl+hs));
            if xi <= s
                h = hl;
                hu = hl*ul;
            else
                h = hs;
                hu = hs*us;
            end
        else
            head = ul -sqrt(g*hl);
            tail = us - sqrt(g*hs);
            if xi <= head
                h = hl;
                hu = hl*ul;
            elseif xi >= tail
                h = hs;
                hu = hs*us;
            else
                h = (((ul + 2*sqrt(g*hl) - xi)/3)^2)/g;
                u = xi + sqrt(g*h);
                hu = h*u;
            end
        end
    else
        if hs > hr
            s = ur + sqrt(0.5*g*hs/hr*(hs+hr));
            if xi <=s
                h = hs;
                hu = hs*us;
            else
                h = hr;
                hu = hr*ur;
            end
        else
            head = ur + sqrt(g*hr);
            tail = us + sqrt(g*hs);
            if xi >= head
                h = hr;
                hu = hr*ur;
            elseif xi <= tail
                h = hr;
                hu = hr*ur;
            else
                h = (((xi-ur+2*sqrt(g*hr))/3)^2)/g;
                u = xi - sqrt(g*h);
                hu = h*u;
            end
        end
    end
    Q = [h;hu./h;hu];
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

