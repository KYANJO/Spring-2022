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