% phi function
function [h] = phi(hs,hlr)
    global g; 
    if (hs>hlr)
        h = sqrt(0.5*g*(hs + hlr)/(hs*hlr))*(hs - hlr);
    else
        h = 2*sqrt(g)*(sqrt(hs) - sqrt(hlr));     
    end
end
