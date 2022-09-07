function lp = logprior(m)

    if (m(1) >= 0) && (m(1) <= 0.01) && (m(2) >= 0) && (m(2) <= 2)
      lp=0;
    else
      lp=-Inf;
    
    end
end
