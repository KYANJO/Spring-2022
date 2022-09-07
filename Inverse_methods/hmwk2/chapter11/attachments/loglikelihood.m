function l = loglikelihood(m)

    global t; global H;
    global sigma; global D; global Q;
    
    fvec = fun(m);
    
    % The log likelihood is (-1/2)*sum(fvec(i)^2,i=1..n);
    
    l=(-1/2)*sum(fvec.^2);
end