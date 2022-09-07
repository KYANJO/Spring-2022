function m = occam(Gm,L,m0,d,delta)

global x;
n = length(m0);
m = m0;
oldm = zeros(size(m));
r = norm(oldm-m) + 1;
iter = 0;
mchi2 = 10;

while (r > 1e-6)

    iter = iter + 1;
    
    % store the old mode to test for convergance
    oldm = m;
    
    % get the current data that would be generated and the jacobian
    
    G = Gm(x,m);
    
    [J,FAC] = numjac(Gm,x,m,G,1e-7*eye(20,1),[]);
    
    % get the dhat that is in equation 10.14
    dhat = d - G + J * m;
    
    alphas = logspace(-12, 0, 50);
    
    for i = 1:length(alphas)
        
        M = J' * J + alphas(i)^2 * L' * L;
        
        % if M is not terribly conditioned
        if (cond(M) < 1.0e15)
            m = M\(J' * dhat);
            
            % store the associated data misfit
            chis(i) =  norm(Gm(x,m)-d)^2;
        else
            % M behaves poorly enough it should not be used
            chis(i) = + Inf;
        end
        
    end
    
    [Y, I] = min(chis);
    
    if (Y > delta^2)  % Improve Chi^2 if the condition is not satisfied
        
        alpha = alphas(I(1));
        
    else
        
        I = find(chis <= delta^2);
        alpha = alphas(max(I));
    end
    
    % store the new model and misfit
    m = (J' * J + alpha^2 * L' * L) \ (J' * dhat);
    r = norm(oldm - m);
    mchi2 = norm(Gm(x,m)-d)^2;
end