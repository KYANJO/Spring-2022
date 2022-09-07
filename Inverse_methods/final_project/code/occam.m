% Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% m=occam(fun,jac,L,d,m0,delta)
%
% INPUT
%   fun   - a function handle that computes the forward problem
%   jac   - a function handle that computes the Jacobian of the forward problem
%   L     - a regularization matrix
%   d     - the data that should be fit
%   m0    - a guess at the model
%   delta - the cutoff to use for the discrepancy principle portion
%
% OUTPUT
%   m - the model found
%
function m = occam(fun, jac, L, d, m0, delta)
m = m0;
oldm = zeros(size(m));
iter = 0;
mchi2 = 100;
nr = 1;

% while we have not converged sufficiently or the data misfit is higher than
% allowed keep iterating
while (nr > 1e-6 || (mchi2 > delta^2 * 1.01))
    % only allow 30 iterations
    iter = iter + 1;
    if (iter > 30)
        return;
    end
    
    % store the old mode to test for convergance
    oldm = m;
    
    % get the current data that would be generated and the jacobian
    Gm = feval(fun, m); 
    J = feval(jac ,m);
  
    % get the dhat that is in equation 10.14
    dhat = d - Gm + J * m; 
   
    % This is a simple brute force way to do the line search.  Much more
    % sophisticated methods are available.  Note: we've restricted the line
    % search to the range from 1.0e-20 to 1.  This seems to work well in
    % practice, but might need to be adjusted for a particular problem.
    
    alphas = logspace(-12,0, 20);
    for i = 1:length(alphas)
       
        M = J' * J + alphas(i)^2 * L' * L; 
        
        % if M is not terribly conditioned
        if (cond(M) < 1.0e15)
            m = M\(J'* dhat);
            
            % store the associated data misfit
            chis(i) = norm(feval(fun, m)-d, 2)^2;
           
        else
            % M behaves poorly enough it should not be used
            chis(i) = + Inf;
        end
    end
    
    [Y, I] = min(chis);
    
    if (Y > delta^2)
        %disp('Improving Chi^2');
        alpha = alphas(I(1));
    else
        %disp('Smoothing m');
        I = find(chis <= delta^2);
        alpha = alphas(max(I));
    end
    
    % store the new model and misfit
    m = (J' * J + alpha^2 * L' * L) \(J'* dhat); 
    mchi2 = norm(feval(fun,m)-d, 2)^2;
    nr = norm(oldm - m);
    mchi2 = norm(fun(m)-d, 2)^2;
end