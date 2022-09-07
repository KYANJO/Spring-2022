function I = simpsons(f,a,b,n,p)

if(nargin == 4)
    if numel(f)>1 % If the input provided is a vector
        h=(b-a)/(n-1);
        I = h/3*(f(1)+2*sum(f(3:2:end-2))+4*sum(f(2:2:end))+f(end));
    else % If the input provided is an anonymous function
        h = (b-a)/(n-1); xi = (a:h:b)';
        I = h/3*(f(xi(1))+2*sum(f(xi(3:2:end-2)))+4*sum(f(xi(2:2:end)))+f(xi(end)));
    end
else
    h = (b-a)/(n-1);
    xi = (a:h:b)';
    
    I = h/3*(f(xi(1),p(1))+2*sum(f(xi(3:2:end-2),p(3:2:end-2)))+4*sum(f(xi(2:2:end),p(2:2:end)))+f(xi(end),p(end)));
end
