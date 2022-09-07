function lr = logproposal(x,y)

    global stepsize
    
    lr=(-1/2)*sum((x-y).^2./stepsize.^2);
end
