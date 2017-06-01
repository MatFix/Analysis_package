function [c] = linFit(x,y,n)
% Linear Polynomial fit

    x = x(:);
    y = y(:);
    
    A = zeros(length(x), n);
    A(:,1) = 1;
    
    for ii = 1:n
        A(:,ii+1) = x.^n;
    end
    
    c = A\y;
    
    c = flip(c');
end

