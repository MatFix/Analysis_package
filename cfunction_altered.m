function ccoeff = cfunction_altered(x,y,ux,uy)

global xx yy Mx My rCoeff

xpts = x;
ypts = y;


% rescaled vectors initialization
m_resc_x = zeros(1,length(xpts));
m_resc_y = m_resc_x;

% magnetization at mesh centers coordinates
for ii = 1: length(xpts)
    
    tmp = abs(xpts(ii) - xx);
    [~, idx] = min(tmp);
    
    tmp = abs(ypts(ii) - yy);
    [~, idy] = min(tmp);
    
    m_resc_x(ii) = Mx(idx,idy);
    m_resc_y(ii) = My(idx,idy);
end

% the |delta u| is required to rescale the gradient to a unit vector
modulus_squared = ux.^2 + uy.^2;

% this avoids problems with derivatives equal to zero - in this case 
% |delta u| gives no contribution so the value to which we fix modulus is
% arbitrary

if mean(modulus_squared) == 0
    modulus_squared = 1;
end

% this calculates the c coefficient in the elliptic pde, i.e. the local 
% conductivity (in units of sigma_0, since the pde is homogeneous anyways)

ccoeff = 1./( 1 + rCoeff./modulus_squared.*(ux.*m_resc_x + uy.*m_resc_y).^2);
end