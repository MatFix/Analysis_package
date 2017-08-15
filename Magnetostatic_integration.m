phi = 0:0.0005:(2*pi);
theta = phi-0.000000001;

[PHI,THETA] = meshgrid(phi,theta);

beta = @(c,PHI) atan2(-c^2*cos(PHI),sin(PHI)) + pi/2;
c = 1:0.1:3;
for i = 1:length(c)
    k = c(i);
    ker = (cos(beta(k,PHI) - PHI).*cos(beta(k,THETA) - THETA))./abs(PHI - THETA);
    Int2(i) = trapz(phi,trapz(phi,ker,2));
end

plot(c,Int2)
hold on
