function [k2_m,k1_m,k3_m] = kf_func(B,mx,h,display)

switch nargin
    case 2
        h = 30e-9;
        display = 'on';
    case 3
        display = 'on';
end
        %% Hysteresis cycle k determination

        Ms = 480e3;
        
        %method 1
        
        k1 = 2/3*pi*h*Ms*B./mx;
        
        %method 2
        
        r = zeros(size(B));
        
        options = optimset('Display','off');
        
        for ii=1:length(mx)
            m = mx(ii);
            r(ii) = fsolve(@(r) r^3 /8 - r + m,0.5,options);
        end
        
        k2 = 2/3*pi*h*Ms*B.*(1./r - 3/8*r);
        k2(k2==0) = NaN;
        
        %method 3
        k3 = 2/3*pi*h*Ms*B.*(1./mx - 3/8*mx);
        
        if strcmp(display,'on')
            figure
            scatter(B,k1,'filled')
            hold on
            scatter(B,k2,'filled')
            scatter(B,k3,'filled')
            xlabel('B [T]')
            ylabel('k_F [J \cdot m^{-2}]')
            ylim([0 1.2*max(k1)])
            legend('Method 1','Method 2','Method 3')
        end
        
        k1_m = mean(k1);
        k2_m = mean(k2(2:end));
        k3_m = mean(k3);
end