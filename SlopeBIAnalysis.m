%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find Eigenvalues of BI on a slope
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function omeg = SlopeBIAnalysis
maxZ= 1000;
Z =  linspace(-maxZ, 0, 1000); % Slope normal 
f=1e-4;
U = 0.0001.*maxZ.*(Z./max(abs(Z)) + 1);
Uz = gradient(U, Z);
By = -f^(1).*Uz;
Bz = (10*f).^2.*ones(size(Z));
alpha = 0;
Ri = Bz./Uz.^2;

0.00001;

kvecs =2*pi./4000;0.1:.1:3;
for i=1:length(kvecs);
omega = (1i).*(2*pi*sqrt(5/54).*f./(sqrt(1+Ri(2)))); %Should probably initialize from flat bottom sol.
k = kvecs(i);
solinit=bvpinit(Z, @yinit, omega);
sol =  bvp4c(@SlopeODE,@SlopeBC,solinit);
omeg(i) = sol.parameters;
end

    function yi = yinit(Z)
        yi = 0.01.*[sin(pi.*Z./maxZ);
               cos(pi.*Z./maxZ)];
    end

    function dydz = SlopeODE(z, y, omega)
       gamma = U - omega./k;
       chi = gamma.^1 -k^(-2).*gamma.^(-1).*f.^2.*(1+k^(-2).*gamma.^(-2).*By.*alpha);
       Beta = Uz + f.*alpha + k.^(-2).*gamma.^(-2).*f.*Bz.*alpha.*(1-k.^(-2).*gamma.^(-2).*By.*alpha);
       chiz = gradient(chi, Z);
       Betaz = gradient(Beta, Z);
       
       chi = interp1(Z, chi, z);
       Beta = interp1(Z, Beta, z);
       chiz = interp1(Z, chiz, z);
       Betaz = interp1(Z, Betaz, z);
       gamma = interp1(Z, gamma, z);
       Byi = interp1(Z, By, z);
       Bzi = interp1(Z, Bz, z);
       
       t1 = -chi;
       t2 = chiz+Beta-k.^(-2).*gamma.^(-2).*f.*Byi.*(1-k^(-2).*gamma.^(-2).*Byi.*alpha);
       t3 = Betaz - Bzi.*(1-k^(-2).*Byi.*alpha);
       
       dydz = [y(2);
            -t2./t1.*y(2) - t3./t1.*y(1)];
    end

    function res = SlopeBC(ya,yb,lambda)
        res = [ya(1);
               yb(1)
               ya(2)-lang];
    end
end