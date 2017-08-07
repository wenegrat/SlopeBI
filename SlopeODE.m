 function dydz = SlopeODE(z, y, omega)
       gamma = U - omega./k;
       chi = gamma.^2 -k^(-2).*gamma.^(-1).*f.^2.*(1+k^(-2).*gamma.^(-2).*By.*alpha);
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
