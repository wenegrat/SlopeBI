%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find Eigenvalues of BI on a slope
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kvecs, omeg, omegas, omegt, out, yout] = ShootEigenvalues(omegainit, alpha, nsteps, maxZ)

warnId = 'MATLAB:ode45:IntegrationTolNotMet';
warnstate = warning('error', warnId);  
plotS = false;

% maxZ= 500;
Z =  linspace(-maxZ, 0, 1000); % Slope normal 
deltaz = maxZ./length(Z);
f=1e-4;


% Ad Hoc Structure
U = 0.0001.*maxZ.*(Z./max(abs(Z)) + 1);
Uz = gradient(U, Z);
By = -f^(1).*Uz;
Ri = 100.*ones(size(Uz));
Bzf = Ri.*Uz.^2;
% Bz = (10*f).^2.*ones(size(Z));
tpoint = 1*250;
% Bz(1:tpoint) = 1.*Bz(1:tpoint)./100; %Variable N^2
% alpha = 0.001;
Bstr = 1/2.*(tanh( (Z-Z(tpoint))./40 ) - 1);
% Bstr = fliplr(Bstr);
Bz = Bzf.*10.^(2*Bstr);
delta = -(By(1).*cos(alpha) + Bz(1).*sin(alpha))./abs(By(1)+1i.*Bz(1));
Ri = Bz./Uz.^2;

% % Arrested Ekman Solution
% Uo = 0.5; % Set interior velocity scale
% S = 0.5;
% N2c = f.^2./alpha.^2.*ones(size(Z)).*S; % Set slope burger number to 1
% By = -abs(alpha).*N2c;
% Uz = -By./f;
% Bz = N2c;
% tpoint = floor(abs(Uo.*alpha./(f.*S))./deltaz);
% disp(num2str(tpoint))
% Bstr = 1/2.*(tanh( (Z-Z(tpoint))./10 ) - 1);
% Bz = Bz.*10.^(3*Bstr);
% Uz = 10.^(-1.*(1+Bstr)).*Uz;
% % By = 10.^(-5.*(1+Bstr)).*Uz;
% U = cumtrapz(Z, Uz);

Ri = Bz./Uz.^2;
% Ric = N2c./(N2c)
out.Ri = Ri;
out.Uz = Uz;
out.Bz = Bz;
out.Z = Z;

% figure
% subplot(1,3,1)
% plot(Uz, Z);
% hold on;
% plot(sqrt(Bz), Z);
% hold off
% 
% subplot(1,3,2)
% plot(Ri, Z);
% 
% subplot(1,3,3)
% plot(U, Z);

Uscale = mean(abs(U(1:tpoint)));
Uscale = abs(U(tpoint));
km = max(Uscale)./f.*sqrt(5/2/(1+Ri(2)));
lm = 2*pi*Uscale./f.*sqrt((1+Ri(2))./(5/2));
disp(['Richardson Number: ', num2str(Ri(2)), '  l = ', num2str(lm),'   alpha = ', num2str(alpha),  '  Delta = ', num2str(delta)]);
epsfinish = 5e-3;
maxI = 10;
% lscales = linspace(lm*10, lm./2, 20);
% lscales = linspace(60e3, 1e3, 100);
% lscales = linspace(lm.*100, 1
% % lscales = 6e3;
% % 
% kvecs =2*pi./lscales;
% 
% kvecs = linspace(2*pi./100e3, 2*pi./1e3, 100);
% lscales = 2*pi./kvecs;

hc = .1;
ls = linspace(1e-3, 1,nsteps);
thetb = 0;
thets = 10;
C = (1-thetb).*sinh(thets.*ls)./(sinh(thets)) + thetb.*(tanh(thets*(ls+1/2))./(2*tanh(1/2*thets)) - 1/2);
S = hc.*ls+ (1-hc).*C;
% lscales = fliplr(60e3.*S);
lscales = fliplr(100e3*S);

lscales = [fliplr(linspace(100e3, lm, nsteps./2)) (linspace(lm+1, lm.*0.5, nsteps./2))];
lscales = [linspace(100e3, 25e3, nsteps./2) linspace(25e3+1, 0.5.*lm, 1*nsteps./2)];
% lscales = linspace(1*lm, 0.5.*lm,nsteps);
kvecs = 2*pi./lscales;


figure
title(num2str(alpha))
omeg = NaN(length(kvecs),1);
omegt = omeg;

for i=1:length(kvecs);

 
% omega = (1i).*(2*pi*sqrt(5/54).*f./(sqrt(1+Ri(2)))); %Should probably initialize from flat bottom sol.
% omega = epsilon;
% if i>1
%     omega = omeg(i-1);
% end
k = kvecs(i);
kn = U./f.*k;
kn = Uscale./f.*k; %Normalize to top of BL.

for t=2;1:2;

    if t==1 % Use interpolated omega if possible
        disp('Using Interpolated Omega');
        mask = isfinite(omegt);
        if (sum(mask)>1)
            disp('Interpolating Omega');
            omega =  interp1(kvecs(mask), omegt(mask), kvecs(i), 'pchip');
%             if imag(omega)<=0; omega = real(omega); end
        else
            omega = NaN;
        end
        
    else % Try to use passed in conditions
        if (isempty(omegainit) | ~isfinite(omegainit(i)))
            disp('Using analytical Omega');
           omega = f.*(kn./2 + 1i.*kn./(2*sqrt(3)).*(1-(2*kn.^2)/15.*(1+Ri(2))));
%            omega = .001.*f.*(1+1i);
        else
            disp('Using passed Omega');
            omega = omegainit(i);
        end
    end
    
    

omegas(i) = imag(omega);
% if ~isempty(omegainit); 
%     if isfinite(omegainit(i))
%         omega = omegainit(i); 
%     end
% end


% omega = 2*pi*max(abs(U)).*k./2+1i.*2*pi*(max(abs(U))).*k./(2*sqrt(3)).*(1-(2*k.^2)/15.*(1+Ri(2)));

disp(['Step: ', num2str(i), '   L-Scale: ', num2str(lscales(i)), '    Linear Theory: ', num2str(lm)]);
opts = odeset('RelTol',1e-11,'AbsTol',1e-12);
epsilon = (1+1i).*1e-6;
efact = 1;
lang = .1;

hold on
plot(kvecs(i), imag(omega), 'o');
hold off

%% NONLINEAR SOLVER
%%%%%%%%%%%%%%%%%%%%%
% if isfinite(omega)
% % options = optimset('TolCon', 1e-3, 'ObjectiveLimit', 1e-3);
% options = optimset( 'TolX', 1e-8, 'TolFun', 1e-6, 'MaxIter', 200, 'OutputFcn', @TolMet);
% [om, wend] = fminsearch(@endSlope, [real(omega) imag(omega)],options);
% omega = om(1)+1i.*om(2);
% if (abs(wend)>epsfinish)
%     omega = NaN;
% end
% omegai = imag(omega);
% end
% disp(num2str(abs(wend)./max(abs(wend))))
% [omegar, wend] = fmincon(@endSlopeR, real(omega), [],[], [],[],[],[],[],[]);
% disp(num2str(wend))
% omega = omegar + 1i.*omegai;
% [t y] = ode45(@SlopeODE, linspace(Z(1), Z(end), 3), [0 lang], opts);
% disp(num2str(abs(y(end,1))));
% end

% %% Direct BVP SOLVER
% %%%%%%%%%%%%%%%%%%%%
% try
%     Bz = Bzf;
% %     if t==2; break; end
% sol=bvpinit(Z, @yinit, omega);
% % opt = bvpset( 'RelTol', epsfinish);
% for k=1:10
%     disp(num2str(k));
%     Bz = Bzf.*10.^(2./(11-k).*Bstr);
% sol =  bvp4c(@SlopeODE,@SlopeBC,sol);
% end
% omega = sol.parameters;
% catch
%     disp('Error')
% end


%% SHOOTER
%%%%%%%%%%%%%%%%%%%
tic
for j=1:maxI
% Main Iteration Loop
bbc = 0; %-0.01.*k.^(-2).*(U(1)-omega./k)^(-1).*f.*lang;

try
[T1, T2, T3] = DefineGIFunctions();

[t, y] = ode113(@SlopeODE, Z, [bbc lang], opts);
yn1 = abs(y(:,1))./max(abs(y(:,1)));
% yn1 = abs(real(y(:,1))./max(abs(real(y(:,1)))));

if (yn1(end) <epsfinish)
    if (i>ceil(nsteps./2) & plotS);
     plotStreamFunction(Z, kvecs(i), y(:,1).', y(:,2).', omega,U, Uz,By, Bz, f, alpha);
    plotS = false;
    drawnow
%     uistack(gca, 'bottom')
end
    break
end
omega1 = omega;
% epsilon = 0.05.*abs(omega);
omega = omega+(1).*epsilon;
omega2=omega;
[T1, T2, T3] = DefineGIFunctions();
[t, y2] = ode113(@SlopeODE, Z, [bbc lang], opts);
yn2 = abs(y2(:,1))./max(abs(y2(:,1)));
% yn2 = abs(real(y2(:,1))./max(abs(real(y2(:,1)))));

if (yn2(end) <epsfinish)
    break
end
catch ME
    disp('Integration Failed. ');
    omega = NaN;
    epsilon = 1e-6;
    break
end
disp(['Iteration:  ', num2str(j), '   Error 1: ', num2str(abs(yn1(end,1))), '   Error 2: ', num2str(abs(yn2(end,1))), '   Omega: ', num2str(imag(omega))]);

% Assume w = m*sigma + c
% w2 = real(y2(end,1));
% w1 = real(y(end,1));
% m = (w2-w1)./((1).*epsilon);
% c = 1/2.*((w1+w2)-m.*(real(omega+omega1)));
% omegar = -c./m;
% w2 = imag(y2(end,1));
% w1 = imag(y(end,1));
% m = (w2-w1)./((1).*epsilon);
% c = 1/2.*((w1+w2)-m.*(imag(omega+omega1)));
% omega = omegar -1i.*c./m;
% epsilon = epsilon/2;

w2 = (y2(end,1));
w1 = (y(end,1));

if ~isfinite(abs(w1+w2))
    omega = NaN;
        epsilon = 1e-6;

    break
end

m = (w2-w1)./(epsilon);
if ( abs(w2)./max(abs(w2)) < abs(w1)./max(abs(w1)))
% if ( abs(real(w2))./max(abs(real(w2))) < abs(real(w1))./max(abs(real(w1))))
c = (w2)-m.*(omega);
else
c = (w1)-m.*omega1;
% efact=efact./2;
end
% c = 1/2.*(w1+w2 - m*(omega+omega1));
omega = -c./m;
omega = real(omega)+1i.*abs(imag(omega));
if (abs(omega)>100*f)
    disp('Resetting Omega')
    omega = NaN;
    break
        omega = omega1.*.1;
        hcounter = hcounter + 1;
        if (hcounter > 5)
            omega = NaN;
            break % Tends to indicate a divergent calculation.
        end
%         break
%         if (i>1)
%             omega = omeg(i-1)
%         end
%         epsilon = epsilon./10;
%         m = m./10;
end
% epsilon = (omega-omega2)./abs(omega-omega2).*1e-2.*abs(omega)*efact;
deltaomr = real(omega - omega1)./real(omega);
deltaomi = imag(omega - omega1)./imag(omega);

% fact = 1e-1.*abs(omega); if ((fact < 1e-7) | ~isfinite(fact)); fact = 1e-7; end
fact = 1e-7;
% if (yn1(end)< 0.1) | yn2(end)<0.1; fact = fact.*10; end
epsilon = (1+1i).*fact;
% epsilon = (deltaomr + 1i.*deltaomi).*fact;

% efact = efact./2;
% deltaOm = -w1./m;
% epsilon = epsilon.*10^(-j/10);
% epsilon=epsilon.*exp(1i*2*pi*rand);
% omega = omega1 + deltaOm;
% epsilon = 0.001.*(omega);%*10^floor(log10(w1));
if ~isfinite(omega)
    omega = NaN;
    break
end

end
if j == maxI;
    omega = NaN;
end
toc
%% END SHOOTER

%%
omeg(i) = imag(omega);
omegt(i) = omega;
hold on
plot(k, omeg(i), 'x');
set(gca, 'ylim', [0 5e-5]);
drawnow
hold off
end

if alpha==-0.005;
   yout(:,:,i) = y; 
end

end

% Save Z structure of most unstable mode.
% [~, I] = max(omeg);
% if ~isempty(I)
%     omega = omegt(I);
%     k = kvecs(I);
%     [T1, T2, T3] = DefineGIFunctions();
%     [t, yout] = ode113(@SlopeODE, Z, [bbc lang], opts);
% else
%     yout = NaN(length(Z), 2);
% end
% toc


%% FUNCTIONS


    function yi = yinit(Z)
        yi = -maxZ.*lang./pi.*[sin(pi.*Z./maxZ);
               pi./maxZ.*cos(pi.*Z./maxZ)];
    end

    function res = SlopeBC(ya,yb, lambda)
        res = [ya(1);
               real(yb(1))
               ya(2)-lang];
    end

    function wend = endSlope(om)
        omega = om(1)+1i.*om(2);
        [T1, T2, T3] = DefineGIFunctions(omega);
        [t y] = ode45(@(x,y) SlopeODE(x, y,T1, T2, T3), Z, [0 lang], opts);
     wend = abs((y(end,1)))./max(abs((y(:,1))));
%      wend = abs(real(y(end,1)));
    end

    function stop = TolMet(x, optim, state)
        disp(optim.fval)
        if (optim.fval<epsfinish)
            stop = true;
        else
            stop = false;
        end
    end

    function wend = endSlopeR(om)
        omega = om+omegai;
        [t y] = ode45(@SlopeODE, linspace(Z(1), Z(end), 3), [0 lang], opts);
     wend = abs(y(end,1));
    end
    function wend = endSlopeI(om)
        omega = real(omega) + 1i.*om;
        [t y] = ode45(@SlopeODE, linspace(Z(1), Z(end), 3), [0 lang], opts);
     wend =imag(y(end,1));
    end


%     function wtop = EulerODE
%         deltah;
%         
%         t1 = chi; % XXX-HACK
%        t2 = chiz - Beta + f.*alpha - k^(-2).*f.*Byi.*zeta.^(-1);
%        t3 = -Betaz + gamma.*Bzi.*zeta.^(-1);
% 
%        for zi = 2:length(z)
%            v(zi) = v(zi-1)+deltah.*(y(zi-1));
%            y(zi) = y(zi-1)+ deltah.*(-t2./t1.*
%        end
% %        y( = [y(2);
% %             -t2./t1.*y(2) - t3./t1.*y(1)];
%     end
    function dydz = SlopeODE(z, y)
%         [T1, T2, T3] = DefineGIFunctions(omega);
%        gamma = U - omega./k;
%        zeta = gamma.^2 - alpha*k.^(-2).*By;
%        chi = gamma - k.^(-2).*gamma.^(-1).*f.^2.*(1+alpha.*k.^(-2).*By.*zeta.^(-1));
%        Beta = Uz + f.*alpha.*(1-k^(-2).*Bz.*zeta.^(-1));
%        chiz = gradient(chi, Z);
%        Betaz = gradient(Beta, Z);
%        
%        t1 = chi;
%        t2 = chiz - Beta + f.*alpha - k.^(-2).*f.*By.*zeta.^(-1);
%        t3 = -Betaz + gamma.*Bz.*zeta.^(-1);
% %        
%        t1 = interp1(Z, t1, z,'pchip');
%        t2 = interp1(Z, t2, z, 'pchip');
%        t3 = interp1(Z, t3, z, 'pchip');
%        chi = interp1(Z, chi, z, 'pchip');
%        Beta = interp1(Z, Beta, z, 'pchip');
%        chiz = interp1(Z, chiz, z, 'pchip');
%        Betaz = interp1(Z, Betaz, z, 'pchip');
%        gamma = interp1(Z, gamma, z, 'pchip');
%        Byi = interp1(Z, By, z, 'pchip');
%        Bzi = interp1(Z, Bz, z, 'pchip');
%        zeta = interp1(Z, zeta, z, 'pchip');
%        Uzi = interp1(Z, Uz, z, 'pchip');
% 
%        t1 = chi; % XXX-HACK
%        t2 = chiz - Beta + f.*alpha - k^(-2).*f.*Byi.*zeta.^(-1);
%        t3 = -Betaz + gamma.*Bzi.*zeta.^(-1);
       
       
         t1 = T1(z);
         t2 = T2(z);
         t3 = T3(z);
         
         if (t1==0)
             disp('here');
         end
%        if (sum(~isfinite(t2./t1+t3./t1))>0)
%            disp('here')
%        end
%        if (sum(~isfinite(y(2)+y(1)))>0)
%            disp('and here')
%        end
       dydz = [y(2);
            -t2./t1.*y(2) - t3./t1.*y(1)];
        
%         dydz = [y(2);
%              (k.^2.*gamma.^2-f^2)^(-1).*(-(2*f.^2*Uzi)./gamma.*y(2) - k.^2.*Bzi .*y(1))];
    end

    function [T1, T2, T3] = DefineGIFunctions()
       gamma = U - omega./k;
       zeta = gamma.^2 - alpha*k.^(-2).*By;
       chi = gamma - k.^(-2).*gamma.^(-1).*f.^2.*(1+alpha.*k.^(-2).*By.*zeta.^(-1));
       Beta = Uz + f.*alpha.*(1-k^(-2).*Bz.*zeta.^(-1));
       chiz = gradient(chi, Z);
       Betaz = gradient(Beta, Z);
       
%        chi = griddedInterpolant(Z, chi, 'pchip');
%        chiz = griddedInterpolant(Z, chiz, 'pchip');
%        
       t1 = chi; % XXX-HACK
       t2 = chiz - Beta + f.*alpha - k^(-2).*f.*By.*zeta.^(-1);
       t3 = -Betaz + gamma.*Bz.*zeta.^(-1);
       
       T1 = griddedInterpolant(Z, t1, 'pchip');
       T2 = griddedInterpolant(Z, t2, 'pchip');
       T3 = griddedInterpolant(Z, t3, 'pchip');
       
       
%        chi = interp1(Z, chi, z, 'pchip');
%        Beta = interp1(Z, Beta, z, 'pchip');
%        chiz = interp1(Z, chiz, z, 'pchip');
%        Betaz = interp1(Z, Betaz, z, 'pchip');
%        gamma = interp1(Z, gamma, z, 'pchip');
%        Byi = interp1(Z, By, z, 'pchip');
%        Bzi = interp1(Z, Bz, z, 'pchip');
%        zeta = interp1(Z, zeta, z, 'pchip');
%        Uzi = interp1(Z, Uz, z, 'pchip');

       
    
    end

% figure
% plot(kvecs, omeg, 'x');
% hold on
% plot(kvecs, omegas, '-');
% hold off

end