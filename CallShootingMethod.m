alpha = 0;
nz = 1000;
alphavec = [0 -0.0025 -0.005 -0.01];
alphavec = fliplr(-0.01:0.0025:-0.00); % For particular k.
% alphavec = [-0.01 -0.02 -0.03];
% alphavec = -0.01;
omeginit = [];
nsteps = 500;
omeg = NaN(nsteps, length(alphavec));
omegat = omeg;
kvecs = NaN(nsteps, length(alphavec));
yout = NaN(nz, 2, length(alphavec));
tic
parfor i=1:length(alphavec)
[kvecsp, omeg(:,i), omegas, omegat(:,i), outstructp, yout(:,:,i)] = ShootEigenvalues(omeginit, alphavec(i), nsteps,nz);
% mask = isfinite(omegat);
kvecs(:,i) = kvecsp;
% outstruct = outstructp
% omeginit = interp1(kvecs(mask), omegat(mask), kvecs, 'pchip');
end
toc
%% GENERATE RI
maxZ= 1000;
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

%%
gap = [0.15 0.05];margh=.12; margw=.2;

om = omeg;
om(om==0) = NaN;
omso = NaN(length(kvecs), length(alphavec));
for i=1:length(alphavec)
    try
[ks, I] = sort(kvecs(:,i));
oms = om(I,i);
mask = isfinite(oms);

omso(:,i) = interp1(ks(mask), oms(mask), ks, 'pchip', NaN);
    catch
        disp('error')
    end
end
figure
% subtightplot(2, 12, 13:14, gap, margh,margw)
% plot(Uz.^2, Z+1000, 'LineWidth', 3)
% hold on;
% plot(Bz, Z+1000, 'LineWidth', 3);
% hold off
% xlabel('Ri'); ylabel('HAB (m)');
% grid on
% % set(gca, 'xlim', [1 2e2])
% set(gca, 'FontSize', 16)

subtightplot(2, 12, 13:15, gap, margh,margw)
semilogx(Ri, Z+1000, 'LineWidth', 3)
xlabel('Ri'); ylabel('HAB (m)');
grid on
set(gca, 'xlim', [1 2e2])
set(gca, 'FontSize', 16)

subtightplot(2,12,17:24, gap, margh, margw)
hold on;
for i=1:1:length(alphavec)
    try
   plot(ks.*0.1./1e-4, smooth(omso(:,i)./1e-4, 1),'-', 'lineWidth', 3);
    catch
    end
end
set(gca, 'xlim', [0.05 0.3], 'ylim', [0 0.04]);
hold off

% legend(num2str(alphavec.'));
grid on
ylabel('$\omega_i/f$'); 
xlabel('$\hat{k}$');
set(gca, 'FontSize', 16)

subtightplot(2, 12, 1:12, gap, margh,margw)
hold on;
for i=1:1:length(alphavec)
    
   plot(ks.*0.1./1e-4, smooth(omso(:,i)./1e-4, 1),'-', 'lineWidth', 3);
end
rectangle('Position', [0.05 0 0.3 0.04], 'LineStyle', '--');
set(gca,'xlim', [0 7])
% plot(kvecs, omegas, '--');
hold off
legend(num2str(alphavec(1:1:end).'));
grid on
ylabel('$\omega_i/f$'); 
xlabel('$\hat{k}$');
set(gca, 'FontSize', 16)
set(gcf, 'Color','w', 'Position', [   675   342   919   632]);
% set(gca, 'ylim', [0 2.5e-5]);

%%
% Plot growth rate as a function of slope
ls = 3.5e3;
ks = 2*pi./ls;
ind = find(kvecs>=ks, 1, 'first');
omn = om(ind,:)./om(ind,1);
mask = isfinite(omn);
omn = interp1(alphavec(mask), omn(mask), alphavec);
figure
plot(abs(alphavec), omn, '-k', 'LineWidth', 2)
xlabel('|\alpha|');
ylabel('\omega_i/\omega_0');
title(['Length Scale: ', num2str(ls), '   k = ', num2str(kvecs(ind))]);
set(gcf, 'color', 'w');
set(gca, 'FontSize', 18);
grid on