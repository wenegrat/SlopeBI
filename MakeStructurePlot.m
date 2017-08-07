% Make Structure Plot
% omega = % pick an omega
% alpha = 0.1;
% lang = .1;
aind = 3;
[~, index] = max(omeg(:,aind));
% 
% [T1, T2, T3] = DefineGIFunctions();
% 
% [t, y] = ode113(@SlopeODE, Z, [bbc lang], opts);
% yn1 = abs(y(:,1))./max(abs(y(:,1)));
% yn1 = abs(real(y(:,1))./max(abs(real(y(:,1)))));
y = yout(:,:,aind);
alpha = alphavec(aind);
index2 = index;
alpha2 = -0.1;
nsteps = 250;


%%
[kvecp, omeg, omegas, omegat, outstructp, yout] = ShootEigenvalues(omeginit, alphavec(aind), nsteps,nz);

%%
index1 = 131;
% index1 = 132;
index2 = 246;
y1 = squeeze(yout(:,:,index1));
y2 = squeeze(yout(:,:,index2));
alpha1 = -0.005;
alpha2 = alpha1;
plotStreamFunctionDouble(Z, kvecp(index1),kvecp(index2), y1(:,1).', y1(:,2).', y2(:,1).', y2(:,2).', omegat(index1), omegat(index2),U, Uz,By, Bz, f, alpha, alpha2);
disp(num2str(2*pi./kvecs(index)));