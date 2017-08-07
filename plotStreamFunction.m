function plotStreamFunction(Z, k, W,Wz, omega, U,Uz, By, Bz, f, alpha)
x = linspace(0, 2*2*pi./k, 1000);

% for i=1:length(Z)
%    psi(:,i) = cumtrapz(x, (W(i).*exp(1i.*k.*x)));
% %    psi(:,i) = 1i.*k.*W(i).*ones(size(x));
% end
% 
% subplot(1,2,1)
% contour(x, Z, real(psi).');
% subplot(1,2,2)
% contour(x, Z, imag(psi).');

gamma = U - omega./k;
b = -1i.*(gamma.^2-k.^(-2).*By.*alpha).^(-1).*k.^(-3).*f.*By.*Wz ...
    + 1i.*(gamma.^2-k.^(-2).*By.*alpha).^(-1).*k.^(-1).*gamma.*W.*Bz;
u = 1i.*k.*Wz;
v = 1i.*k.^(-1).*gamma.^(-1).*(f.*u - alpha.*b);
for i=1:length(Z)
    bt(:,i) = real(b(i).*exp(1i.*k.*x));
    ut(:,i) = real(u(i).*exp(1i.*k.*x));
    vt(:,i) = real( v(i).*exp(1i.*k.*x));
    wt(:,i) = real(W(i).*exp(1i.*k.*x));
    wbp1(:,i) = real(b(i).*W(i).*exp(1i.*2*k.*x));
    vbp1(:,i) = real(By(i)./Bz(i).*b(i).*v(i).*exp(1i.*2*k.*x));
%     vrs(:,i) = real(Uz(i).*W(i).*u(i).*exp(1i.*2*k.*x));
end
wbp = real(bt.*wt);
vbp = repmat(By./Bz, [length(x) 1]).*real(vt.*bt);
Bzr = repmat(Bz, [length(x) 1]);
Byr = repmat(By, [length(x) 1]);
wbpc = real(bt.*(vt.*sin(alpha)+wt.*cos(alpha)));
vbpc = real(bt.*(vt.*cos(alpha)-wt.*sin(alpha))).*(Byr.*cos(alpha)-Bzr.*sin(alpha))./(Bzr.*cos(alpha)+Byr.*sin(alpha));
ind = find(real(gamma)>=0, 1, 'first');     
stlevel = Z(ind);
h = figure;
subplot(1,3,1)
contour(x.*k./(2*pi), Z-Z(1) , bt.');
title('b')
set(gca, 'FontSize', 14);
xlabel('$\hat{x}$');
ylabel('HAB (m)');

subplot(1,3,2)
contour(x.*k./(2*pi), Z-Z(1) , ut.');
title('u')
set(gca, 'FontSize', 14);
xlabel('$\hat{x}$');
% subplot(1,4,3)
% contour(x, Z, vt.');
% title('v');
ylabel('HAB (m)');

subplot(1,3,3)
plot(-mean(wbpc)./max(abs(mean(wbpc))), Z-Z(1), 'LineWidth', 2);
hold on
% plot(mean(vrs), Z);
% plot(mean(vbpc), Z);
xt =get(gca,'xtick');
plot(xt, stlevel.*ones(size(xt))-Z(1),'--');
hold off
% legend('<wb>', 'U_z<wv>')
title('$-\langle w b\rangle_c$');
set(gca, 'FontSize', 14);
grid on
ylabel('HAB (m)');

set(gcf, 'Color','w', 'Position', [675   569   851   405]);
% subplot(1,4,3)
% contour(x, Z , vt.');
% title('w''b''')
% 
% subplot(1,4,4)
% contour(x, Z , wbp.');
% title('w''b''')
% uistack(h, 'bottom');
end