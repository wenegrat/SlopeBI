function plotStreamFunctionDouble(Z, k1, k2, W1,W1z,W2, W2z, omega1, omega2, U,Uz, By, Bz, f, alpha1, alpha2)


% for i=1:length(Z)
%    psi(:,i) = cumtrapz(x, (W(i).*exp(1i.*k.*x)));
% %    psi(:,i) = 1i.*k.*W(i).*ones(size(x));
% end
% 
% subplot(1,2,1)
% contour(x, Z, real(psi).');
% subplot(1,2,2)
% contour(x, Z, imag(psi).');
h = figure;

for m=1:2
    if m==1
        alpha = alpha1; k=k1; omega=omega1; W=W1; Wz = W1z;
    else
        alpha = alpha2; k=k2; omega=omega2; W = W2; Wz = W2z;
    end
    l = 2*pi./k; if m==1; l1=l; else l2 = l; end
    x = linspace(0, 2*2*pi./k, 1000);
gamma = U - omega./k;
b = -1i.*(gamma.^2-k.^(-2).*By.*alpha).^(-1).*k.^(-3).*f.*By.*Wz ...
    + 1i.*(gamma.^2-k.^(-2).*By.*alpha).^(-1).*k.^(-1).*gamma.*W.*Bz;
u = 1i.*k.^(-1)*Wz;
v = 1i.*k.^(-1).*gamma.^(-1).*(f.*u - alpha.*b);
for i=1:length(Z)
    bt(:,i) = real(b(i).*exp(1i.*k.*x));
    ut(:,i) = real(u(i).*exp(1i.*k.*x));
    vt(:,i) = real(v(i).*exp(1i.*k.*x));
    wt(:,i) = real(W(i).*exp(1i.*k.*x));
    wbp1(:,i) = real(b(i).*W(i).*exp(1i.*2*k.*x));
    vbp1(:,i) = real(By(i)./Bz(i).*b(i).*v(i).*exp(1i.*2*k.*x));
%     vrs(:,i) = real(Uz(i).*W(i).*u(i).*exp(1i.*2*k.*x));
end
ke = (ut.^2 + vt.^2 + wt.^2)./2;

wbp = real(bt.*wt);
vbp = repmat(By./Bz, [length(x) 1]).*real(vt.*bt);
Bzr = repmat(Bz, [length(x) 1]);
Byr = repmat(By, [length(x) 1]);
wbpc = real(bt.*(vt.*sin(alpha)+wt.*cos(alpha)));
kekm = mean(real(ut.*(vt.*sin(alpha)+wt.*cos(alpha)))).*Uz;
vbpc = real(bt.*(vt.*cos(alpha)-wt.*sin(alpha))).*(Byr.*cos(alpha)-Bzr.*sin(alpha))./(Bzr.*cos(alpha)+Byr.*sin(alpha));
if m==1;
    wbcn1 = mean(wbpc)./mean(mean(ke));
else
    wbcn2 = mean(wbpc)./mean(mean(ke));
end

ind = find(real(gamma)>=0, 1, 'first');     
stlevel = Z(ind);
gap = [0.1 0.075]; margh = .075; margw = 0.075;

subtightplot(4,3,[1 4] +6*(m-1), gap, margh, margw)
contour(x.*k./(2*pi), Z-Z(1) , bt.');
title(['b     ($\lambda=$', num2str(l./1000, 2),' km)'])
set(gca, 'FontSize', 14);
xlabel('$\hat{x}$');
ylabel('HAB (m)');
% t= text(.35, 930,0, ['$\lambda = $', num2str(l./1000, 2), ' km'], 'FontSize', 12, 'Color', 'k', ...
%     'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');


subtightplot(4,3,[2 5] +6*(m-1), gap, margh, margw)
contour(x.*k./(2*pi), Z-Z(1) , ut.');
title(['u     ($\lambda=$ ', num2str(l./1000, 2),' km)'])
set(gca, 'FontSize', 14);
xlabel('$\hat{x}$');
% subplot(1,4,3)
% contour(x, Z, vt.');
% title('v');
ylabel('HAB (m)');
% t= text(.35, 930,0, ['$\lambda = $', num2str(l./1000, 2), ' km'], 'FontSize', 12, 'Color', 'k', ...
%     'HorizontalAlignment', 'center', 'BackgroundColor','w', 'EdgeColor', 'k');
% if m==2
% subtightplot(4,3,[6 9], gap, margh, margw)
% plot(-mean(wbpc)./max(abs(mean(wbpc))), Z-Z(1), 'LineWidth', 2);
% hold on
% % plot(mean(vrs), Z);
% % plot(mean(vbpc), Z);
% xt =get(gca,'xtick');
% % plot(xt, stlevel.*ones(size(xt))-Z(1),'--');
% hold off
% % legend('<wb>', 'U_z<wv>')
% title('$-\langle w b\rangle_c$');
% set(gca, 'FontSize', 14);
% grid on
% ylabel('HAB (m)');
% xlabel('');
% end

set(gcf, 'Color','w', 'Position', [ 487          87        1067         823]);

end

subtightplot(4,3,[6 9], gap, margh, margw)
plot(wbcn2./f, Z-Z(1), 'LineWidth', 2);
hold on
plot(wbcn1./f, Z-Z(1), 'LineWidth', 2);

% plot(mean(vrs), Z);
% plot(mean(vbpc), Z);
xt =get(gca,'xtick');
% plot(xt, stlevel.*ones(size(xt))-Z(1),'--');
hold off
% legend('<wb>', 'U_z<wv>')
title('$\widehat{\langle w b\rangle_c}$');
set(gca, 'FontSize', 14);
grid on
ylabel('HAB (m)');
set(gca, 'xlim',[-1 3])
% xlabel('$s^{-1}$');
legend(['$\lambda =$', num2str(l2./1000,2),' km'], ['$\lambda =$', num2str(l1./1000,2),' km'], 'Location', 'NorthEast')
% subplot(2,3,[3 6])
% % plot(-wbcn1, Z-Z(1), 'LineWidth', 2);
% hold on
% plot(-wbcn2, Z-Z(1), 'LineWidth', 2);
% % plot(mean(vrs), Z);
% % plot(mean(vbpc), Z);
% xt =get(gca,'xtick');
% % plot(xt, stlevel.*ones(size(xt))-Z(1),'--');
% hold off
% % legend('<wb>', 'U_z<wv>')
% title('$-\langle w b\rangle_c$');
% set(gca, 'FontSize', 14);
% grid on
% ylabel('HAB (m)');
% subplot(1,4,3)
% contour(x, Z , vt.');
% title('w''b''')
% 
% subplot(1,4,4)
% contour(x, Z , wbp.');
% title('w''b''')
% uistack(h, 'bottom');
end