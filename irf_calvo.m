% irfplot;

load calvo.mat

gamma  = M_.params(1);
chi    = M_.params(2);
varphi = M_.params(3);
gbar   = M_.params(9);

y      = oo_.steady_state(3);    
c      = y - gbar;
n      = y;    
w      = chi*n^varphi*c^gamma; 
    
figure(3)
subplot(3,3,1)
plot(100*[oo_.irfs.z_uz' oo_.irfs.g_ug'/gbar])
title('shocks')
legend('z','g')

subplot(3,3,2)
plot(100*[oo_.irfs.pi_uz' oo_.irfs.pi_ug'])
title('pi')
ylim([-1 1])

subplot(3,3,3)
plot(100*[oo_.irfs.r_uz' oo_.irfs.r_ug'])
title('R')
ylim([-0.25 0.25])

subplot(3,3,4)
plot(100*[oo_.irfs.pstar_uz oo_.irfs.pstar_ug])
title('pstar')
ylim([-1 1])

subplot(3,3,5)
plot(100*[oo_.irfs.Delta_uz' oo_.irfs.Delta_ug'])
title('Delta')
ylim([-1 1])

subplot(3,3,6)
plot(100*[oo_.irfs.r_uz(1:end-1)'-oo_.irfs.pi_uz(2:end)' oo_.irfs.r_ug(1:end-1)'-oo_.irfs.pi_ug(2:end)'])
title('rr')
ylim([-0.25 0.25])

subplot(3,3,7)
plot(100*[oo_.irfs.c_uz' oo_.irfs.c_ug']/c)
title('C')
ylim([-0.5 1])

subplot(3,3,8)
plot(100*[oo_.irfs.n_uz' oo_.irfs.n_ug']/n)
title('N')
ylim([-1 1])

subplot(3,3,9)
plot(100*[oo_.irfs.w_uz' oo_.irfs.w_ug']/w)
title('w')
ylim([-0.5 1.5])


figure(2)
subplot(3,3,1)
plot(100*oo_.irfs.cp_shock_ucp')
title('Cost push shock')
legend boxoff

subplot(3,3,8)
plot(100*oo_.irfs.pi_ucp)
title('Inflation')
xlabel('Quarters')

subplot(3,3,3)
plot(100*oo_.irfs.r_ucp)
title('Nominal interest rate')

subplot(3,3,7)
plot(100*oo_.irfs.pstar_ucp)
title('Optimal price')
xlabel('Quarters')

subplot(3,3,9)
plot(100*oo_.irfs.Delta_ucp)
title('Price dispersion')
xlabel('Quarters')

subplot(3,3,2)
plot(100*(oo_.irfs.r_ucp(1:end-1)-oo_.irfs.pi_ucp(2:end)))
title('Real interest rate')

subplot(3,3,4)
plot(100*oo_.irfs.c_ucp/c)
title('Consumption')

subplot(3,3,5)
plot(100*oo_.irfs.n_ucp/n)
title('Hours worked')

subplot(3,3,6)
plot(100*oo_.irfs.w_ucp/w)
title('Real wage')



% figure
% subplot(3,1,1)
% plot(100*pi_ucp)
% title('pi')
% subplot(3,1,2)
% plot(100*pstar_ucp)
% title('pstar')
% subplot(3,1,3)
% plot(100*Delta_ucp)
% title('delta')
