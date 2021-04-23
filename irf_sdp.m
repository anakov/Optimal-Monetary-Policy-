% irfplot;

load sdp_results.mat
load fp.mat

gamma  = M_.params(1);
chi    = M_.params(2);
varphi = M_.params(3);
gbar   = M_.params(8);

y      = oo_.steady_state(3);
c      = y - gbar;    
n      = y;    
w      = chi*n^varphi*c^gamma; 
    
figure(1)
subplot(3,3,1)
plot(100*[oo_.irfs.z_uz' oo_.irfs.g_ug'/gbar])
title('Shocks')
legend('Tech. process','Govt spending')
legend boxoff

subplot(3,3,8)
plot(100*[oo_.irfs.pi_uz' oo_.irfs.pi_ug'])
title('Inflation')
ylim([-1 1])
xlabel('Quarters')

subplot(3,3,3)
plot(100*[oo_.irfs.y_uz'-yfp_uz' oo_.irfs.y_ug'-yfp_ug'])
title('Output gap')
ylim([-1 1])

subplot(3,3,7)
plot(100*[oo_.irfs.pstar_uz' oo_.irfs.pstar_ug'])
title('Optimal price')
ylim([-1 1])
xlabel('Quarters')

subplot(3,3,9)
plot(100*[oo_.irfs.Delta_uz' oo_.irfs.Delta_ug'])
title('Price dispersion')
ylim([-1 1])
xlabel('Quarters')

subplot(3,3,2)
plot(100*[oo_.irfs.r_uz(1:end-1)'-oo_.irfs.pi_uz(2:end)' oo_.irfs.r_ug(1:end-1)'-oo_.irfs.pi_ug(2:end)'])
title('Real interest rate')
ylim([-0.1 0.05])

subplot(3,3,4)
plot(100*[oo_.irfs.c_uz' oo_.irfs.c_ug']/c)
title('Consumption')
ylim([-0.5 1])

subplot(3,3,5)
plot(100*[oo_.irfs.n_uz' oo_.irfs.n_ug']/n)
title('Hours worked')
ylim([-0.5 0.5])


subplot(3,3,6)
plot(100*[oo_.irfs.w_uz' oo_.irfs.w_ug']/w)
title('Real wage')
ylim([-0.5 1.5])
xlabel('Quarters')
