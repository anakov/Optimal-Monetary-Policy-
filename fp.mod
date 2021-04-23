// "Optimal Monetary Policy with State-Dependent Pricing",
// International Journal of Central Banking, 2014 
// By A. Nakov and C. Thomas


var    c n y w;
var    z  g cp_shock;
varexo uz ug ucp; 
 
parameters  gamma chi varphi epsilon rhoz rhog gbar rhocp    ; 
 
gamma       = 2; 
chi         = 6; 
varphi      = 1; 
epsilon     = 7; 
rhoz        = 0.95;
rhog        = 0.9;
rhocp       = 0.8;
gbar        = 0.1;
 
model; 
c^(-gamma)*w = chi*n^varphi; 
z*n = y; 
w = (epsilon-1)/epsilon*z/cp_shock;
y = c + g;
 
// Exogenous processes
z = z(-1)^rhoz*exp(uz); 
g/gbar = (g(-1)/gbar)^rhog*exp(ug);
cp_shock = cp_shock(-1)^rhocp*exp(ucp);
end; 
 
shocks; 
var uz; stderr .01; 
var ug; stderr .01; 
var ucp; stderr .01; 
end; 
 
yss = fsolve(@(y) epsilon/(epsilon-1)*chi*y^varphi*(y-gbar)^gamma - 1, .3);
initval; 
y       = yss; 
g       = gbar;
c       = y - g;    
n       = y;    
z       = 1;
w       = (epsilon-1)/epsilon;
cp_shock = 1;
end; 

steady(solve_algo=1);
stoch_simul(order=1,nograph,noprint,nomoments,nocorr,nofunctions);

yfp_uz = oo_.irfs.y_uz;
yfp_ug = oo_.irfs.y_ug;
yfp_ucp = oo_.irfs.y_ucp;

save fp yfp_uz yfp_ug yfp_ucp
