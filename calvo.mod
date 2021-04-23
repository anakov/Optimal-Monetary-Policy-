// "Optimal Monetary Policy with State-Dependent Pricing",
// International Journal of Central Banking, 2014 
// By A. Nakov and C. Thomas


var    c n y r pi pstar w Delta K F;
var    lm1 lm2 lm3 lm4 lm5 lm6 lm7 lm8 lm9;
var    z  g cp_shock;
varexo uz ug ucp; 
 
parameters  gamma chi varphi epsilon beta rhoz rhog rhocp lambda gbar; 
 
gamma       = 2; 
chi         = 6; 
varphi      = 1; 
epsilon     = 7; 
beta        = (1/1.04)^(1/4);
rhoz        = 0.95;
rhog        = 0.9;
rhocp       = 0.8;
lambda      = 1/3; 
gbar        = 0.1;
 
model; 
c^(-gamma)*w = chi*n^varphi; 
1 = r*beta*(c(+1)/c)^(-gamma)/pi(+1); 
z*n = y*Delta; 
y = c + g;

1 = lambda*pstar^(1-epsilon) + (1-lambda)*pi^(epsilon-1); 
Delta = lambda*pstar^(-epsilon) + (1-lambda)*pi^epsilon*Delta(-1); 
pstar = K/F*cp_shock; 
K = c^(-gamma)*epsilon/(epsilon-1)*w/z*y + beta*(1-lambda)*pi(+1)^epsilon*K(+1); 
F = c^(-gamma)*y + beta*(1-lambda)*pi(+1)^(epsilon-1)*F(+1); 
 
// Policymaker's First-Order Conditions 
c^(1-gamma)/c-lm1*c^(-gamma)*gamma/c*w-lm2*r*beta*(c(+1)/c)^(-gamma)*gamma/c/pi(+1)+lm2(-1)*r(-1)*(c/c(-1))^(-gamma)*gamma/c/pi-lm4+lm8*c^(-gamma)*gamma/c*epsilon/(epsilon-1)*w/z*y+lm9*c^(-gamma)*gamma/c*y; 
-lm3*y+lm6+beta*lm6(+1)*(-1+lambda)*pi(+1)^epsilon; 
lm7*K/F^2+lm9-lm9(-1)*(1-lambda)*pi^(epsilon-1); 
-lm7/F+lm8-lm8(-1)*(1-lambda)*pi^epsilon; 
-chi*n^(1+varphi)/n-lm1*chi*n^varphi*varphi/n+lm3*z; 
lm2(-1)*r(-1)*(c/c(-1))^(-gamma)/pi^2+lm5*(-1+lambda)*pi^(epsilon-1)*(epsilon-1)/pi+lm6*(-1+lambda)*pi^epsilon*epsilon/pi*Delta(-1)-lm8(-1)*(1-lambda)*pi^epsilon*epsilon/pi*K-lm9(-1)*(1-lambda)*pi^(epsilon-1)*(epsilon-1)/pi*F; 
-lm5*lambda*pstar^(1-epsilon)*(1-epsilon)/pstar+lm6*lambda*pstar^(-epsilon)*epsilon/pstar+lm7; 
-lm2*beta*(c(+1)/c)^(-gamma)/pi(+1); 
lm1*c^(-gamma)-lm8*c^(-gamma)*epsilon/(epsilon-1)/z*y; 
-lm3*Delta+lm4-lm8*c^(-gamma)*epsilon/(epsilon-1)*w/z-lm9*c^(-gamma); 
 
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
r       = 1/beta;   
pi      = 1;    
pstar   = 1;    
w       = chi*n^varphi*c^gamma; 
Delta   = 1;    
z       = 1;
cp_shock = 1;
K       = c^(-gamma)*epsilon/(epsilon-1)*w*y/(1-beta*lambda); 
F       = c^(-gamma)*y/(1-beta*lambda); 
lm1 = 0;
lm2 = 0;
lm3 = 0;
lm4 = 0;
lm5 = 0;
lm6 = 0;
lm7 = 0;
lm8 = 0;
lm9 = 0;
end; 

steady(solve_algo=1);
stoch_simul(order=1,nograph,noprint,nomoments,nocorr,nofunctions);

save calvo.mat
%irf_calvo

delete calvo.m
delete calvo_static.m
delete calvo_dynamic.m
