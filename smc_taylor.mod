// "Optimal Monetary Policy with State-Dependent Pricing",
// International Journal of Central Banking, 2014 
// By A. Nakov and C. Thomas

var
    // We first declare the variables that are common to all the models that we are going to consider
    c           // consumption
    n           // labor hours
    y           // output
    r           // gross nominal interest rate
    pi          // gross inflation rate
    pstar       // relative optimal price
    w           // real wage
    Delta       // price dispersion
    z           // aggregate productivity process
    g           // government consumption

    // We now declare the variables that are specific to this particular model  
    @#for j in 1:24
        psi@{j}
    @#endfor
    @#for j in 0:(24-1)
        v@{j}
    @#endfor
    @#for j in 1:(24-1)
        lambda@{j}
    @#endfor
    @#for j in 1:(24-1)
        piacc@{j}
    @#endfor
    @#for j in 1:(24-1)
        thetaacc@{j}
    @#endfor

    // Average menu costs in each cohort
    @#for j in 1:(24-1)
        Xi@{j}
    @#endfor
;    

varexo
    ug          // shock to government consumption
    uz          // shock to the aggregate productivity process
    ur;

parameters
    // We first declare the parameters that are common to all the models that we are going to consider
    gamma       // coefficient of RRA
    chi         // scale parameter of labor disutility
    varphi      // convexity of labor disutility
    epsilon     // elasticity of substitution across consumption varieties
    beta        // subjective discount factor
    rhoz rhog   // autocorrelation coefficients of 'z' and 'g'
    gbar        // steady state value of 'g'
    delta       // (one minus) retirement rate
    // We now declare the parameters that are specific to this particular model
    alpha miniepsilon
    piacc0 thetaacc0 lambda24
    pibar phir phipi;


gamma       = 2;
chi         = 6;
varphi      = 0;
epsilon     = 7;
beta        = (1/1.04)^(1/4);
rhoz        = .95;
rhog        = .9;
delta       = 1;//.95;//
pibar       = (1.02)^(1/4);
phir        = .8;
phipi       = 1.5;
miniepsilon = 0.0000000001;
alpha       = .0005; 
piacc0      = 1;
thetaacc0   = 1;
lambda24    = 1;
gbar = .1;

save calibrationq gamma chi varphi epsilon beta pibar miniepsilon alpha gbar


model;
    // We first write the equations that are common to all the models that we are going to consider

        c^(-gamma)*w = chi*n^varphi;

        1 = r*beta*(c(+1)/c)^(-gamma)/pi(+1);

        n = y*Delta/z 
                       @#for j in 1:(24-1)
                           + psi@{j}*Xi@{j}/w
                       @#endfor
                        ;

        y = c + g;

        z = z(-1)^rhoz*exp(uz);

        g/gbar = (g(-1)/gbar)^rhog*exp(ug);

    // We now write the equations that are are specific to this particular model

        pstar = epsilon/(epsilon-1)*(c^(-gamma)*y*w/z
                                                    @#for j in 1:(24-1)
                                                        + beta^@{j}*c(+@{j})^(-gamma)*thetaacc@{j}(+@{j})*piacc@{j}(+@{j})^epsilon*y(+@{j})*w(+@{j})/z(+@{j})
                                                    @#endfor
                                     )/(c^(-gamma)*y	 
                                                    @#for j in 1:(24-1)
                                                        + beta^@{j}*c(+@{j})^(-gamma)*thetaacc@{j}(+@{j})*piacc@{j}(+@{j})^(epsilon-1)*y(+@{j})
                                                    @#endfor
                                        );

        @#for j in 1:(24-1)
            lambda@{j} = (miniepsilon + (v0-v@{j})/w)/(alpha + (v0-v@{j})/w); // 'lambda' corresponds to the CDF of the menu cost distribution, evaluated at each cohort's adjustment gain
        @#endfor

        thetaacc1 = 1 - lambda1;
        @#for j in 2:(24-1)
            thetaacc@{j} = (1 - lambda@{j})*thetaacc@{j-1}(-1); 
        @#endfor

        piacc1 = pi;
        @#for j in 2:(24-1)
            piacc@{j} = pi*piacc@{j-1}(-1); 
        @#endfor

        1 = pstar^(1-epsilon)*(
                               @#for j in 1:24
                                  + lambda@{j}*psi@{j} 
                               @#endfor
                               )  
                                       @#for j in 1:(24-1)
                                          + (pstar(-@{j})/piacc@{j})^(1-epsilon)*(1-lambda@{j})*psi@{j}
                                       @#endfor
                                        ;

        Delta = pstar^(-epsilon)*(
                               @#for j in 1:24
                                  + lambda@{j}*psi@{j} 
                               @#endfor
                               )  
                                       @#for j in 1:(24-1)
                                            + (pstar(-@{j})/piacc@{j})^(-epsilon)*(1-lambda@{j})*psi@{j}
                                       @#endfor
                                        ;

        psi1 = 1 
                @#for j in 2:24
                    - psi@{j} 
                @#endfor
                ;

        @#for j in 2:24
            psi@{j} = (1-lambda@{j-1}(-1))*psi@{j-1}(-1);
        @#endfor

        @#for j in 0:(24-2)
            v@{j} = (pstar(-@{j})/piacc@{j}-w/z)*(pstar(-@{j})/piacc@{j})^(-epsilon)*y + beta*(c(+1)/c)^(-gamma)*(lambda@{j+1}(+1)*v0(+1) + (1-lambda@{j+1}(+1))*v@{j+1}(+1) - Xi@{j+1}(+1));
        @#endfor
        v@{24-1} = (pstar(-@{24-1})/piacc@{24-1}-w/z)*(pstar(-@{24-1})/piacc@{24-1})^(-epsilon)*y + beta*(c(+1)/c)^(-gamma)*(delta*v0(+1) + (1-delta)*y(+1)/epsilon/(1-beta));

        @#for j in 1:(24-1)                        
            Xi@{j}/w = (alpha-miniepsilon)*log(1+(v0-v@{j})/w/alpha) - (v0-v@{j})/w*(1-lambda@{j}); //See file 'integration of menu costs.pdf' for analytical derivations of 'Xi@{j}' for different menu cost distributions
        @#endfor

    // Monetary policy rule

        r/(pibar/beta) = (r(-1)/(pibar/beta))^phir*((pi/pibar)^phipi)^(1-phir)*exp(ur);
end;


shocks;
    var ug; stderr .01; 
    var uz; stderr .01;
end;


J = 24; % Number of cohorts

lambdass = rand(1,J-1);
theta =  1-lambdass;
Xiss = zeros(J,1);

vecdiff = 1;
while vecdiff > .000001;
    
    A = zeros(J,J);
    A(1,2:J) = -1;
    for j=2:J
        A(j,j-1) = 1-lambdass(j-1);
    end;
    B = [1; zeros(J-1,1)];
    psiss = inv(eye(J)-A)*B;

    pstarss = ([lambdass 1]*psiss + sum((pibar^(epsilon-1)).^(1:1:J-1).*(1-lambdass).*psiss(1:J-1)'))^(1/(epsilon-1));

    Deltass = ([lambdass 1]*psiss + sum((pibar^epsilon).^(1:1:J-1).*(1-lambdass).*psiss(1:J-1)'))*pstarss^(-epsilon);

    cumprodtheta = cumprod(theta);
    wss = (epsilon-1)/epsilon*pstarss*(1 + sum((beta*pibar^(epsilon-1)).^(1:1:J-1).*cumprodtheta)) ...
                                 /(1 + sum((beta*pibar^epsilon).^(1:1:J-1).*cumprodtheta));

    yss = fsolve(@(yss) chi*(yss*Deltass + psiss(1:J-1)'*Xiss(1:J-1)/wss)^varphi*(yss-gbar)^gamma - wss, .5);    %(wss*Deltass^gamma/chi)^(1/(varphi+gamma));
    nss = yss*Deltass + psiss(1:J-1)'*Xiss(1:J-1)/wss;
    css = yss - gbar;

    F = (pstarss./pibar.^(0:1:J-1) - wss).*(pstarss./pibar.^(0:1:J-1)).^(-epsilon)*yss;
    A = zeros(J,J);
    A(:,1) = [lambdass'; 1];
    for j=1:J-1
        A(j,j+1) = 1-lambdass(j);
    end;
    v = (eye(J)-beta*A)\(F'-beta*Xiss);
    v0ss = v(1);
    vss = v(2:J); 
    
    for j=1:J-1
        Xiss(j) = wss*quadl(@(k) k.*(1./(alpha+k)-(miniepsilon+k)./(alpha+k).^2), 0, (v0ss-vss(j))/wss);
    end;

    for j=1:J-1
        lambdaa(j) = (miniepsilon + (v0ss-vss(j))/wss)/(alpha + (v0ss-vss(j))/wss);
    end;
    vecdiff = max(abs(lambdaa-lambdass));
    jump = .5;
    lambdass = jump*lambdaa + (1-jump)*lambdass;
    
    for j=1:J-1
        theta(j) = 1 - lambdass(j);
    end;

end;


initval;
    y       = yss;
    c       = css;
    n       = nss;
    r       = pibar/beta;
    pi      = pibar;
    pstar   = pstarss;
    w       = wss;
    Delta   = Deltass;
    z       = 1;
    g       = gbar;

    @#for j in 1:(24-1)
        lambda@{j} = lambdass(@{j}); 
    @#endfor
    
    @#for j in 1:(24-1)
        piacc@{j} = pibar^@{j}; 
    @#endfor

    thetaacc1 = 1-lambda1;
    @#for j in 2:(24-1)
        thetaacc@{j} = (1-lambda@{j})*thetaacc@{j-1};
    @#endfor

    @#for j in 1:24
        psi@{j} = psiss(@{j}); 
    @#endfor

    v0 = v0ss;
    @#for j in 1:(24-1)
        v@{j} = vss(@{j}); 
    @#endfor

    @#for j in 1:(24-1)
        Xi@{j} = Xiss(@{j});
    @#endfor
end;


steady(solve_algo=0);
stoch_simul(nograph,order=1);