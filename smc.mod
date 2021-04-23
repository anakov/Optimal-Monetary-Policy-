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
    cp_shock    // cost-push factor

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

    // Density function of menu costs evalated at each cohort's adjustment gain
    @#for j in 1:(24-1)
        gL@{j}
    @#endfor

    // Lagrange multipliers of the Ramsey problem
    lm_w lm_pstar lm_n lm_pi lm_Delta
    @#for j in 1:(24-1)
        lm_lambda@{j}
    @#endfor
    @#for j in 1:24
        lm_psi@{j}
    @#endfor
    @#for j in 0:(24-1)
        lm_v@{j}
    @#endfor
    @#for j in 1:(24-1)
        lm_piacc@{j}
    @#endfor
    @#for j in 1:(24-1)
        lm_thetaacc@{j}
    @#endfor
;    

varexo
    ug          // shock to government consumption
    uz          // shock to the aggregate productivity process
    ucp         // cost-push innovation
;

parameters
    // We first declare the parameters that are common to all the models that we are going to consider
    gamma       // coefficient of RRA
    chi         // scale parameter of labor disutility
    varphi      // convexity of labor disutility
    epsilon     // elasticity of substitution across consumption varieties
    beta        // subjective discount factor
    rhoz rhog rhocp  // autocorrelation coefficients of 'z', 'g', 'cp shock'
    gbar        // steady state value of 'g'
    delta       // (one minus) retirement rate
    // We now declare the parameters that are specific to this particular model
    alpha miniepsilon
    piacc0 thetaacc0 lambda24
;


gamma       = 2;
chi         = 6;
varphi      = 1;
epsilon     = 7;
beta        = (1/1.04)^(1/4);
rhoz        = .95;
rhog        = .9;
rhocp       = .8;
delta       = 1;//.95;
miniepsilon = 0.0000000001;
alpha       = .0006; 
piacc0      = 1;
thetaacc0   = 1;
lambda24    = 1;
gbar = .1;


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

        cp_shock = cp_shock(-1)^rhocp*exp(ucp);

    // We now write the equations that are are specific to this particular model

        pstar = cp_shock*epsilon/(epsilon-1)*(c^(-gamma)*y*w/z
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

        @#for j in 1:(24-1)
            gL@{j} = 1/(alpha + (v0-v@{j})/w) - (miniepsilon + (v0-v@{j})/w)/(alpha + (v0-v@{j})/w)^2; // PDF of the menu cost distribution, evaluated at each cohort's adjustment gain
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

    // FOCs of the Ramsey problem

    0 = c^(-gamma) + lm_w*(-gamma)*c^(-gamma-1)*w + ((-gamma)*c^(-gamma-1)*y + c^(-gamma))*(
                                                    @#for j in 0:(24-1)
                                                       + lm_pstar(-@{j})*(pstar(-@{j})/piacc@{j} - epsilon/(epsilon-1)*w/z)*thetaacc@{j}*piacc@{j}^epsilon
                                                    @#endfor
                                                     ) - lm_n*Delta/z
                                                    @#for j in 0:(24-1)
                                                       + lm_v@{j}*((pstar(-@{j})/piacc@{j} - w/z)*(pstar(-@{j})/piacc@{j})^(-epsilon)*((-gamma)*c^(-gamma-1)*y + c^(-gamma)) - (-gamma)*c^(-gamma-1)*v@{j})
                                                    @#endfor
                                                    + (-gamma)*c^(-gamma-1)*(
                                                    @#for j in 0:(24-2)
                                                       lm_v@{j}(-1)*(lambda@{j+1}*v0 + (1-lambda@{j+1})*v@{j+1} - Xi@{j+1}) +
                                                    @#endfor
                                                       lm_v@{24-1}(-1)*(delta*v0 + (1-delta)*y/epsilon/(1-beta)) ) + lm_v@{24-1}(-1)*(1-delta)/epsilon/(1-beta)*c^(-gamma);

    0 = - chi*n^varphi - lm_w*chi*varphi*n^(varphi-1) + lm_n;

    0 = - lm_n*y/z - lm_Delta;

    0 = lm_pstar*(c^(-gamma)*y
                  @#for j in 1:(24-1)
                     + beta^@{j}*thetaacc@{j}(+@{j})*(piacc@{j}(+@{j}))^(epsilon-1)*c(+@{j})^(-gamma)*y(+@{j})
                  @#endfor
                  ) + (lm_pi*(
                       @#for j in 1:24
                          + lambda@{j}*psi@{j}
                       @#endfor
                       )
                        @#for j in 1:(24-1)
                           + beta^@{j}*lm_pi(+@{j})*(piacc@{j}(+@{j}))^(epsilon-1)*(1-lambda@{j}(+@{j}))*psi@{j}(+@{j})
                        @#endfor
                        )*(1-epsilon)*pstar^(-epsilon) + (lm_Delta*(
                                                                    @#for j in 1:24
                                                                       + lambda@{j}*psi@{j}
                                                                    @#endfor
                                                                    )
                                                                     @#for j in 1:(24-1)
                                                                        + beta^@{j}*lm_Delta(+@{j})*(piacc@{j}(+@{j}))^epsilon*(1-lambda@{j}(+@{j}))*psi@{j}(+@{j})
                                                                     @#endfor
                                                                     )*(-epsilon)*pstar^(-epsilon-1) 
        + lm_v0*(epsilon*w/z - (epsilon-1)*pstar)*pstar^(-epsilon-1)*y*c^(-gamma)
        @#for j in 1:(24-1) 
            + beta^@{j}*lm_v@{j}(+@{j})*(epsilon*w(+@{j})/z(+@{j}) - (epsilon-1)*pstar/piacc@{j}(+@{j}))*pstar^(-epsilon-1)*(piacc@{j}(+@{j}))^epsilon*y(+@{j})*c(+@{j})^(-gamma)
        @#endfor
        ;

    0 = lm_w*c^(-gamma) - epsilon/(epsilon-1)*c^(-gamma)*y/z*(
                                                              @#for j in 0:(24-1)
                                                                 + lm_pstar(-@{j})*thetaacc@{j}*piacc@{j}^epsilon
                                                              @#endfor
                                                              ) + lm_n*(
                                                                        @#for j in 1:(24-1)
                                                                            + psi@{j}*((v0-v@{j})/w)^2/w*gL@{j}
                                                                        @#endfor
                                                                        )
                                                                - c^(-gamma)*y/z*(
                                                                                  @#for j in 0:(24-1)
                                                                                     + lm_v@{j}*(pstar(-@{j})/piacc@{j})^(-epsilon)
                                                                                  @#endfor
                                                                                  ) 
        @#for j in 1:(24-1)
            + lm_lambda@{j}*gL@{j}*(v0-v@{j})/w^2    
        @#endfor
        @#for j in 0:(24-2)
            - lm_v@{j}(-1)*c^(-gamma)*(Xi@{j+1}/w - ((v0-v@{j+1})/w)^2*gL@{j+1})
        @#endfor
        ;

    0 = - lm_piacc1
                    @#for j in 2:(24-1)
                        - lm_piacc@{j}*piacc@{j-1}(-1)
                    @#endfor
                     ;

    @#for j in 1:(24-2)
        0 = lm_pstar(-@{j})*thetaacc@{j}*piacc@{j}^(epsilon-1)*((epsilon-1)*pstar(-@{j})/piacc@{j} - epsilon^2/(epsilon-1)*w/z)*y*c^(-gamma)
            + (lm_pi*(epsilon-1)*pstar(-@{j})/piacc@{j} + lm_Delta*epsilon)*pstar(-@{j})^(-epsilon)*piacc@{j}^(epsilon-1)*(1-lambda@{j})*psi@{j}
            + lm_v@{j}*((epsilon-1)*pstar(-@{j})/piacc@{j} - epsilon*w/z)*pstar(-@{j})^(-epsilon)*piacc@{j}^(epsilon-1)*y*c^(-gamma) + lm_piacc@{j} - beta*lm_piacc@{j+1}(+1)*pi(+1);
    @#endfor
        0 = lm_pstar(-@{24-1})*thetaacc@{24-1}*piacc@{24-1}^(epsilon-1)*((epsilon-1)*pstar(-@{24-1})/piacc@{24-1} - epsilon^2/(epsilon-1)*w/z)*y*c^(-gamma)
            + (lm_pi*(epsilon-1)*pstar(-@{24-1})/piacc@{24-1} + lm_Delta*epsilon)*pstar(-@{24-1})^(-epsilon)*piacc@{24-1}^(epsilon-1)*(1-lambda@{24-1})*psi@{24-1}
            + lm_v@{24-1}*((epsilon-1)*pstar(-@{24-1})/piacc@{24-1} - epsilon*w/z)*pstar(-@{24-1})^(-epsilon)*piacc@{24-1}^(epsilon-1)*y*c^(-gamma) + lm_piacc@{24-1};

    @#for j in 1:(24-2)
        0 = lm_pstar(-@{j})*piacc@{j}^epsilon*(pstar(-@{j})/piacc@{j} - epsilon/(epsilon-1)*w/z)*y*c^(-gamma) + lm_thetaacc@{j} - beta*lm_thetaacc@{j+1}(+1)*(1-lambda@{j+1}(+1));
    @#endfor
    0 = lm_pstar(-@{24-1})*piacc@{24-1}^epsilon*(pstar(-@{24-1})/piacc@{24-1} - epsilon/(epsilon-1)*w/z)*y*c^(-gamma) + lm_thetaacc@{24-1};

    @#for j in 1:(24-1)
        0 = lm_pi*(pstar^(1-epsilon) - (pstar(-@{j})/piacc@{j})^(1-epsilon))*psi@{j} + lm_Delta*(pstar^(-epsilon) - (pstar(-@{j})/piacc@{j})^(-epsilon))*psi@{j}
            + lm_lambda@{j} + beta*lm_psi@{j+1}(+1)*psi@{j} + lm_v@{j-1}(-1)*c^(-gamma)*(v0-v@{j}) + lm_thetaacc@{j}*thetaacc@{j-1}(-1);
    @#endfor

    @#for j in 2:(24-1)
        0 = - lm_n*Xi@{j}/w + lm_pi*(pstar^(1-epsilon)*lambda@{j} + (pstar(-@{j})/piacc@{j})^(1-epsilon)*(1-lambda@{j})) +
            lm_Delta*(pstar^(-epsilon)*lambda@{j} + (pstar(-@{j})/piacc@{j})^(-epsilon)*(1-lambda@{j})) + lm_psi@{j} - beta*lm_psi@{j+1}(+1)*(1-lambda@{j}) + lm_psi1;
    @#endfor
    0 = - lm_n*Xi1/w + lm_pi*(pstar^(1-epsilon)*lambda1 + (pstar(-1)/piacc1)^(1-epsilon)*(1-lambda1)) +
        lm_Delta*(pstar^(-epsilon)*lambda1 + (pstar(-1)/piacc1)^(-epsilon)*(1-lambda1)) + lm_psi1 - beta*lm_psi2(+1)*(1-lambda1);
    0 = lm_pi*pstar^(1-epsilon) + lm_Delta*pstar^(-epsilon) + lm_psi@{24} + lm_psi1;

    @#for j in 1:(24-1)
        0 = lm_n*psi@{j}*(v0-v@{j})/w^2*gL@{j} + lm_lambda@{j}*gL@{j}/w - lm_v@{j}*c^(-gamma) + lm_v@{j-1}(-1)*c^(-gamma)*(1 - lambda@{j} + (v0-v@{j})/w*gL@{j});
    @#endfor
    0 = - lm_n*(
                @#for j in 1:(24-1)
                    + psi@{j}*(v0-v@{j})/w^2*gL@{j}
                @#endfor
                )
        @#for j in 1:(24-1)
            - lm_lambda@{j}*gL@{j}/w      
        @#endfor
        - lm_v0*c^(-gamma) + lm_v@{24-1}(-1)*c^(-gamma)*delta
                            @#for j in 0:(24-2)
                                + lm_v@{j}(-1)*c^(-gamma)*(lambda@{j+1} - (v0-v@{j+1})/w*gL@{j+1})
                            @#endfor
                            ;
end;


shocks;
    var ug; stderr .01; 
    var uz; stderr .01;
    var ucp; stderr .01;
end;


lambdabar = miniepsilon/alpha; % This is the value of 'lambda' at L=0

A = zeros(24,24);
A(1,2:24) = -1;
for j=2:24
    A(j,j-1) = 1-lambdabar;
end;
B = [1; zeros(24-1,1)];
psiss = inv(eye(24)-A)*B;

yss = fsolve(@(y) epsilon/(epsilon-1)*chi*y^varphi*(y-gbar)^gamma - 1, .24);

initval;
    y       = yss;  //((epsilon-1)/epsilon/chi)^(1/(varphi+gamma));     //nss/Deltass;
    c       = y - gbar;    //nss/Deltass;
    n       = y;    //nss;
    r       = 1/beta;   //pibar/beta;
    pi      = 1;    //pibar;
    pstar   = 1;    //pstarss;
    w       = chi*n^varphi*c^gamma; //wss;
    Delta   = 1;    //Deltass;
    z       = 1;
    cp_shock = 1;
    g       = gbar;

    @#for j in 1:(24-1)
        lambda@{j} = lambdabar;         //lambdass(@{j}); 
    @#endfor
    
    @#for j in 1:(24-1)
        piacc@{j} = 1;              //pibar^@{j}; 
    @#endfor

    thetaacc1 = 1-lambda1;
    @#for j in 2:(24-1)
        thetaacc@{j} = (1-lambda@{j})*thetaacc@{j-1};
    @#endfor

    @#for j in 1:24
        psi@{j} = psiss(@{j}); 
    @#endfor

    @#for j in 0:(24-1)
        v@{j} = y/epsilon/(1-beta);     //vss(@{j+1}); 
    @#endfor

    @#for j in 1:(24-1)
        Xi@{j} = 0;
    @#endfor

    @#for j in 1:(24-1)
        gL@{j} = 1/alpha - miniepsilon/alpha^2;
    @#endfor

    lm_w = (c^(-gamma) - chi*n^varphi)/(chi*varphi*n^(varphi-1) - (-gamma)*c^(-gamma-1)*w);
    lm_n = chi*n^varphi + lm_w*chi*varphi*n^(varphi-1);
    lm_Delta = - lm_n*y;
    lm_pstar = (epsilon-1)/epsilon*lm_w/y/(
                                           @#for j in 0:(24-1)
                                                + (1-lambdabar)^@{j}
                                           @#endfor
                                           );
    lm_pi = 1/(epsilon-1)*(lm_pstar*c^(-gamma)*y*(
                                                   @#for j in 0:(24-1)
                                                        + (beta*(1-lambdabar))^@{j}
                                                   @#endfor
                                                  )/(
                                                   @#for j in 1:(24-1)
                                                        (lambdabar + (1-lambdabar)*beta^@{j})*psi@{j} +
                                                   @#endfor
                                                     psi@{24}) - epsilon*lm_Delta);

    @#for j in 1:(24-1)
        lm_lambda@{j} = 0;
    @#endfor
    lm_psi1 = -(lm_pi + lm_Delta);
    @#for j in 2:24
        lm_psi@{j} = 0;
    @#endfor
    @#for j in 0:(24-1)
        lm_v@{j} = 0;
    @#endfor
    @#for j in 1:(24-1)
        lm_piacc@{j} = 0;
    @#endfor
    @#for j in 1:(24-1)
        lm_thetaacc@{j} = 0;
    @#endfor
end;


steady(solve_algo=0);
stoch_simul(nograph,order=1);
save smc