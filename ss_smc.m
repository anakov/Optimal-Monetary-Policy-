% This code computes the steady state with trend inflation (pibar) in the
% stochastic menu cost model


close all
clear
clc

gamma       = 2;
chi         = 6;
varphi      = 1;
epsilon     = 7;
beta        = (1/1.04)^(1/4);
pibar       = (1.02)^(1/4);
miniepsilon = 0.0000000001;
alpha       = .0006; 
gbar = .1;


% This algorithm starts with an initial guess for the vectors of lambda's and Xi's

J = 24; % Number of cohorts

lambda = ones(1,J-1);
theta =  1-lambda;
Xi = zeros(J,1);

vecdiff = 1;
while vecdiff > .000001;
    
    A = zeros(J,J);
    A(1,2:J) = -1;
    for j=2:J
        A(j,j-1) = 1-lambda(j-1);
    end;
    B = [1; zeros(J-1,1)];
    psi = inv(eye(J)-A)*B;

    pstar = ([lambda 1]*psi + sum((pibar^(epsilon-1)).^(1:1:J-1).*(1-lambda).*psi(1:J-1)'))^(1/(epsilon-1));

    Delta = ([lambda 1]*psi + sum((pibar^epsilon).^(1:1:J-1).*(1-lambda).*psi(1:J-1)'))*pstar^(-epsilon);

    cumprodtheta = cumprod(theta);
    w = (epsilon-1)/epsilon*pstar*(1 + sum((beta*pibar^(epsilon-1)).^(1:1:J-1).*cumprodtheta)) ...
                                 /(1 + sum((beta*pibar^epsilon).^(1:1:J-1).*cumprodtheta));

    y = fsolve(@(y) chi*(y*Delta + psi(1:J-1)'*Xi(1:J-1)/w)^varphi*(y-gbar)^gamma - w, .5);    %(w*Delta^gamma/chi)^(1/(varphi+gamma));
    n = y*Delta + psi(1:J-1)'*Xi(1:J-1)/w;
    c = y - gbar;

    F = (pstar./pibar.^(0:1:J-1) - w).*(pstar./pibar.^(0:1:J-1)).^(-epsilon)*y;
    A = zeros(J,J);
    A(:,1) = [lambda'; 1];
    for j=1:J-1
        A(j,j+1) = 1-lambda(j);
    end;
    v = (eye(J)-beta*A)\(F'-beta*Xi);
    v0 = v(1);
    v = v(2:J); 
    
    for j=1:J-1
        Xi(j) = w*quadl(@(k) k.*(1./(alpha+k)-(miniepsilon+k)./(alpha+k).^2), 0, (v0-v(j))/w);    %lambdabar/(lambdabar+(1-lambdabar)*(alpha/(max(0.0000001,v0-v(j))/w))^xi);
    end;

    for j=1:J-1
        lambdaa(j) = (miniepsilon + (v0-v(j))/w)/(alpha + (v0-v(j))/w);    %lambdabar/(lambdabar+(1-lambdabar)*(alpha/(max(0.0000001,v0-v(j))/w))^xi);
    end;
    vecdiff = max(abs(lambdaa-lambda));
    jump = .5;
    lambda = jump*lambdaa + (1-jump)*lambda;
    
    for j=1:J-1
        theta(j) = 1 - lambda(j);% - lambda(j)*(1-lambda(j));     %lambda(j)^2*xi*(1-lambdabar)/lambdabar*(alpha/(max(0.0000001,v0-v(j))/w))^xi;
    end;

    
end;

Util = c^(1-gamma)/(1-gamma) - chi*n^(1+varphi)/(1+varphi);
Elambda = [lambda 1]*psi;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%save endogvarss c n y pstar w Delta psi lambda theta v0 v pibar