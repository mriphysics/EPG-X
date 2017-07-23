%% 12-6-17: SSFP steady state, BM equations
function [Mss] = ssSSFP_BM(alpha, TR, dw,  T1x, T2x, f,ka)

M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;
R1a = 1/T1x(1);
R1b = 1/T1x(2);
R2a = 1/T2x(1);
R2b = 1/T2x(2);

%%% Full system matrix including off resonance in "+" basis
Ea = diag([-R2a -R2a -R1a]);
Eb = diag([-R2b -R2b -R1b]);
Ka = ka*eye(3);
Kb = kb*eye(3);
W  = diag([-1i*dw 1i*dw 0]);

A = [[Ea-Ka Kb];[Ka Eb-Kb]] + kron(eye(2),W);
eA = expm(A*TR);
C = [0 0 R1a*M0a 0 0 R1b*M0b].';

% Assume RF is rotation around x-axis. Phase doesn't matter here
T = RF_rot(alpha,0);

I = eye(6);

% Z-rotation 180
D = diag([-1 -1 1 -1 -1 1]);

% Components to sum before returning
S = [1 0 0 1 0 0];

% SS solution
Mss = S* T * inv(D - eA*T) * (eA - I) * inv(A) * C;


    %%% NORMAL EPG transition matrix but replicated twice 
    % As per Weigel et al JMR 2010 276-285 
    function T = RF_rot(a,p)
        T = zeros([3 3]);
        T(1) = cos(a/2).^2;
        T(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*1i*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));
        T(5) = T(1);
        T(6) = 0.5*1i*exp(1i*p)*sin(a);
        T(7) = -1i*exp(1i*p)*sin(a);
        T(8) = 1i*exp(-1i*p)*sin(a);
        T(9) = cos(a);
        T = kron(eye(2),T);
    end

end