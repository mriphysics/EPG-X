%% 25-6-17: steady state SINGLE pool
function [Mss] = ssSSFP(alpha, TR, dw,  T1, T2)


R1 = 1/T1;
R2 = 1/T2;

%%% Full system matrix - same notation as CX version
E = diag([-R2 -R2 -R1]);
W  = diag([-1i*dw 1i*dw 0]);
A = E + W;
eA = expm(A*TR);
C = [0 0 R1].';

% Assume RF is rotation around x-axis. Phase doesn't matter here
T = RF_rot(alpha,0);

I = eye(3);

% Z-rotation 180
D = diag([-1 -1 1]);

% SS solution
Mss = T * inv(D - eA*T) * (eA - I) * inv(A) * C;

%%% Just return the Mxy component 
Mss = (Mss(1));

    %%% NORMAL EPG transition matrix 
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
    end

end