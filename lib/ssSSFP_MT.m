%% 13-6-17: steady state with MT, balanced -- now matrix version, not analytic
function [Mss] = ssSSFP_MT(alpha,b1sqrdtau, TR,dw, T1x,T2,f,ka,G)

M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;
R1a = 1/T1x(1);
R1b = 1/T1x(2);
R2a = 1/T2;

G = G *1e-3; %<- convert from us to ms
gam = 267.5221 *1e-3; %< rad /ms /uT
WT = pi*gam^2*b1sqrdtau*G; 

%%% Full system matrix including off resonance in "+" basis: BUT MT case
Ea = diag([-R2a -R2a -R1a]);
Eb = diag([0 0 -R1b]);
Ka = diag([0 0 ka]);%<--- No transverse coupling
Kb = diag([0 0 kb]);
W  = diag([-1i*dw 1i*dw 0]);

A = [[Ea-Ka Kb];[Ka Eb-Kb]] + kron(eye(2),W);
% delete rows and columns for Mxy in pool b
A([4 5],:)=[];
A(:,[4 5])=[];

eA = expm(A*TR);
C = [0 0 R1a*M0a R1b*M0b].';

% Assume RF is rotation around x-axis. Phase doesn't matter here
T = RF_rot(alpha,0,WT);

I = eye(4);

D = diag([-1 -1 1 1]);

% Component to return
S = [1 0 0 0];

% SS solution
Mss = S*T * inv(D - eA*T) * (eA - I) * inv(A) * C;



    %%% NORMAL EPG transition matrix but add MT 
    % As per Weigel et al JMR 2010 276-285 
    function T = RF_rot(a,p,WT)
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
        T = blkdiag(T,exp(-WT));
    end

end