%% 22-7-17: steady state SPGR with BM equations
function [Mss] = ssSPGR_BM(alpha, TR, T1x, f,ka)

M0a = (1-f);
M0b = f;
kb = ka * M0a/M0b;
R1a = 1/T1x(1);
R1b = 1/T1x(2);

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1a-ka kb];[ka -R1b-kb]];
Xi_L = expm(TR*Lambda_L);


% Use same notation as paper
I = eye(2);
T = cos(alpha)*I;
C_L = [R1a*M0a;R1b*M0b];

Mss = [1 1]*sin(alpha)* inv(I-Xi_L*T) * (Xi_L - I) * inv(Lambda_L) * C_L;


end