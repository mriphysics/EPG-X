%% 22-7-2017 Steady state SGPR for MT
function Mss = ssSPGR_MT(theta,b1sqrdtau, TR, T1x,f,ka,G)

M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;
R1a = 1/T1x(1);
R1b = 1/T1x(2);

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1a-ka kb];[ka -R1b-kb]];
Xi_L = expm(TR*Lambda_L); 

% For inhomogeneous solution also need:
C_L = [M0a*R1a;M0b*R1b];

%%% Compute saturation terms for RF pulse
G = G *1e-3; %<- convert from us to ms
gam = 267.5221 *1e-3; %< rad /ms /uT
WT = pi*gam^2*b1sqrdtau*G; 


T = diag([cos(theta) exp(-WT)]);
I = eye(2);
S = [sin(theta) 0];
Mss = S * inv(I-Xi_L*T) * (Xi_L - I) * inv(Lambda_L) * C_L;


end