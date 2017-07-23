%% Gloor SSFP function 23-5-2017
% from Gloor, M., Scheffler, K. & Bieri, O. Quantitative magnetization transfer 
% imaging using balanced SSFP. Magn. Reson. Med. 60, 691–700 (2008).
function Mss = ssSSFP_Gloor(alpha,b1sqrdtau, TR, T1x,T2,f,ka,G)

M0a = (1-f);
M0b = f;
kb = ka * M0a/M0b;
R1a = 1/T1x(1);
R1b = 1/T1x(2);
R2 = 1/T2;

G = G *1e-3; %<- convert from us to ms
gam = 267.5221 *1e-3; %< rad /ms /uT
WT = pi*gam^2*b1sqrdtau*G; 

F = M0b/M0a;
fw = exp(-WT);
fk = exp(-(ka+kb)*TR);

E1b = exp(-R1b*TR);
E2 = exp(-R2*TR);
E1a = exp(-R1a*TR);

A = 1+F-fw*E1b*(F+fk);
B = 1+fk*(F-fw*E1b*(F+1));
C = F*(1-E1b)*(1-fk);

Mss = M0a*sin(alpha)*((1-E1a)*B+C)/(A-B*E1a*E2-(B*E1a-A*E2)*cos(alpha));

end