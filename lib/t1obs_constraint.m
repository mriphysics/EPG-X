function [c,ceq] = t1obs_constraint(x,T1obs,eps)
% Function to enforce inequality constraint during model fitting.
% x is parameter vector [k*1000 f T1a T1b Gsf] Gsf is scaling factor applied to G0, 
% T1a,b in seconds
%
% T1obs is the observed T1 that the fit parameters must be consistent with
% eps is the allowed deviation from T1obs
% 
% use linear inequality constraint (T1obs(pars)-T1obs)^2<eps^2

R1f = 1e-3/x(3);
kf = x(1)*1e-3;
R1b = 1e-3/x(4);
kb = kf * (1-x(2))/x(2);

if kf==0
    kb=0;
end
%%% the idea is that the observerd T1 should match the combined version
%%% from R1f, R1b, kf and f
R1obs = 0.5*(R1f + kf + R1b + kb)-0.5*sqrt((R1f + kf + R1b + kb).^2 ...
    -4*(R1f*R1b + R1f*kb + R1b*kf));

ceq = [];
c = ((1/R1obs - T1obs).^2) - eps^2; %<- eps is the tolerance on the solution

end