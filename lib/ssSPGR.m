%%% 22-7-17: Ernst function for ssSPGR
function sig = ssSPGR(theta,TR,T1)
E1=exp(-TR./T1);
sig = sin(theta).*(1-E1)./(1-E1.*cos(theta));
end