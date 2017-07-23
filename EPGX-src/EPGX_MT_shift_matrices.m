function S = EPGX_MT_shift_matrices(Nmax)
% Generates shift matrices up to order Nmax for 4 component 2 pool MT case
% Layout of state vector is [F0A F0*A Z0A Z0B F1A F-1A Z1A Z1B ... etc]
% i.e. 4 states per n-value
%
%   Shaihan Malik 2017-07-20

N = (Nmax+1) * 4;
S = zeros([N N]);

%%% F(k>=1) look @ states just BEFORE np+1 pulse
kidx = 5:4:N; % miss out F1+
sidx = kidx-4;%<- Negative states shift 4 back
idx = kidx + N*(sidx-1); % linear indices of S(kidx,sidx)
S(idx)=1;

%%% F(k<1) <--- start at F-1 
kidx = 6:4:N;
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+4;%<- Positive states shift 4 forward
ix = kidx + N*(sidx-1); % linear indices of S(kidx,sidx)
S(ix)=1;

%%% Z states, these don't shift
% Za states
kidx = 3:4:N ;
ix = kidx + N*(kidx-1);
S(ix)=1;
% Zb states
kidx = 4:4:N ;
ix = kidx + N*(kidx-1);
S(ix)=1;


%%% finally F0+ - relates to F-1
S([1 2],6)=1; % also F0* - this isn't used in practice


end