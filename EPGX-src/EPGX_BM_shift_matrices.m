function S = EPGX_BM_shift_matrices(Nmax)
% Generates shift matrices up to order Nmax for 6 component 2 pool Bloch
% McConnell case
% Layout of state vector is [F0A F0*A Z0A F0B F0*B Z0B F1A F-1A Z1A F1B F-1B Z1B ... etc]
% i.e. 6 states per n-value
%
%   Shaihan Malik 2017-07-19

N = (Nmax+1) * 6;

S = zeros([N N]);

%%% In each case below, sidx is the ...

%%% FA(k>=1) look @ states just BEFORE pulse
kidx = 7:6:N; % F+A states starting at F1A
sidx = kidx-6;% <- Positive states are shifted from 6 back
idx = kidx + N*(sidx-1); % linear indices of S(kidx,sidx)
S(idx)=1;

%%% FB(k>=1) look @ states just BEFORE pulse
kidx = 10:6:N; % F+B states starting at F1B
sidx = kidx-6;% <- Positive states are shifted from 6 back
idx = kidx + N*(sidx-1);
S(idx)=1;

%%% FA(k<1) <--- start at F-1A
kidx = 8:6:N; % F-A states starting from F-1A
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+6;%<- Negative states come from more negaative state with higher index
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% FB(k<1) <--- start at F-BA
kidx = 11:6:N; % F-B states starting from F-1B
kidx(end)=[];% most negative state relates to nothing; related row is empty
sidx = kidx+6;%<- Negative states come from more negative state with higher index
ix = kidx + N*(sidx-1);
S(ix)=1;

%%% Z states, these don't shift
% ZA states
kidx = 3:6:N ;
ix = kidx + N*(kidx-1);
S(ix)=1;
% ZB states
kidx = 6:6:N ;
ix = kidx + N*(kidx-1);
S(ix)=1;


%%% finally F0A - relates to F-1A
S([1 2],8)=1;
% Same for F0B
S([4 5],11)=1;

end