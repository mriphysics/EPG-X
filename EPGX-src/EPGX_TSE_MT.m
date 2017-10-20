function [F0,Fn,Zn,F] = EPGX_TSE_MT(theta,B1SqrdTau,ESP,T1x,T2a,f,ka,G,varargin)
%%%% [F0,Fn,Zn,F] = EPGX_TSE_MT(theta,B1SqrdTau,ESP,T1x,T2a,f,ka,G,varargin)
%
%   EPG-X for Pulsed MT systems w/ gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               B1SqrdTau:  vector of RF pulse energies, uT^2 ms
%               ESP:        echo spacing, ms
%               T1x:        [T1a T1b], ms
%               T2a:        T2a, ms (only compartment a has appreciable T2)
%               f:          fraction of compartment b
%               ka:         exchange rate from a->b (units ms^-1)
%               G:          absorption line value at the frequency of
%                           interest. Units us (microseconds)
%
%   optional arguments (use string then value as next argument)
%
%               kmax:       maximum EPG order to include. Can be used to
%                           accelerate calculation. 
%                           Setting kmax=inf ensures ALL pathways are
%                           computed
%              diff:        structure with fields:
%                           G    - Gradient amplitude(s)
%                           tau  - Gradient durations(s)
%                           D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%               * experimental *
%
%               zinit:      User specified initial state of Z0 ([Z0a Z0b])
%
%   Outputs:                
%               F0:         signal (F0 state) at each echo time
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0a F0a* Z0a Z0b F1a F-1a* Z1a Z1b F2a F-2a* Z2a Z2b ...] etc
%
%   Shaihan Malik 2017-07-21


%% Extra variables

for ii=1:length(varargin)
    
    % Kmax = this is the maximum EPG 'order' to consider
    % If this is infinity then don't do any pruning
    if strcmpi(varargin{ii},'kmax')
        kmax = varargin{ii+1};
    end
    
    % Diffusion - structure contains, G, tau, D
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
    end
    
    % Zinit - allow the initial state to vary, but only allow the Z0 term
    % to be different, assume all other terms are zero
    if strcmpi(varargin{ii},'zinit')
        zinit=varargin{ii+1};% 2 element vector [Z0a Z0b]
    end
    
end


%% Calculation is faster when considering only appropriate EPG orders 

%%% The maximum order varies through the sequence. This can be used to speed up the calculation    
np = length(theta); % number of pulses, this includes exciatation

% if not defined, assume want max
if ~exist('kmax','var')
    kmax = 2*(np - 1);  % kmax at last echo time, this is twice # refocusing pulses
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = 2*(np - 1); % this is maximum value
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = 2*(2:np);
else
    kmax_per_pulse = 2*[1:ceil(np/2) (floor(np/2)):-1:1]+1;
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
     
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end

%%% Number of states is 4x(kmax +1) -- +1 for the zero order
N=4*(kmax+1);

%%
%%% Sort out RF pulse amplitude and phases
% split magnitude and phase
alpha = abs(theta);phi=angle(theta);

% add CPMG phase
phi(2:end) = phi(2:end) + pi/2;

    
%% Dependent variables for exchange case
M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;

R1a = 1/T1x(1);
R1b = 1/T1x(2);
R2a = 1/T2a;

%%% handle no exchange case
if (f==0)||(ka==0)
    ka = 0;
    kb = 0;
    M0b = 0;
    R1b = 1e-3;%<- doesn't matter what it is, just avoid singular matrices
end

%%% Build Shift matrix, S
S = EPGX_MT_shift_matrices(kmax);
S = sparse(S);



%% Set up matrices for Relaxation and Exchange

%%% Relaxation-exchange matrix for transverse components
Lambda_T = diag([-R2a -R2a]);
Xi_T = expm(0.5*ESP*Lambda_T); %<-- operator for time period ESP/2

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1a-ka kb];[ka -R1b-kb]];
Xi_L = expm(0.5*ESP*Lambda_L); 

% For inhomogeneous solution also need:
C_L = [M0a*R1a;M0b*R1b];
Zoff = (Xi_L - eye(2))*(Lambda_L\C_L);

% Now build full matrix Xi
Xi = blkdiag(Xi_T,Xi_L);

%%% regrowth
b = zeros([N 1]);
b([3 4]) = Zoff;%<--- just applies to the Z0b and Z0b states, not Z1,2 etc

%%% Add in diffusion at this point - this is an experimental feature
if exist('diff','var')
   Xi = Xi_diff_MT(Xi_T,Xi_L,diff,kmax,N);
else
   % If no diffusion, Xi is the same for all EPG orders
   Xi = kron(eye(kmax+1),Xi);
end
    

%%% Composite exchange-relax-shift
XS=S*Xi;
XS=sparse(XS);

%%% Compute saturation terms for RF pulse
G = G *1e-3; %<- convert from us to ms
gam = 267.5221 *1e-3; %< rad /ms /uT
WT = pi*gam^2*B1SqrdTau*G; 


%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);
% store the indices of the top 4x4 corner, this helps build_T
i1 = [];
for ii=1:4
    i1 = cat(2,i1,sub2ind(size(T),1:4,ii*ones(1,4)));
end


%% F matrix - records the state at each echo time
F = zeros([N np-1]); % there are np-1 echoes because of excitation pulse

%%% Initial State - can be specified by user
FF = zeros([N 1]);
if ~exist('zinit','var')
    FF(3)=1-f;
    FF(4)=f;
else
    FF(3)=zinit(1);
    FF(4)=zinit(2);
end
    

%% Now handle excitation 

A = RF_rot_sat(alpha(1),phi(1),WT(1));
kidx = 1:4;
FF(kidx) = A*FF(kidx); %<---- state straight after excitation, zero order only


%%% Now simulate the dephase gradient & evolution, half the readout
kidx=1:8;
FF(kidx) = XS(kidx,kidx)*FF(kidx)+b(kidx);


%% Now simulate the refocusing pulses

for jj=2:np 
    A = RF_rot_sat(alpha(jj),phi(jj),WT(jj));
    build_T(A);%<- replicate A to make large T matrix
    
    % variable maximum EPG order - accelerate calculation
    kidx = 1:4*kmax_per_pulse(jj);
    
    % Apply RF pulse to current state
    FF(kidx)=T(kidx,kidx)*FF(kidx);

    % Now evolve for half echo spacing, store this as the echo
    F(kidx,jj-1) = XS(kidx,kidx)*FF(kidx)+b(kidx);
    % Deal with complex conjugate after shift
    F(1,jj-1)=conj(F(1,jj-1)); %<---- F0 comes from F-1 so conjugate
    
    if jj==np
        break
    end
    
    % Finally, evolve again up to next RF pulse
    FF(kidx) = XS(kidx,kidx)*F(kidx,jj-1)+b(kidx);
    FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate 
    
end

%%% Return signal
F0 = F(1,:)*1i;

%%% Construct Fn and Zn
idx=[fliplr(6:4:size(F,1)) 1 5:4:size(F,1)]; 
kvals = -kmax:kmax;
%%% Remove the lowest two negative states since these are never populated
%%% at echo time
idx(1:2)=[];
kvals(1:2)=[];

%%% Now reorder
Fn = F(idx,:);
%%% Conjugate
Fn(kvals<0,:)=conj(Fn(kvals<0,:));

%%% Similar for Zn
ZnA = F(3:4:end,:);
ZnB = F(4:4:end,:);
Zn = cat(3,ZnA,ZnB);



    %%% Modified EPG transition matrix including saturation
    % As per Weigel et al JMR 2010 276-285 plus additional 4th row/col
    function Tap = RF_rot_sat(a,p,WT)
        Tap = zeros([4 4]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(5) = conj(Tap(2));
        Tap(6) = Tap(1);
        Tap(7) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(9) = -1i*exp(1i*p)*sin(a);
        Tap(10) = 1i*exp(-1i*p)*sin(a);
        Tap(11) = cos(a);
        Tap(16) = exp(-WT); %<--- this value is bound pool saturation
    end

    function build_T(AA)
        ksft = 4*(4*(kmax+1)+1);
        for i2=1:16
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
