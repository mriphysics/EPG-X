function [F0,Fn,Zn,F] = EPGX_TSE_BM(theta,ESP,T1x,T2x,f,ka,varargin)
%%%% [F0,Fn,Zn,F] = EPGX_TSE_BM(theta,ESP,T1x,T2x,f,ka,varargin)
%
%   EPG-X for Bloch McConnell coupled systems w/ TSE sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               ESP:        echo spacing, ms
%               T1x:        [T1a T1b], ms
%               T2x:        [T2a T2b], ms
%               f:          fraction of compartment b
%               ka:         forward exchange rate from a->b (units ms^-1, i.e. kHz)
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
%           * Diffusion is same in both compartments, this is experimental *
%
%              delta:       Frequency offset for compartment b, kHz
%
%   Outputs:                
%               F0:         signal (F0 state) at each echo time
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0a F0a* Z0a F0b F0b* Z0b F1a F-1a* Z1a F1b F-1b* Z1b ...] etc
%
%
%   Shaihan Malik 2017-10-13


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
    
    % chemical shift of pool b, for phase gain during evolution period
    if strcmpi(varargin{ii},'delta')
        delta = 2*pi*varargin{ii+1};
    end
end

if ~exist('delta','var')
    delta=0;
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

%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=6*(kmax+1);

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
R2a = 1/T2x(1);
R2b = 1/T2x(2);

%%% handle no exchange case
if (f==0)||(ka==0)
    ka = 0;
    kb = 0;
    M0b = 0;
    R1b = 1e-3;%<- doesn't matter what it is, just avoid singular matrices
end

%%% Build Shift matrix, S
S = EPGX_BM_shift_matrices(kmax);
S = sparse(S);



%% Set up matrices for Relaxation and Exchange

%%% Relaxation-exchange matrix for transverse components
Lambda_T = diag([-R2a-ka -R2a-ka -R2b-kb-1i*delta -R2b-kb+1i*delta]);
Lambda_T(1,3) = kb;
Lambda_T(2,4) = kb;
Lambda_T(3,1) = ka;
Lambda_T(4,2) = ka;
Xi_T = expm(0.5*ESP*Lambda_T); %<-- operator for time period ESP/2

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1a-ka kb];[ka -R1b-kb]];
Xi_L = expm(0.5*ESP*Lambda_L); 

% For inhomogeneous solution also need:
C_L = [M0a*R1a;M0b*R1b];
Zoff = (Xi_L - eye(2))*(Lambda_L\C_L);

% Now build full matrix Xi
Xi = zeros(6);
Xi([1 2 4 5],[1 2 4 5])=Xi_T;
Xi([3 6],[3 6])=Xi_L;


%%% regrowth
b = zeros([N 1]);
b([3 6]) = Zoff;%<--- just applies to the Z0b and Z0b states, not Z1,2 etc

%%% Add in diffusion at this point - this is an experimental feature
if exist('diff','var')
   Xi = Xi_diff_BM(Xi_T,Xi_L,diff,kmax,N);
else
   % If no diffusion, Xi is the same for all EPG orders
   Xi = kron(eye(kmax+1),Xi);
end
    

%%% Composite exchange-relax-shift
XS=S*Xi;
XS=sparse(XS);

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 6x6 corner, this helps build_T
i1 = [];
for ii=1:6
    i1 = cat(2,i1,sub2ind(size(T),1:6,ii*ones(1,6)));
end


%% F matrix - records the state at each echo time
F = zeros([N np-1]); % there are np-1 echoes because of excitation pulse

%%% Initial State
FF = zeros([N 1]);
FF(3)=1-f;   %Z0a
FF(6)=f;     %Z0b

%% Now handle excitation 

A = RF_rot(alpha(1),phi(1));
kidx = 1:6;
FF(kidx) = A*FF(kidx); %<---- state straight after excitatio, zero order only


%%% Now simulate the dephase gradient & evolution, half the readout
kidx=1:12;
FF(kidx) = XS(kidx,kidx)*FF(kidx)+b(kidx);


%% Now simulate the refocusing pulses

for jj=2:np 
    A = RF_rot(alpha(jj),phi(jj));
    build_T(A);%<- replicate A to make large T matrix
    
    % variable maximum EPG order - accelerate calculation
    kidx = 1:6*kmax_per_pulse(jj);
    
    % Apply RF pulse to current state
    FF(kidx)=T(kidx,kidx)*FF(kidx);

    % Now evolve for half echo spacing, store this as the echo
    F(kidx,jj-1) = XS(kidx,kidx)*FF(kidx)+b(kidx);
    % Deal with complex conjugate after shift
    F([1 4],jj-1)=conj(F([1 4],jj-1)); %<---- F0 comes from F-1 so conjugate (do F0a and F0b)
    
    if jj==np
        break
    end
    
    % Finally, evolve again up to next RF pulse
    FF(kidx) = XS(kidx,kidx)*F(kidx,jj-1)+b(kidx);
    FF([1 4])=conj(FF([1 4])); %<---- F0 comes from F-1 so conjugate (do F0a and F0b)
    
end

%%% Return summed signal
F0 = sum(F([1 4],:),1)*1i;

%%% Construct Fn and Zn
idx=[fliplr(8:6:size(F,1)) 1 7:6:size(F,1)]; 
kvals = -kmax:kmax;
%%% Remove the lowest two negative states since these are never populated
%%% at echo time
idx(1:2)=[];
kvals(1:2)=[];

%%% Now reorder
FnA = F(idx,:);
%%% Conjugate
FnA(kvals<0,:)=conj(FnA(kvals<0,:));

%%% Repeat for Fnb - shift indices by 3
FnB = F(idx+3,:);
%%% Conjugate
FnB(kvals<0,:)=conj(FnB(kvals<0,:));

Fn = cat(3,FnA,FnB);

%%% Similar for Zn
ZnA = F(3:6:end,:);
ZnB = F(6:6:end,:);
Zn = cat(3,ZnA,ZnB);


    %%% NORMAL EPG transition matrix but replicated twice 
    % As per Weigel et al JMR 2010 276-285 
    function Tap = RF_rot(a,p)
        Tap = zeros([3 3]);
        Tap(1) = cos(a/2).^2;
        Tap(2) = exp(-2*1i*p)*(sin(a/2)).^2;
        Tap(3) = -0.5*1i*exp(-1i*p)*sin(a);
        Tap(4) = conj(Tap(2));
        Tap(5) = Tap(1);
        Tap(6) = 0.5*1i*exp(1i*p)*sin(a);
        Tap(7) = -1i*exp(1i*p)*sin(a);
        Tap(8) = 1i*exp(-1i*p)*sin(a);
        Tap(9) = cos(a);
        Tap = kron(eye(2),Tap);
    end

    function build_T(AA)
        ksft = 6*(6*(kmax+1)+1);
        for i2=1:36
            T(i1(i2):ksft:end)=AA(i2);
        end
    end

end
