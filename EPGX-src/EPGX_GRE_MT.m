function [F0,Fn,Zn,F] = EPGX_GRE_MT(theta,phi,B1SqrdTau,TR,T1x,T2a,f,ka,G,varargin)
%   [F0,Fn,Zn,F] = EPGX_GRE_MT(theta,phi,B1SqrdTau,TR,T1x,T2a,f,ka,G,varargin)
%
%   EPG-X for Pulsed MT systems w/ gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               B1SqrdTau:  vector of RF pulse energies, uT^2 ms
%               TR:         repetition time, ms
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
%
%               prep:       can be used to simulate prep pulse prior to
%                           gradient echo. Assume all transverse
%                           magnetization after prep is spoiled.
%                           structure with fields:
%                           flip    -   flip angle, rad
%                           t_delay -   time delay, ms
%                           B1SqrdTau - uT^2 ms (for RF saturation)
%
%   Outputs:                
%               F0:         signal (F0 state) directly after each
%                           excitation
%               Fn:         full EPG diagram for all transverse states
%               Zn:         full EPG diagram for all longitudinal states
%               F:          full state matrix. Each column is arranged as
%                           [F0a F0a* Z0a Z0b F1a F-1a* Z1a Z1b F2a F-2a* Z2a Z2b ...] etc
%
%   Shaihan Malik 2017-07-20


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
    
    % Prep pulse - struct contains flip (rad), t_delay
    if strcmpi(varargin{ii},'prep')
        prep = varargin{ii+1};
    end
end

%%% The maximum order varies through the sequence. This can be used to speed up the calculation    
np = length(theta);
% if not defined, assume want max
if ~exist('kmax','var')
    kmax = np - 1;
end

if isinf(kmax)
    % this flags that we don't want any pruning of pathways
    allpathways = true;
    kmax = np - 1; % this is maximum value
else
    allpathways = false;
end

%%% Variable pathways
if allpathways
    kmax_per_pulse = 0:kmax;
else
    kmax_per_pulse = [1:ceil(np/2) (floor(np/2)):-1:1];
    kmax_per_pulse(kmax_per_pulse>kmax)=kmax;
    
    if max(kmax_per_pulse)<kmax
        kmax = max(kmax_per_pulse);
    end
end

%%% Number of states is 4x(kmax +1) -- +1 for the zero order
N=4*(kmax+1);

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
Xi_T = expm(TR*Lambda_T); %<-- operator for time period TR

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1a-ka kb];[ka -R1b-kb]];
Xi_L = expm(TR*Lambda_L); 

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


%% F matrix (many elements zero, not efficient)
F = zeros([N np]); %%<-- records the state after each RF pulse 

%%% Initial State
FF = zeros([N 1]);
FF(3)=1-f;   %Z0a
FF(4)=f;     %Z0b


%% Prep pulse - execute here
if exist('prep','var')
    %%% Assume the prep pulse leaves NO transverse magnetization, or that
    %%% this is spoiled so that it cannot be refocused. Only consider
    %%% z-terms
    
    % RF rotation just cos(theta) on both Mz terms
    prep.WT = pi*gam^2*prep.B1SqrdTau*G; 

    R = diag([cos(prep.flip) exp(-prep.WT)]);
    zidx = [3 4];
    FF(zidx)=R*FF(zidx);
    
    % Now apply time evolution during delay period
    Xi_L_prep = expm(prep.t_delay*Lambda_L); 
    Zoff_prep = (Xi_L_prep - eye(2))*(Lambda_L\C_L);
    FF([3 4]) = Xi_L_prep * FF([3 4]) + Zoff_prep;
    
end

%% Main body of gradient echo sequence, loop over TRs 

for jj=1:np 
    %%% RF transition matrix
    A = RF_rot_sat(theta(jj),phi(jj),WT(jj));
   
    %%% Variable order of EPG, speed up calculation
    kmax_current = kmax_per_pulse(jj);
    kidx = 1:4*(kmax_current+1); %+1 because states start at zero
    
    %%% Replicate A to make large transition matrix
    build_T(A);
    
    %%% Apply flip and store this
    F(kidx,jj)=T(kidx,kidx)*FF(kidx);
    
    if jj==np
        break
    end
    
    %%% Now deal with evolution
    FF(kidx) = XS(kidx,kidx)*F(kidx,jj)+b(kidx);
    
    % Deal with complex conjugate after shift
    FF(1)=conj(FF(1)); %<---- F0 comes from F-1 so conjugate 
end


%%% Return FID (F0) signal
F0=F(1,:);

%%% phase demodulate
F0 = F0(:) .* exp(-1i*phi(:)) *1i;

%%% Construct Fn and Zn
idx=[fliplr(6:4:size(F,1)) 1 5:4:size(F,1)]; 
kvals = -kmax:kmax;

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
