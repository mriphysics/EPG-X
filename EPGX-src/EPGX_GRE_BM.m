function [F0,Fn,Zn,F] = EPGX_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,varargin)
%   [F0,Fn,Zn,F] = EPGX_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,varargin)
%
%
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

%%% Number of states is 6x(kmax +1) -- +1 for the zero order
N=6*(kmax+1);

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
Lambda_T = diag([-R2a-ka -R2a-ka -R2b-kb -R2b-kb]);
Lambda_T(1,3) = kb;
Lambda_T(2,4) = kb;
Lambda_T(3,1) = ka;
Lambda_T(4,2) = ka;
Xi_T = expm(TR*Lambda_T); %<-- operator for time period TR

%%% Relaxation-exchange matrix for Longitudinal components
Lambda_L = [[-R1a-ka kb];[ka -R1b-kb]];
Xi_L = expm(TR*Lambda_L); 

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
XS=Xi*S;
XS=sparse(XS);

%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

% store the indices of the top 6x6 corner, this helps build_T
i1 = [];
for ii=1:6
    i1 = cat(2,i1,sub2ind(size(T),1:6,ii*ones(1,6)));
end


%% F matrix (many elements zero, use sparse represenatation)
F = zeros([N np]); %%<-- records the state after each RF pulse 

%%% Initial State
FF = zeros([N 1]);
FF(3)=1-f;   %Z0a
FF(6)=f;     %Z0b


%% Prep pulse - execute here
if exist('prep','var')
    %%% Assume the prep pulse leaves NO transverse magnetization, or that
    %%% this is spoiled so that it cannot be refocused. Only consider
    %%% z-terms
    
    % RF rotation just cos(theta) on both Mz terms
    R = diag([cos(prep.flip) cos(prep.flip)]);
    zidx = [3 6];
    FF(zidx)=R*FF(zidx);
    
    % Now apply time evolution during delay period
    Xi_L_prep = expm(prep.t_delay*Lambda_L); 
    Zoff_prep = (Xi_L_prep - eye(2))*(Lambda_L\C_L);
    FF([3 6]) = Xi_L_prep * FF([3 6]) + Zoff_prep;
    
end

%% Main body of gradient echo sequence, loop over TRs 

for jj=1:np 
    %%% RF transition matrix
    A = RF_rot(theta(jj),phi(jj));
   
    %%% Variable order of EPG, speed up calculation
    kmax_current = kmax_per_pulse(jj);
    kidx = 1:6*(kmax_current+1); %+1 because states start at zero
    
    %%% Replicate A to make large transition matrix
    build_T(A);
    
    %%% Apply flip and store this: splitting these large matrix
    %%% multiplications into smaller ones might help
    F(kidx,jj)=T(kidx,kidx)*FF(kidx);
    
    if jj==np
        break
    end
    
    %%% Now deal with evolution
    FF(kidx) = XS(kidx,kidx)*F(kidx,jj)+b(kidx);
    
    % Deal with complex conjugate after shift
    FF([1 4])=conj(FF([1 4])); %<---- F0 comes from F-1 so conjugate (do F0a and F0b)
end


%%% Return summed signal
F0=sum(F([1 4],:),1);

%%% phase demodulate
F0 = F0(:) .* exp(-1i*phi(:)) *1i;


%%% Construct Fn and Zn
idx=[fliplr(8:6:size(F,1)) 1 7:6:size(F,1)]; 
kvals = -kmax:kmax;

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