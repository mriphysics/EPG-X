function [s,mxy,M] = isochromat_GRE(theta,phi,TR,T1,T2,Niso,varargin)
%   [F0,Fn,Zn,F] = isochromat_GRE(theta,phi,TR,T1,T2,Niso,varargin)
%
%   Single pool ISOCHROMAT simulation for gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               TR:         repetition time, ms
%               T1:         T1, ms
%               T2:         T2, ms
%               Niso:       number of isochromats, dephasing angles linearly spaced -pi to pi
%
%   optional arguments (use string then value as next argument)
%
%              diff:        structure with fields:
%                           G    - Gradient amplitude(s)
%                           tau  - Gradient durations(s)
%                           D    - Diffusion coeff m^2/s (i.e. expect 10^-9)
%
%   Outputs:                
%               s:          signal (ensemble average) directly after each
%                           excitation
%               ...
%
%
%   Shaihan Malik 2017-09-04
%                 2018-07-17: added diffusion handling  

%% Extra variables

diffusion_calc = false;

for ii=1:length(varargin)

    % Prep pulse - struct contains flip (rad), t_delay
    if strcmpi(varargin{ii},'prep')
        prep = varargin{ii+1};
    end
    
    % Diffusion - structure contains, G, tau, D - added 17/7/18
    if strcmpi(varargin{ii},'diff')
        diff = varargin{ii+1};
        diffusion_calc = true;
    end
    
    
end

%%% Fix number of TR periods
Npulse = length(theta);

%%% Isochromat phase distribution
psi = 2*pi*(0:fix(Niso)-1)/fix(Niso);

%%% Number of variables (each isochromat has Mx,My,Mz)
N = 3*Niso;

%%% Gradient dephasing matrices
rg={};
for jj=1:Niso
    rg{jj} = rotmat([0 0 psi(jj)]);
end
Rg = blkdiag(rg{:});
    
    
%%% Relaxation matrices 
E1=exp(-TR/T1);
E2=exp(-TR/T2);
E=eye(N);
ii = 1:(N/3);
E(3*N*(ii-1)+3*(ii-1)+1)=E2;
E(3*N*ii-2*N+3*(ii-1)+2)=E2;
E(3*N*ii-N+3*(ii-1)+3)=E1;

%%% composite matrix for rotation, dephasing
Reg = E*Rg;
Reg=sparse(Reg);


if diffusion_calc
    %%% Build a k-space filter function
    kk = -floor((Niso)/2):floor((Niso-1)/2); 

    kk = abs(kk); %<-- kk is the ORDER of the EPG state, so it must be positive (i.e. a single k accounts for the positive and negative states of that order) 
    [bDL, bDT] = EPG_diffusion_weights(diff.G,diff.tau,diff.D,kk);
    
    bDL=ifftshift(bDL);
    bDT=ifftshift(bDT);
end

%%% Now run the sim
M = zeros([N Npulse]);

xidx = 1:3:N; % indices of all Mx
yidx = 2:3:N; % indices of all My
zidx = 3:3:N; % indices of all Mz 

% Initialize Mz=1 for all isochromats
M0 = zeros([N 1]);
M0(zidx)=1;


%%% Initialize RF rotation matrix, which is modified but not re-declared
%%% each time
T=sparse(zeros([N N]));

% Initialise with thermal equilibrium
M(:,1) = M0; 
zf = (1-E1);% thermal recovery of Mz in TR period

% Loop over pulses: flip, record FID then dephase
for jj=1:Npulse
    
    % Get matrix ready for flip
    build_T_matrix_sub_implicit(rotmat(theta(jj)*[cos(phi(jj)) sin(phi(jj)) 0]));
    
    % Flip
    M(:,jj) = T*M(:,jj);
    
    % Dephase and T1 recovery
    if jj<Npulse
        M(:,jj+1) = Reg*M(:,jj);
        M(zidx,jj+1) = M(zidx,jj+1) + zf;
        
        if diffusion_calc
            
            % Now apply diffusion filtering
            mxy = M(xidx,jj+1)+1i*M(yidx,jj+1);
            mz = M(zidx,jj+1);
            
            % Apply
            mxyf = ifft(fft(mxy(:)).*bDT(:));
            mzf = ifft(fft(mz(:)).*bDL(:));
            
            % Put these back in the array
            M(xidx,jj+1)=real(mxyf);
            M(yidx,jj+1)=imag(mxyf);
            M(zidx,jj+1)=real(mzf);
        end
       
    
    end
end

% Now generate signal and demodulate it
mxy = M(xidx,:) + 1i*M(yidx,:);
% Get signal from mean
s = 1i*mean(mxy,1);
% demodulate this
s = s .* exp(-1i*phi(1:Npulse));


    function build_T_matrix_sub_implicit(AA)
        %%% This function operates on the existing T matrix, rather than
        %%% re-declare it each time. This is much faster 
        ix = 1:(N/3);
        T(3*N*(ix-1)+3*(ix-1)+1)=AA(1);
        T(3*N*(ix-1)+3*(ix-1)+2)=AA(2);
        T(3*N*(ix-1)+3*(ix-1)+3)=AA(3);
        T(3*N*ix-2*N+3*(ix-1)+1)=AA(4);
        T(3*N*ix-2*N+3*(ix-1)+2)=AA(5);
        T(3*N*ix-2*N+3*(ix-1)+3)=AA(6);
        T(3*N*ix-N+3*(ix-1)+1)=AA(7);
        T(3*N*ix-N+3*(ix-1)+2)=AA(8);
        T(3*N*ix-N+3*(ix-1)+3)=AA(9);
    end

    % Rotation matrix function
    function R = rotmat(u)
        
        % check zero input
        if any(u)
            th=norm(u);
            u=u/th;
            ct=cos(th);
            st=sin(th);
            R = [[ct + u(1)^2*(1-ct) u(1)*u(2)*(1-ct)-u(3)*st u(1)*u(3)*(1-ct)+u(2)*st];...
                [u(2)*u(1)*(1-ct)+u(3)*st ct+u(2)^2*(1-ct) u(2)*u(3)*(1-ct)-u(1)*st];
                [u(3)*u(1)*(1-ct)-u(2)*st u(2)*u(3)*(1-ct)+u(1)*st] ct+u(3)^2*(1-ct)];
            
        else
            R=eye(3);
        end
        
    end

    
end