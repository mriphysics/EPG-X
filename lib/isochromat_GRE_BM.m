function [s,mxy] = isochromat_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,Niso)
%   [F0,Fn,Zn,F] = isochromat_GRE_BM(theta,phi,TR,T1x,T2x,f,ka,Niso)
%
%   Isochromat based simualtion of Bloch McConnell coupled systems w/ gradient echo sequences
%
%   arguments:
%               theta:      vector of flip angles (rad) - length = #pulses
%               phi:        phase per pulse. This function can hence be
%                           used to simulate RF spoiling or balanced
%                           sequences depending on how phase is cycled
%                           see function RF_phase_cycle()
%               TR:         repetition time, ms
%               T1x:        [T1a T1b], ms
%               T2x:        [T2a T2b], ms
%               f:          fraction of compartment b
%               ka:         forward exchange rate from a->b (units ms^-1)
%               Niso:       Number of isochromats for simulation
%
%
%   Outputs:                
%
%
%   Shaihan Malik 2017-09-04


%% Set up variables


%%% Fix number of TR periods
Npulse = length(theta);

%%% Isochromat phase distribution
%psi = @(n)(2*pi*(0:fix(n)-1)/fix(n)); %<- distributes phases equivalently to FFT
psi = 2*pi*(0:fix(Niso)-1)/fix(Niso);

%%% Number of variables (each isochromat has Mxa,Mya,Mza,Mxb,Myb,Mzb)
N = 6*Niso;

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


%% Set up matrices for Relaxation and Exchange

%%% Relaxation-exchange matrix
La = diag(-[R2a R2a R1a]);
Lb = diag(-[R2b R2b R1b]);
Ka = diag([ka ka ka]);
Kb = diag([kb kb kb]);
A = [[La-Ka Kb];[Ka Lb-Kb]];
C = [0 0 R1a*M0a 0 0 R1b*M0b]';

Xi = expm(TR*A); %<-- operator for time period TR


% For inhomogeneous solution also need:
Zoff = (Xi - eye(6))*(A\C);
Zoff = repmat(Zoff(:),[Niso 1]); % replicate for all isochromats

% Multiply this up for all isochromats
Xi = kron(eye(Niso),Xi);
Xi = sparse(Xi);    

%%% Gradient dephasing matrices
rg={};
for jj=1:Niso
    rg{jj} = kron(eye(2),rotmat([0 0 psi(jj)])); %<- copy matrix for both pools
end
Rg = blkdiag(rg{:});
Rg = sparse(Rg);



%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);

%%% Matrix to hold all results
M = zeros([N Npulse]);

% Initialize Mz=1 for all isochromats
M0 = zeros([N 1]);
M0(3:6:N)=M0a;
M0(6:6:N)=M0b;
M(:,1) = M0;

%% Main body of gradient echo sequence, loop over TRs 

% Loop over pulses: flip, record FID then dephase
for jj=1:Npulse
    
    % Get matrix ready for flip
    build_T_matrix_sub_implicit(rotmat(theta(jj)*[cos(phi(jj)) sin(phi(jj)) 0]));
    
    % Flip
    M(:,jj) = T*M(:,jj);
    
    % Dephase and T1 recovery
    if jj<Npulse
        % First apply gradient dephasing
        M(:,jj+1) = Rg*M(:,jj);
        % now evolution due to relaxation and exchange
        M(:,jj+1) = Xi*M(:,jj+1)+Zoff;
    end
     
end


% Now generate signal and demodulate it
mxya = M(1:6:end,:) + 1i*M(2:6:end,:);
mxyb = M(4:6:end,:) + 1i*M(5:6:end,:);
% Get signal from mean
sa = 1i*mean(mxya,1);
sb = 1i*mean(mxyb,1);
% demodulate this
sa = sa .* exp(-1i*phi(1:Npulse));
sb = sb .* exp(-1i*phi(1:Npulse));

s=sa+sb;%return total signal
mxy=cat(3,mxya,mxyb);


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
