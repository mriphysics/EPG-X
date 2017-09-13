function [s,mxy] = isochromat_GRE_MT(theta,phi,B1SqrdTau,TR,T1x,T2a,f,ka,G,Niso) 
%   [F0,Fn,Zn,F] = isochromat_GRE_MT(theta,phi,B1SqrdTau,TR,T1x,T2a,f,ka,G,Niso) 
%   Isochromat based simualtion of Pulsed MT systems w/ gradient echo sequences
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
%
%   Outputs:                
%
%
%   Shaihan Malik 2017-09-05


%% Set up variables


%%% Fix number of TR periods
Npulse = length(theta);

%%% Isochromat phase distribution
%psi = @(n)(2*pi*(0:fix(n)-1)/fix(n)); %<- distributes phases equivalently to FFT
psi = 2*pi*(0:fix(Niso)-1)/fix(Niso);

%%% Number of variables (each isochromat has Mxa,Mya,Mza,Mzb)
N = 4*Niso;

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

% Multiply this up for all isochromats
Xi = kron(eye(Niso),Xi);
Xi = sparse(Xi);    

Zoff = repmat([0;0;Zoff(:)],[Niso 1]); % replicate for all isochromats

%%% Gradient dephasing matrices
rg={};
for jj=1:Niso
    rg{jj} = rotmat([0 0 psi(jj)]);
    rg{jj}(4,4)=1;% no z-rotation in either pool from gradient
end
Rg = blkdiag(rg{:});
Rg = sparse(Rg);


%%% Pre-allocate RF matrix
T = zeros(N,N);
T = sparse(T);
% store the indices of the top 4x4 corner, this helps build_T
i1 = [];
for ii=1:4
    i1 = cat(2,i1,sub2ind(size(T),1:4,ii*ones(1,4)));
end

%%% Compute saturation terms for RF pulse
G = G *1e-3; %<- convert from us to ms
gam = 267.5221 *1e-3; %< rad /ms /uT
WT = pi*gam^2*B1SqrdTau*G; 


%%% Matrix to hold all results
M = zeros([N Npulse]);

% Initialize Mz=1 for all isochromats
M0 = zeros([N 1]);
M0(3:4:N)=M0a;
M0(4:4:N)=M0b;
M(:,1) = M0;

%% Main body of gradient echo sequence, loop over TRs 

% Loop over pulses: flip, record FID then dephase
for jj=1:Npulse
    
    % Get matrix ready for flip
    rfmat=rotmat(theta(jj)*[cos(phi(jj)) sin(phi(jj)) 0]);
    rfmat(4,4)=exp(-WT(jj));
    build_T(rfmat);
    
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
mxy = M(1:4:end,:) + 1i*M(2:4:end,:);

% Get signal from mean
s = 1i*mean(mxy,1);
% demodulate this
s = s .* exp(-1i*phi(1:Npulse));


    function build_T(AA)
        %ksft = 4*(4*(kmax+1)+1);
        ksft=4*(4*Niso+1);
        for i2=1:16
            T(i1(i2):ksft:end)=AA(i2);
        end
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
