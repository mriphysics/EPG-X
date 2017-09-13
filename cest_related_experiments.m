%% Look at effect on SPGR of varying the TR with BM model, such that the TR is a multiple of exchange time

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%% SPGR (RF spoiling)
% Use simple parameters

%%% Sequences
alpha = 15;
phi0 = 117;

%%% Relaxation parameters: ChX
T1x = [1000 500];
T2x = [100 20];    
f = 0.2; %<-- small pool fraction
k = 50e-3; % Exchange rate from free to bound (large to small)


% EPG simulations
npulse=200;
phi = RF_phase_cycle(npulse,phi0);

% range of TRs
tr = 1:100;
ntr=length(tr);
S = zeros([npulse ntr 2]);

for ii=1:ntr
    [sx,Fnx,Znx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,tr(ii),T1x,T2x,f,k,'offres',1);
    S(:,ii,1) = sx(:);
    [sx,Fnx,Znx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,tr(ii),T1x,T2x,f,k,'offres',0.4,'delta',-0.4);
    S(:,ii,2) = sx(:);
end
%%
figfp(1)
for ii=1:2
subplot(1,3,ii)
imagesc(tr,1:200,abs(S(:,:,ii)),[0 0.2])
colorbar
end
subplot(1,3,3)
imagesc(tr,1:200,abs(S(:,:,1))-abs(S(:,:,2)),[-0.02 0.05])
colorbar

%% Look at off resonance profile

T1x = [2000 1000];
T2x = [50 1000];
f = 0.05/1.05;


TR = 2.025;

d_ppm = 1; %delta in ppm
d_kHz = d_ppm*1e-6 * 3 * 42e3; %<- 3T, 42kHz/T

k_Hz_ba = 1. * d_kHz*1e3;% solute to water exchange rate specified in paper

% k_Hz_ba = 0.02;
% d_kHz = 0.02;

k_Hz_ab = k_Hz_ba * f/(1-f); % EPG-X requires water to solute rate

alpha = 10;
npulse = 2000;
phi = RF_phase_cycle(npulse,'balanced');
alpha = d2r(alpha)*ones(npulse,1);

d = struct;
d.G = [-10 5 -10]; % mT/m 
d.tau = [0.5 1 0.5]; %ms %balanced readout
d.D = 2.3e-9;

[sx,Fnx,Znx] = EPGX_GRE_BM(alpha,phi,TR,T1x,T2x,...
    f,k_Hz_ab*1e-3,'offres',1,'delta',d_kHz,'kmax',250,'diff',d);

Niso=101;

[si,mxy] = isochromat_GRE_BM(alpha,phi,TR,T1x,T2x,...
    f,k_Hz_ab*1e-3,Niso,'delta',-d_kHz);

%M = fftshift(ifft(ifftshift(Fnx,1),[],1),1);
M = size(Fnx,1)*ifft(ifftshift(Fnx,1),[],1);
M = cat(1,M,M,M);
pp = linspace(-3,3,size(M,1));
pp2 = 2*pi*(0:fix(Niso)-1)/fix(Niso)-pi;
figfp(1)
plot(pp,squeeze(abs(M(:,end,1))))
hold
plot(pp2/pi,squeeze(abs(mxy(:,end,1))))
grid on
xlim([-2 2])

%% Generate some cases from the Zhang paper

%%% Parameters that are fixed
T1x = [2000 1000];
T2x = [50 1000];
f = 0.05/1.05;

%%% Sequence parameters
TR = 2.025;
alpha = 10;
npulse = 1500;
phi = RF_phase_cycle(npulse,'balanced');
alpha = d2r(alpha)*ones(npulse,1);
d = struct;
d.G = [-10 5 -10]; % mT/m 
d.tau = [0.5 1 0.5]; %ms %balanced readout
d.D = 2.3e-9;
tau_rf = 25e-6;%<- RF duration

%%%% Parameters to loop over
d_ppm = [0 1 3 9];
k_scaling = [0.5 1 2 3 5]; % k_ba is this factor multiplied by delta in Hz

n1 = length(d_ppm);
n2 = length(k_scaling);

%%% Variable to store results
S={};
Siso={};

%%% Loop over the parameters
for ii=1:n1
    d_kHz = d_ppm(ii)*1e-6 * 3 * 42.6e3; %<- 3T, 42kHz/T
    
    for jj=1:n2
        
        k_Hz_ba = k_scaling(jj) * d_kHz*1e3;% solute to water exchange rate specified in paper
        
        k_Hz_ab = k_Hz_ba * f/(1-f); % EPG-X requires water to solute rate
        
        
        
%         [sx,Fnx,Znx] = EPGX_GRE_BM(alpha,phi,TR,T1x,T2x,...
%             f,k_Hz_ab*1e-3,'offres',sinc(tau_rf*d_kHz*1e3),'delta',-d_kHz,'kmax',50,'diff',d);%<-- flip is same for both pools (v short pulses)
        tic
        [sx,Fnx,Znx] = EPGX_GRE_BM(alpha,phi,TR,T1x,T2x,...
            f,k_Hz_ab*1e-3,'delta',-d_kHz,'kmax',50,'diff',d);%<-- flip is same for both pools (v short pulses)
        tt1=toc;
        tic
        [~,mxy] = isochromat_GRE_BM(alpha,phi,TR,T1x,T2x,...
            f,k_Hz_ab*1e-3,101,'delta',-d_kHz);
        tt2=toc;
        fprintf(1,'%d %d EPGX = %1.3f s iso = %1.3f s\n',ii,jj,tt1,tt2);
        M = size(Fnx,1)*ifft(ifftshift(Fnx,1),[],1);
        M = cat(1,M,M,M);
        
        S{ii,jj} = squeeze(M(:,end,:));
        
        mxy = cat(1,mxy,mxy,mxy);
        Siso{ii,jj} = squeeze(mxy(:,end,:));
        
    end
end

%% plots
figfp(21)

pp = linspace(-3,3,size(M,1));
pp = pp / (2*TR*1e-3);

for ii=1:3
    subplot(2,3,ii)
    hold on
    plot(pp,abs(S{1,1}(:,1)))
    for jj=1:n2
        plot(pp,abs(S{ii+1,jj}(:,1)))
        %plot(pp,abs(sum(S{ii+1,jj},2)))
    end
    grid on
   xlim([-500 500])
   title(sprintf('EPG-X %d ppm',d_ppm(ii+1)))
end
    
    
for ii=1:3
    subplot(2,3,ii+3)
    hold on
    plot(pp,abs(Siso{1,1}(:,1)))
    for jj=1:n2
        plot(pp,abs(Siso{ii+1,jj}(:,1)))
        %plot(pp,abs(sum(Siso{ii+1,jj},2)))
    end
    grid on
    xlim([-500 500])
    title(sprintf('Isochromat %d ppm',d_ppm(ii+1)))
end
    