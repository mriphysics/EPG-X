%% Replicate result from Zhang et al:
%%% Zhang S, Liu Z, Grant A, Keupp J, Lenkinski RE, Vinogradov E. 
%%% Balanced Steady-State Free Precession (bSSFP) from an effective field 
%%% perspective: Application to the detection of chemical exchange (bSSFPX). 
%%% J. Magn. Reson. 2017;275:55?67.

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));


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


%%%% Parameters to loop over
d_ppm = [0 1 3 9];
k_scaling = [0.5 1 2 3 5]; % k_ba is this factor multiplied by delta in Hz

n1 = length(d_ppm);
n2 = length(k_scaling);

%%% Variable to store results
S={};

%%% Loop over the parameters
for ii=1:n1
    d_kHz = d_ppm(ii)*1e-6 * 3 * 42.58e3; %<- 3T, 42.58kHz/T
    
    for jj=1:n2
        
        k_Hz_ba = k_scaling(jj) * d_kHz*1e3;% solute to water exchange rate specified in paper
        
        k_Hz_ab = k_Hz_ba * f/(1-f); % EPG-X requires water to solute rate
        
 
        tic
        [sx,Fnx,Znx] = EPGX_GRE_BM(alpha,phi,TR,T1x,T2x,...
            f,k_Hz_ab*1e-3,'delta',d_kHz,'kmax',50);%<-- flip is same for both pools (v short pulses)
        tt1=toc;

        fprintf(1,'%d %d EPGX = %1.3f s \n',ii,jj,tt1);
        M = fft(ifftshift(Fnx,1),[],1);
        M = cat(1,M,M,M);
        
        S{ii,jj} = squeeze(M(:,end,:));
        
    end
end

    
%% Figure with asymmetry
figure(22)
clf
pp = linspace(-3,3,size(M,1));
pp = pp / (2*TR*1e-3);

for ii=1:3
    subplot(2,3,ii)
    hold on
    plot(pp,abs(sum(S{1,1},2)),'LineWidth',2)
    for jj=1:n2
        plot(pp,abs(sum(S{ii+1,jj},2)),'LineWidth',2)
    end
    grid on
    xlim([-500 500])
    title(sprintf('EPG-X %d ppm',d_ppm(ii+1)))
    set(gca,'XTick',-500:250:500)
end
    
xl = {[-250 0],[-500 -250],[-250 0]};
yl = {[0 0.8],[-1 0.8],[0 0.8]};
pidx={};
for ii=1:3
    pidx{ii} = (pp<xl{ii}(1))|(pp>xl{ii}(2));
end
%%% Asymmetry
for ii=1:3
    subplot(2,3,ii+3)
    hold on
    mxy=abs(sum(S{1,1},2));
    mxyneg=circshift(flipud(mxy),[1 0]);
    mxy_asymm=(mxyneg-mxy)./mxyneg;
    mxy_asymm(pidx{ii})=0;
    plot(pp,mxy_asymm,'LineWidth',2)
    for jj=1:n2
        mxy=abs(sum(S{ii+1,jj},2));
        mxyneg=circshift(flipud(mxy),[1 0]);
        mxy_asymm=(mxyneg-mxy)./mxyneg;
        mxy_asymm(pidx{ii})=0;
        plot(pp,mxy_asymm,'LineWidth',2)
    end
    grid on
    xlim([-500 500])
    ylim(yl{ii})
    title(sprintf('bSSFPX assymetry %d ppm',d_ppm(ii+1)))
    set(gca,'XTick',-500:250:500)
end
setpospap([100 300 1000 400])
print -dpng -r300 bin/bSSFPXfigure.png


%% 2021-3-11: think about a bSSFP sequence with variable flips on each side

TR = 2.025;
alpha = 10;
npulse = 1500;
phi = RF_phase_cycle(npulse,'balanced');
alpha = d2r(alpha)*ones(npulse,1);

% different flip and phase for CSET pool
alpha2 = 5;
phi2 = 0*phi;

% Run sims
IX = 
d_kHz = d_ppm(ii)*1e-6 * 3 * 42.58e3; %<- 3T, 42.58kHz/T
    
    for jj=1:n2
        
        k_Hz_ba = k_scaling(jj) * d_kHz*1e3;% solute to water exchange rate specified in paper
        
        k_Hz_ab = k_Hz_ba * f/(1-f); % EPG-X requires water to solute rate
        
 
        tic
        [sx,Fnx,Znx] = EPGX_GRE_BM(alpha,phi,TR,T1x,T2x,...
            f,k_Hz_ab*1e-3,'delta',d_kHz,'kmax',50);%<-- flip is same for both pools (v short pulses)
        tt1=toc;

        fprintf(1,'%d %d EPGX = %1.3f s \n',ii,jj,tt1);
        M = fft(ifftshift(Fnx,1),[],1);
        M = cat(1,M,M,M);