%% Test 1: Steady state gradient echo vs EPG

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%% SPGR (RF spoiling)
% Use simple parameters

%%% Sequences
TR = 10;
alpha = 15;
phi0 = 117;


%%% Relaxation parameters: single pool
T1=1000;
T2=100;

%%% Relaxation parameters: MT
T1free = 1000;
T1bound = 1000;
f = 0.2; %<-- bound fraction
k = 10e-3; % Exchange rate from free to bound (large to small)

%%% Relaxation parameters: ChX
T1x = [1000 500];
T2x = [100 20];    
% exchange parameters are as above

%%% RF saturation factor for MT
G = 15.1;% us 
b1 = 13; % uT
gam = 267.5221 *1e-3; %< rad /ms /uT
trf = d2r(alpha)/(gam*b1);% ms
b1sqrdtau = b1^2*trf;

% EPG simulations
npulse=200;
phi = RF_phase_cycle(npulse,phi0);

% single pool
[s0,Fn0,Zn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2);

% MT
[smt,Fnmt,Znmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,[T1free T1bound],T2,f,k,G);

% BM
[sx,Fnx,Znx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,f,k);

ss0 = ssSPGR(d2r(alpha),TR,T1);
ssx = ssSPGR_BM(d2r(alpha),TR,T1x,f,k);
ssmt= ssSPGR_MT(d2r(alpha),b1sqrdtau,TR,[T1free T1bound],f,k,G);



%% Compute whole spoiling curve
phi_arr = 1:0.5:180;
nphi=length(phi_arr);
npulse = floor(5*T1/TR);

AA = d2r(alpha)*ones(npulse,1);

if 0 % load stored result if 0
Sig = zeros(nphi,3);
figure(1)
clf
for ii=1:nphi
   
    % Compute RF spoling phase cycles
    phi = RF_phase_cycle(npulse,phi_arr(ii));

    % single pool
    tmp1 = EPG_GRE(AA,phi,TR,T1,T2);

    % MT
    tmp2 = EPGX_GRE_MT(AA,phi,b1sqrdtau*ones(npulse,1),...
    TR,[T1free T1bound],T2,f,k,G);

    % BM
    tmp3 = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,f,k);

    Sig(ii,1) = abs(tmp1(end));
    Sig(ii,2) = abs(tmp2(end));
    Sig(ii,3) = abs(tmp3(end));
    
    disp([ii nphi])
    %save Sig Sig
    plot(Sig,'-')
    drawnow
    pause(0.0001)
    
end
    save bin/SigAll Sig
else
    load bin/SigAll
end


%% Combined figure
figure(1)
clf
subplot(211)
hold on
p1=plot(abs(s0),'linewidth',1);
p2=plot(abs(smt),'linewidth',1);
p3=plot(abs(sx),'linewidth',1);
plot(abs(ss0)*ones(size(s0)),'-.','Color',p1.Color* 0.6)
plot(abs(ssmt)*ones(size(s0)),'-.','Color',p2.Color* 0.6)
plot(abs(ssx)*ones(size(s0)),'-.','Color',p3.Color* 0.6)
grid on
xlabel('TR number')
ylabel('Signal / M_0')
title('SPGR approach to steady-state for \Phi_0=117')
legend('EPG','EPG-X(MT)','EPG-X(BM)','Steady-state (Ernst)',...
    'Steady-state MT','Steady-state BM','location','northeast')
set(gca,'fontsize',13)

subplot(212)
plot(phi_arr,Sig)
grid on
hold on
xlim([0 180])
xlabel('RF spoil phase')
ylabel('Signal / M_0')
title('SPGR steady-state dependence on \Phi_0')
set(gca,'fontsize',13)

%%% Add steady state solutions
plot(abs(ss0)*ones(size(s0)),'-.','Color',p1.Color * 0.6)
plot(abs(ssmt)*ones(size(s0)),'-.','Color',p2.Color* 0.6)
plot(abs(ssx)*ones(size(s0)),'-.','Color',p3.Color* 0.6)

text(-25,0.12,'(a)','fontsize',20,'fontweight','bold')
text(-25,0.03,'(b)','fontsize',20,'fontweight','bold')

setpospap([360   174   457   524])

print -dpng -r300 bin/Figure2.png

%% bSSFP with same parameters as above

% EPG simulations
npulse=5*T1/TR;
phi = RF_phase_cycle(npulse,'balanced');

% single pool
[~,fn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2,'kmax',inf);
mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);

% MT
[~,fnmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,[T1free T1bound],T2,f,k,G,'kmax',inf);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);

% BM
[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,f,k,'kmax',inf);
mxyx = size(fnmt,1)*ifftshift(ifft(ifftshift(fnx,1),[],1),1);
mxyx = sum(mxyx,3); % sum compartments

dw=0;% on resonance 
ss0 = ssSSFP(d2r(alpha),TR,dw,T1,T2);
ssx = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,f,k);
ssmt= ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,dw,[T1free T1bound],T2,f,k,G);


%%% steady state
phiTR = linspace(-pi/TR,pi/TR,npulse);
ssfpxSS = zeros(size(phiTR));
ssfpmtSS= zeros(size(phiTR));
ssfpss = zeros(size(phiTR));

for ii=1:npulse
    ssfpss(ii) = abs(ssSSFP(d2r(alpha),TR,phiTR(ii),T1,T2));
    ssfpxSS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,f,k));
    ssfpmtSS(ii) = abs(ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,phiTR(ii),[T1free T1bound],T2,f,k,G));
end
% analytic form
mxyGloor = ssSSFP_Gloor(d2r(alpha),b1sqrdtau, TR, [T1free T1bound],T2,f,k,G);

%%
%
figure(6)
clf

subplot(321)
imagesc(1:npulse,phiTR*TR,abs(mxys),[0 0.2])
title('M_+(\psi,t): EPG')
ylabel('\psi / rad')
xlabel('TR number')

subplot(323)
imagesc(1:npulse,phiTR*TR,abs(mxyx),[0 0.2])
title('M_+(\psi,t): EPG-X (BM)')
xlabel('TR number')
ylabel('\psi / rad')

subplot(325)
imagesc(1:npulse,phiTR*TR,abs(mxymt),[0 0.2])
title('M_+(\psi,t): EPG-X (MT)')
xlabel('TR number')
ylabel('\psi / rad')

colormap gray

subplot(322)
pp=plot(linspace(-pi,pi,npulse),ssfpss,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxys,1)),squeeze(abs(mxys(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('Dephasing / TR (\psi / rad)')
legend( 'Direct Steady-state','EPG (TR 500)','location','south')
title('bSSFP single compartment')


subplot(324)
pp=plot(linspace(-pi,pi,npulse),ssfpxSS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxyx,1)),squeeze(abs(mxyx(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('Dephasing / TR (\psi / rad)')
legend( 'Direct Steady-state','EPG-X (BM) (TR 500)','location','south')
title('bSSFP: two compt with BM')

subplot(326)
pp=plot(linspace(-pi,pi,npulse),ssfpmtSS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxymt,1)),squeeze(abs(mxymt(:,end))),'--','linewidth',2)
plot(0,mxyGloor,'*r','MarkerSize',12,'Linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('Dephasing / TR (\psi / rad)')
legend( 'Direct Steady-state','EPG-X (MT) (TR 500)','Analytic (Gloor 2008)','location','north')
title('bSSFP: two compt with MT')

setpospap([300 100 500 600])
print -dpng -r300 bin/Figure3.png

%% Additional figure - approach to ss with isochromats

%%% SPGR (RF spoiling)
% Use simple parameters

%%% Sequences
TR = 10;
alpha = 15;
phi0 = 117;


%%% Relaxation parameters: single pool
T1=1000;
T2=100;

%%% Relaxation parameters: MT
T1free = 1000;
T1bound = 1000;
f = 0.2; %<-- bound fraction
k = 10e-3; % Exchange rate from free to bound (large to small)

%%% Relaxation parameters: ChX
T1x = [1000 500];
T2x = [100 20];    
% exchange parameters are as above

%%% RF saturation factor for MT
G = 15.1;% us 
b1 = 13; % uT
gam = 267.5221 *1e-3; %< rad /ms /uT
trf = d2r(alpha)/(gam*b1);% ms
b1sqrdtau = b1^2*trf;


%%% Simulation parameters as above
npulse=200;
phi = RF_phase_cycle(npulse,phi0);

% single pool EPG
[s0,Fn0,Zn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2);
% EPG-X(BM)
[sx,Fnx,Znx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,f,k);
% EPG-X(MT)
[smt,Fnmt,Znmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,[T1free T1bound],T2,f,k,G);


% Isochromat simulations
Niso = [10:10:200 300 400 1000];
Siso = {};
rmse = zeros(length(Niso),3);

for jj=1:length(Niso)
    
    % Single pool
    siso = isochromat_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2,Niso(jj));
    Siso{jj,1} = siso;
    rmse(jj,1) = norm(abs(siso(:))-abs(s0(:)))/norm(s0(:));
    
    % BM
    siso = isochromat_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,f,k,Niso(jj));
    Siso{jj,2} = siso;
    rmse(jj,2) = norm(abs(siso(:))-abs(sx(:)))/norm(sx(:));
    
    % MT
    siso = isochromat_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1)...
        ,TR,[T1free T1bound],T2,f,k,G,Niso(jj));
    Siso{jj,3} = siso;
    rmse(jj,3) = norm(abs(siso(:))-abs(smt(:)))/norm(smt(:));
end

%% plot these
figfp(1)
nr=3;nc=2;

sig_epg = {s0,sx,smt};
titls = {'Single Pool','BM','MT'}
for ii=1:3
    subplot(nr,nc,ii*2-1)
    
    p1=plot(abs(sig_epg{ii}),'r');
    hold
    pp=[];
    for jj=1:length(Niso)
        pp(jj)=plot(abs(Siso{jj,ii}),'^-','color',[1 1 1]*(0.3*(length(Niso)-jj)/length(Niso)+0.6),'linewidth',0.1);
    end
    set(pp,'markersize',2)
    legend('EPG-X','isochromat predictions')
    ylim([0 0.3])
    uistack(p1,'top');
    
    
    grid on
    xlabel('RF pulse number')
    ylabel('Predicted signal')
    title(titls{ii})
    
    subplot(nr,nc,ii*2)
    semilogy(Niso,rmse(:,ii))
    grid on
    xlabel('Number of isochromats')
    ylabel('RMS difference')
    title(sprintf('Difference between methods, %s',titls{ii}))
    ylim([10^-16 1])
end
