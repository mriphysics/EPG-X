% Test 2: transient SPGR

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%%% White matter model (see Gloor 2008)
f = 0.1166;  %% F=0.132 = f/(1-f) => f=0.1166
kf = 4.3e-3;
kb = kf * (1-f)/f;
R1f = 1/779; % ms^-1
R1b = 1/779; %<- Gloor fix as T1f
R2f = 1/45;
R1obs = 0.5*(R1f + kf + R1b + kb)-0.5*sqrt((R1f + kf + R1b + kb).^2 ...
    -4*(R1f*R1b + R1f*kb + R1b*kf));


%%% Flip angle variation
ss = [0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15]/15;
ss = [ss fliplr(ss)];
s0 = zeros([1 16]);
alphas = [s0 s0 ss ss s0 s0 ss]; 
alphas=d2r(sin(alphas*pi/2)*40);

npulse = length(alphas);

%%% B1sqrdtau quoted from scanner
b1sqrdtau = 54.3 * alphas.^2; % ms uT^2
G = 15.1; % us

%%% Define a prep pulse and inversion delay
prep=struct;
prep.flip=pi;
prep.B1SqrdTau = 433.3;%<-- from scanner (min dur block pulse is 158)
prep.t_delay=0;

%%% Sequences
TR=12;
phi_spoiled = RF_phase_cycle(npulse,150); %<-Philips scanner
phi_bSSFP   = RF_phase_cycle(npulse,'balanced');

%%% Spoiled
% Standard EPG
[s0,Fn0,Zn0] = EPG_GRE(alphas,phi_spoiled,TR,1/R1obs,1/R2f,'prep',prep);
% EPGX MT
[s,Fn,Zn] = EPGX_GRE_MT(alphas,phi_spoiled,b1sqrdtau,TR,[1/R1f 1/R1b],1/R2f,f,kf,G,'prep',prep);

%%% balanced SSFP
% Standard EPG
[s1,Fn1,Zn1] = EPG_GRE(alphas,phi_bSSFP,TR,1/R1obs,1/R2f,'prep',prep,'kmax',inf);
% EPGX MT
[s2,Fn2,Zn2] = EPGX_GRE_MT(alphas,phi_bSSFP,b1sqrdtau,TR,[1/R1f 1/R1b],1/R2f,f,kf,G,'prep',prep,'kmax',inf);

%%% iFFT to get off-resonance profile
ssfp_mt = ifftshift(size(Fn2,1)*(ifft(ifftshift(Fn2,1),[],1)),1);
ssfp_nomt = ifftshift(size(Fn1,1)*(ifft(ifftshift(Fn1,1),[],1)),1);

%%
figure(1)
clf
subplot(3,1,1)
plot(r2d(alphas),'-')
ylabel('Flip angle (degrees)')
grid on
title('Flip Angles')
xlim([1 npulse])
subplot(3,1,2)
plot(abs(s0))
hold on
plot(abs(s)/(1-f),'-')
ylabel('$$\tilde{F}_0^a / M_0^a $$','Interpreter','latex')

legend('EPG','EPG-X (MT)','location','northwest')
grid on
title('Signal')

xlim([1 npulse])
subplot(3,1,3)
plot(squeeze(Zn0(1,:,1)))
hold on
plot(squeeze(Zn(1,:,1)))
plot(squeeze(Zn(1,:,2)))
xlabel('RF pulse number')
ylabel('$$\tilde{Z}_0^{a,b} / M_0$$','interpreter','latex')

legend('EPG: Z_0','EPG-X (MT): Z_0^a','EPG-X (MT): Z_0^b','location','southeast')
grid on
title('Longitudinal Magnetization, Z_0 state')
xlim([1 npulse])
gg=get(gcf,'children');
for ii=[2 4 5]
    gg(ii).XTick = 32:32:256;
end
setpospap([300    88   474   617])
set(gg([2 4 5]),'FontSize',13)
gg(2).Position = [0.1300    0.0665    0.7750    0.2494];
gg(4).Position = [0.1300    0.4118    0.7750    0.2134];
gg(5).Position = [0.1300    0.7114    0.7750    0.2134];
% print -dpng -r300 IRTFE_median_brain.png

%%% text
axes(gg(5))
text(-32,-5,'(a)','fontsize',17,'fontweight','bold')
text(-32,-60,'(b)','fontsize',17,'fontweight','bold')
text(-32,-125,'(c)','fontsize',17,'fontweight','bold')
print -dpng -r300 bin/figure4.png

%% Balanced display

phiTR= linspace(-pi,pi,size(ssfp_mt,1));
[~,idx0] = min(abs(phiTR));
[~,idx50] = min(abs(phiTR-pi/2));

figure(1);clf;
nr=2;nc=3;

subplot(nr,nc,2)
plot(abs(ssfp_nomt(idx0,:)))
hold on
plot(abs(ssfp_mt(idx0,:))/(1-f),'-')
ylabel('$$\tilde{F}_0^a / M_0^a $$','Interpreter','latex')

legend('EPG','EPG-X(MT)','location','northwest')
grid on
title('bSSFP signal \psi = 0')
set(gca,'XTick',0:32:256)
xlim([0 npulse])

subplot(nr,nc,5)
plot(abs(ssfp_nomt(idx50,:)))
hold on
plot(abs(ssfp_mt(idx50,:))/(1-f),'-')
xlabel('RF pulse number')
ylabel('$$\tilde{F}_0^a / M_0^a $$','Interpreter','latex')

legend('EPG','EPG-X(MT)','location','northwest')
grid on
title('bSSFP signal \psi = \pi/2')
set(gca,'XTick',0:32:256)
xlim([0 npulse])


subplot(nr,nc,[3 6])
plot(squeeze(Zn1(1,:,1)))
hold on
plot(squeeze(Zn2(1,:,1)))
plot(squeeze(Zn2(1,:,2)))
xlabel('RF pulse number')
ylabel('$$\tilde{Z}_0^{a,b} / M_0$$','interpreter','latex')
set(gca,'XTick',0:32:256)

legend('EPG: Z_0','EPG-X(MT): Z_0^a','EPG-X(MT): Z_0^b','location','southeast')
grid on
title('Longitudinal Magnetization, Z_0 state')
xlim([0 npulse])
set(gca,'XTick',0:32:256)

subplot(nr,nc,[1 4])
imagesc(1:npulse,phiTR,abs(ssfp_mt),[0 0.2])
colormap gray
hold
plot([1 npulse],[0 0],'--r','linewidth',2)
plot([1 npulse],[1 1]*pi/2,'--r','linewidth',2)
set(gca,'XTick',0:32:256)
title('Off-resonance sensitivity (with MT)')
xlabel('RF pulse number')
ylabel('Dephasing per TR period: \psi /rad')


setpospap([102         177        1072         370])

gg=get(gcf,'Children');
gap=0.075;ww=0.25;

gg(1).Position = [gap    0.2500    ww    0.5];
gg(3).Position = [3*gap+2*ww    0.1100    ww    0.8150];
gg(5).Position = [gap*2+ww   0.1100    ww    0.3412];
gg(7).Position = [gap*2+ww   0.5838    ww    0.3412];
set(gg([1 3 5 7]),'FontSize',13)

%%% text
axes(gg(3))
text(-700,-0.9,'(a)','fontsize',18,'fontweight','bold')
text(-380,-1.1,'(b)','fontsize',18,'fontweight','bold')
text(-48,-1.1,'(c)','fontsize',18,'fontweight','bold')

print -dpng -r300 bin/SF2.png
