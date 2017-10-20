%% Test 1: Steady state gradient echo vs EPG

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%% SPGR (RF spoiling)
% Use simple parameters

%%% Sequences
TR = 5;
alpha = 10;
phi0 = 117;


%%% Relaxation parameters: single pool, copy MT model
T1=779;
T2=45;

%%% Relaxation parameters: MT
T1_MT = [779 779];
f_MT = 0.117;
k_MT = 4.3e-3;
T2_MT = 45;

%%% Relaxation parameters: ChX
T1x = [1000 500];
T2x = [100 20];    
kx = 2e-3;
fx = 0.2;

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
    TR,T1_MT,T2_MT,f_MT,k_MT,G);

% BM
[sx,Fnx,Znx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx);

ss0 = ssSPGR(d2r(alpha),TR,T1);
ssx = ssSPGR_BM(d2r(alpha),TR,T1x,fx,kx);
ssmt= ssSPGR_MT(d2r(alpha),b1sqrdtau,TR,T1_MT,f_MT,k_MT,G);



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
    TR,T1_MT,T2_MT,f_MT,k_MT,G);

    % BM
    tmp3 = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,kx);

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
p1=plot(abs(s0),'linewidth',1.5);
p2=plot(abs(smt),'linewidth',1.5);
p3=plot(abs(sx),'linewidth',1.5);

ylim([0 0.175])
grid on
xlabel('TR number')
ylabel('Signal / M_0')
title('SPGR approach to steady-state for \Phi_0=117')
% legend('EPG','EPG-X(MT)','EPG-X(BM)','Steady-state (Ernst)',...
%     'Steady-state MT','Steady-state BM','location','northeast')
legend('EPG','EPG-X(MT)','EPG-X(BM)')
set(gca,'fontsize',13)

%%% Add steady state solutions
[xx,yy]=ds2nfu([0 10],abs(ss0)*[1 1]);
aa1 = annotation('arrow',xx,yy,'Color',p1.Color);
[xx,yy]=ds2nfu([210 201],abs(ss0)*[1 1]);
aa2 = annotation('arrow',xx,yy,'Color',p1.Color);
[xx,yy]=ds2nfu([0 10],abs(ssmt)*[1 1]);
aa3 = annotation('arrow',xx,yy,'Color',p2.Color);
[xx,yy]=ds2nfu([210 201],abs(ssmt)*[1 1]);
aa4 = annotation('arrow',xx,yy,'Color',p2.Color);
[xx,yy]=ds2nfu([0 10],abs(ssx)*[1 1]);
aa5 = annotation('arrow',xx,yy,'Color',p3.Color);
[xx,yy]=ds2nfu([210 201],abs(ssx)*[1 1]);
aa6 = annotation('arrow',xx,yy,'Color',p3.Color);
set([aa1 aa2 aa3 aa4 aa5 aa6],'HeadStyle','cback2','HeadWidth',5,'HeadLength',5)


subplot(212)
plot(phi_arr,Sig,'linewidth',1.5)
grid on
hold on
xlim([0 180])
ylim([0.04 0.07])
xlabel('RF spoil phase')
ylabel('Signal / M_0')
title('SPGR steady-state dependence on \Phi_0')
set(gca,'fontsize',13)

%%% Add steady state solutions
[xx,yy]=ds2nfu([0 5],abs(ss0)*[1 1]);
aa1 = annotation('arrow',xx,yy,'Color',p1.Color);
[xx,yy]=ds2nfu([190 181],abs(ss0)*[1 1]);
aa2 = annotation('arrow',xx,yy,'Color',p1.Color);
[xx,yy]=ds2nfu([0 5],abs(ssmt)*[1 1]);
aa3 = annotation('arrow',xx,yy,'Color',p2.Color);
[xx,yy]=ds2nfu([190 181],abs(ssmt)*[1 1]);
aa4 = annotation('arrow',xx,yy,'Color',p2.Color);
[xx,yy]=ds2nfu([0 5],abs(ssx)*[1 1]);
aa5 = annotation('arrow',xx,yy,'Color',p3.Color);
[xx,yy]=ds2nfu([190 181],abs(ssx)*[1 1]);
aa6 = annotation('arrow',xx,yy,'Color',p3.Color);
set([aa1 aa2 aa3 aa4 aa5 aa6],'HeadStyle','cback2','HeadWidth',5,'HeadLength',5)

text(-25,0.08,'(a)','fontsize',20,'fontweight','bold')
text(-25,0.035,'(b)','fontsize',20,'fontweight','bold')

setpospap([360   174   460   460])

print -dpng -r300 bin/Figure2.png

%% EPG-X(BM) with variable delta
phi_arr = 1:1:180;
delta_arr = linspace(0,2,32);%linspace(0,3,32);
nd=length(delta_arr);
nphi=length(phi_arr);
npulse = floor(3*T1/TR);

AA = d2r(alpha)*ones(npulse,1);

if 0 % load stored result if 0
Sig = zeros(nphi,nd);
figure(1)
clf
for ii=1:nphi
   for jj=1:nd
    % Compute RF spoling phase cycles
    phi = RF_phase_cycle(npulse,phi_arr(ii));

    % BM
    tmp = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,kx,'delta',delta_arr(jj)*1e-6*3*42.6e3,'kmax',30);

    Sig(ii,jj) = abs(tmp(end));
   
    
    disp([ii jj])
    
   end
end
    save bin/SigAll_delta Sig
else
    load bin/SigAll_delta
end

%%% Try with isochromat simulation
% 
% %%% Just for phi=0
% phi = RF_phase_cycle(npulse,123);
% 
% Siso = zeros([nd 1]);
% Sepg = zeros([nd 1]);
% for jj=1:nd
%     % EPG-X
%     tic;
%     tmp = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,kx,'delta',delta_arr(jj)*1e-6*3*42.6e3,'kmax',50);
%     toc
%     Sepg(jj)=abs(tmp(end));
%     % Isochromat
%     tic;
%     tmp = isochromat_GRE_BM(AA,phi,TR,T1_MT,T2x,f,k,400,'delta',delta_arr(jj)*1e-6*3*42.6e3,'phase_range',4*pi);
%     toc
%     Siso(jj)=abs(tmp(end));
%     disp(jj)
% end


%% Figure
figure(1)
clf
mm=mesh(delta_arr,phi_arr,Sig);
view([-80 30])
xlabel('\delta_b / ppm')
ylabel('\Phi_0 / deg')
zlabel('Signal/M_0')
set(gca,'fontsize',14)
zlim([0.045 0.075])
setpospap([50 300 1300 440])
colormap parula
print -dpng -r300 bin/SPGR_spoiling_chemshift.png


%% bSSFP with same parameters as above

%%% Do four cases: EPG, EPG-X(MT), EPG-X(BM) with delta=0 and 0.1ppm
delta = 0.1e-6 * 3 * 42.6e3;

% EPG simulations
npulse=5*T1/TR;
phi = RF_phase_cycle(npulse,'balanced');

% single pool
[~,fn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2,'kmax',inf);
mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);

% MT
[~,fnmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,T1_MT,T2_MT,f_MT,k_MT,G,'kmax',inf);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);

% BM
[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,'kmax',inf);
mxyx = ifftshift(fft(ifftshift(fnx,1),[],1),1);
mxyx1 = sum(mxyx,3); % sum compartments

[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,'kmax',inf,'delta',delta);
mxyx = ifftshift(fft(ifftshift(fnx,1),[],1),1);
mxyx2 = sum(mxyx,3); % sum compartments

dw=0;% on resonance 
ss0 = ssSSFP(d2r(alpha),TR,dw,T1,T2);
ssx1 = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,fx,kx,0);
ssx2 = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,fx,kx,delta);
ssmt= ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,dw,T1_MT,T2_MT,f_MT,k_MT,G);


%%% steady state
phiTR = linspace(-pi/TR,pi/TR,npulse);
ssfpx1SS = zeros(size(phiTR));
ssfpx2SS = zeros(size(phiTR));
ssfpmtSS= zeros(size(phiTR));
ssfpss = zeros(size(phiTR));

for ii=1:npulse
    ssfpss(ii) = abs(ssSSFP(d2r(alpha),TR,phiTR(ii),T1,T2));
    ssfpx1SS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,fx,kx,0));
    ssfpx2SS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,fx,kx,delta));
    ssfpmtSS(ii) = abs(ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,phiTR(ii),T1_MT,T2_MT,f_MT,k_MT,G));
end
% analytic form
mxyGloor = ssSSFP_Gloor(d2r(alpha),b1sqrdtau, TR, T1_MT,T2_MT,f_MT,k_MT,G);

%% Figure just with plots

figure(6)
clf

subplot(221)
pp=plot(linspace(-pi,pi,npulse),ssfpss,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxys,1)),squeeze(abs(mxys(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG','location','south')
title('Single compartment')
set(gca,'fontsize',14)

subplot(222)
pp=plot(linspace(-pi,pi,npulse),ssfpmtSS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxymt,1)),squeeze(abs(mxymt(:,end))),'--','linewidth',2)
plot(0,mxyGloor,'*r','MarkerSize',12,'Linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG-X (MT)','Analytic (Gloor 2008)','location','north')
title('White matter MT model')
set(gca,'fontsize',14)


subplot(223)
pp=plot(linspace(-pi,pi,npulse),ssfpx1SS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxyx1,1)),squeeze(abs(mxyx1(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG-X (BM)','location','south')
title('Myelin water exchange model, \delta_b=0')
set(gca,'fontsize',14)

subplot(224)
pp=plot(linspace(-pi,pi,npulse),ssfpx2SS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxyx1,1)),squeeze(abs(mxyx2(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG-X (BM)','location','south')
title('Myelin water exchange model, \delta_b=0.1ppm')
set(gca,'fontsize',14)

setpospap([300 100 700 550])

gg=get(gcf,'children');
ii=[8 6 4 2];
lbls={'(a)','(b)','(c)','(d)'};
for jj=1:4
    axes(gg(ii(jj)))
    text(-4,-0.02,lbls{jj},'fontsize',20,'fontweight','bold')
end
%
print -dpng -r300 bin/Figure3v2.png


%% Supporting figure with approach to steady state
figure(6)
clf

subplot(221)
imagesc(1:npulse,phiTR*TR,abs(mxys),[0 0.2])
title('Single compartment')
ylabel('\psi / rad')
xlabel('TR number')

subplot(222)
imagesc(1:npulse,phiTR*TR,abs(mxymt),[0 0.2])
title('White matter MT model')
xlabel('TR number')
ylabel('\psi / rad')

subplot(223)
imagesc(1:npulse,phiTR*TR,abs(mxyx1),[0 0.2])
title('Myelin water exchange model, \delta_b=0')
xlabel('TR number')
ylabel('\psi / rad')

subplot(224)
imagesc(1:npulse,phiTR*TR,abs(mxyx2),[0 0.2])
title('Myelin water exchange model, \delta_b=0.1ppm')
xlabel('TR number')
ylabel('\psi / rad')

colormap gray

setpospap([300 100 700 550])
print -dpng -r300 bin/SuppFig_ssbSSFP.png

%% Cross validation with isochromat simulations

%%% SPGR (RF spoiling)
% Use simple parameters

%%% Sequences
TR = 5;
alpha = 10;
phi0 = 117;
npulse=200;
phi = RF_phase_cycle(npulse,phi0);

% single pool EPG
[s0,Fn0,Zn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2);
% EPG-X(BM)
[sx,Fnx,Znx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx);

% EPG-X(MT)
[smt,Fnmt,Znmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,T1_MT,T2_MT,f_MT,k_MT,G);


% Isochromat simulations
Niso = [10:20:190 200 300 400 1000];
Siso = {};
rmse = zeros(length(Niso),3);
maxerr=rmse;
% consider error in steady state
erridx = 1:200;

for jj=1:length(Niso)
    
    % Single pool
    siso = isochromat_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2,Niso(jj));
    Siso{jj,1} = siso(:);
    rmse(jj,1) = norm(Siso{jj,1}(erridx)-s0(erridx))/norm(s0(erridx));
    maxerr(jj,1)=max(abs((Siso{jj,1}(erridx)-s0(erridx))./s0(erridx)));
    
    % MT
    siso = isochromat_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1)...
        ,TR,T1_MT,T2_MT,f_MT,k_MT,G,Niso(jj));
    Siso{jj,2} = siso(:);
    rmse(jj,2) = norm(Siso{jj,2}(erridx)-smt(erridx))/norm(smt(erridx));
    maxerr(jj,2)=max(abs((Siso{jj,2}(erridx)-smt(erridx))./smt(erridx)));
    
    % BM
    siso = isochromat_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,Niso(jj));
    Siso{jj,3} = siso(:);
    rmse(jj,3) = norm(Siso{jj,3}(erridx)-sx(erridx))/norm(sx(erridx));
    maxerr(jj,3)=max(abs((Siso{jj,3}(erridx)-sx(erridx))./sx(erridx)));
    
   
end

%% plot these
figfp(1)
nr=3;nc=2;
cm=colormap(lines);
sig_epg = {s0,smt,sx};
titls = {'Single Compartment','White Matter MT model','Myelin water exchange model'};
name2 = {'EPG','EPG-X(MT)','EPG-X(BM)'};
for ii=1:3
    subplot(nr,nc,ii*2-1)
    
    p1=plot(abs(sig_epg{ii}),'k','linewidth',2);
    hold
    pp=[];
    for jj=1:length(Niso)
        pp(jj)=plot(abs(Siso{jj,ii}),'^-','color',cm(jj,:),'linewidth',0.2);
    end
    set(pp,'markersize',2)
    legend(name2{ii},'isochromat predictions')
    uistack(p1,'top');
    xlim([1 200])
    
    grid on
    xlabel('RF pulse number')
    ylabel('Predicted signal')
    title(titls{ii})
    
    subplot(nr,nc,ii*2)
    semilogy(Niso,rmse(:,ii))
    grid on
    xlabel('Number of isochromats (N_{iso})')
    ylabel('RMS difference')
    title(sprintf('%s',titls{ii}))
    ylim([10^-16 1])
    set(gca,'Ytick',10.^(-16:4:0))
    legend('RMS difference')

end
setpospap([100 100 600 650])
gg=get(gcf,'Children');
axes(gg(12))
text(-35,-0.02,'(a)','fontsize',18,'fontweight','bold')
text(-35,-0.3,'(c)','fontsize',18,'fontweight','bold')
text(-35,-0.58,'(e)','fontsize',18,'fontweight','bold')
text(215,-0.02,'(b)','fontsize',18,'fontweight','bold')
text(215,-0.3,'(d)','fontsize',18,'fontweight','bold')
text(215,-0.58,'(f)','fontsize',18,'fontweight','bold')
%
print -dpng -r300 bin/SF1.png