%% Analyze experimental IR-TFE data (MRF style experiment)

addpath(genpath('EPGX-src'))
addpath(genpath('lib'))

%% Load in data
% Contains im1 = calibration experiment, a series of echoes made with
% increasing flip angles with a large gap in between. Fit this to
% M0sin(alpha*tx) where tx = transmit sensitivity, M0 is effective M0
% including T2*.
%
% im2 is the experimental data from the 'MRF' style experiment

load bin/test2_expt_data

%% Calibration measurements

nf=2;
tx ={};     %<-- B1 scaling
m0={};      %<-- effective M0 (includes receiver coil, and also T2* - readouts are matched)
name = {'MnCl2','BSA'};

flips = (1:15)'*120/15;
xdata = im1{1}(1:16:end-16,:);
sig = @(x)(x(1)*abs(sind(flips*x(2))));

for ii=1:nf
    tx{ii} = zeros([256 1]);
    m0{ii} = zeros([256 1]);
    xdata = im1{ii}(1:16:end-16,:);
    for jj=50:200
        
        xl = xdata(:,jj);
        cf = @(x)(norm(abs(xl)-sig(x)).^2);
        x0 = [16e3 1];
        xs = fminsearch(cf,x0);
        
        tx{ii}(jj) = xs(2);
        m0{ii}(jj) = xs(1);
    end
end

figure(1)
for ii=1:nf
    subplot(2,2,ii)
    plot(tx{ii})
    grid on
    title(sprintf('B1 scaling factor, %s',name{ii}))
    ylim([0.8 1.1])
    
    subplot(2,2,ii+2)
    plot(m0{ii})
    grid on
    title(sprintf('Effective M_0, %s',name{ii}))
    
    
end



%% Make forward EPG predictions with measured calibration values

%%% Pulse sequence properties - flip angle train
ss = [0 0 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15];
ss = [ss fliplr(ss)];
s0 = zeros([1 16]);
alphas = [s0 s0 ss ss s0 s0 ss]; 
alphasc=d2r(sin((alphas/15)*pi/2)*40);
npulse = length(alphas);

%%% Sequence TR, spoiling, inversion pulse
TR=12;
phi_0= [150 117]; %<-- two different spoiling phases

prep=struct;
prep.flip=pi;
prep.t_delay=0; 

%%% Gradients were different for each tube expt (different FOVs)
d={};
d{1}=struct;
d{1}.G = [-5.9 4.5 8.1*sqrt(2)]; % mT/m 
d{1}.tau = [1.4 6.6 3]; %ms
d{2}=struct;
d{2}.G = [-5.6 3.5 sqrt(9.6^2+18.9^2)]; % mT/m 
d{2}.tau = [1.4 8.4 1.1]; %ms



%%% Choose which lines to take from data
xlines = 100+(-2:2);
nx=length(xlines);

xd={};
b1sf = [];
sf=[];
for ii=1:2
    sf(ii) = 1/mean(m0{ii}(xlines));
    
    for jj=1:2
        xd{ii,jj} = sf(ii)*sum(im2{ii,jj}(:,xlines),2)/nx;
    end
    
    b1sf(ii) = mean(tx{ii}(xlines));
end

%%% Sample properties, measured elsewhere
T1 = [899 1290];
T2 = [92 92];
adc = [2.35 1.9]*1e-9;

%%% Single pool EPG simulations
S = {};%<-- with diffusion
S0 = {};%<- without diffusion
for ii=1:2
    d{ii}.D = adc(ii);
    for jj=1:2
        tmp = EPG_GRE(b1sf(ii)*alphasc, RF_phase_cycle(npulse,phi_0(jj)),TR, T1(ii), T2(ii),'prep',prep,'diff',d{ii});
        S{ii,jj} = abs(tmp);
        tmp = EPG_GRE(b1sf(ii)*alphasc, RF_phase_cycle(npulse,phi_0(jj)),TR, T1(ii), T2(ii),'prep',prep);
        S0{ii,jj} = abs(tmp);
    end
end

%% Figure for paper, comparing with and without diffusion

cm=colormap(colorcube);
clrs = {cm(44,:),cm(25,:),cm(48,:)};

figure(3)
clf
for ii=1:2%<-- loop over samples
    for jj=1:2 % <-- loop over spoiling phase
        subplot(4,2,jj*2-1 + ii-1)
        hold on
        plot(xd{ii,jj},'k-*','markersize',2,'markerfacecolor',[0 0 0])
        plot(S0{ii,jj},'color',clrs{3})
        grid on
        xlim([1 256])
        legend('Data','EPG prediction','location','north')
        title(sprintf('\\Phi_0=%d     nrmse = %1.1f %%',phi_0(jj),100*norm(S0{ii,jj}-xd{ii,jj})/norm(xd{ii,jj})))
    end
end

%

for ii=1:2
    for jj=1:2
        subplot(4,2,4+jj*2-1 + ii-1)
        hold on
        plot(xd{ii,jj},'k-*','markersize',2,'markerfacecolor',[0 0 0])
        plot(S{ii,jj},'color',clrs{3})
        grid on
        xlim([1 256])
        legend('Data','EPG prediction','location','north')
        title(sprintf('\\Phi_0=%d     nrmse = %1.1f %%',phi_0(jj),100*norm(S{ii,jj}-xd{ii,jj})/norm(xd{ii,jj})))
    end
end

gg=get(gcf,'Children');
axes(gg(16));

%%% text labels
text(75,0.27,'0.1 mM MnCl_2','fontsize',16,'fontweight','bold')
text(430,0.27,'10% BSA','fontsize',16,'fontweight','bold')
text(-50,-0.2,'EPG w/o Diffusion','fontsize',16,'fontweight','bold','rotation',90)
text(-50,-0.75,'EPG with Diffusion','fontsize',16,'fontweight','bold','rotation',90)

%%% Arrows
aa=[];
aa(1) = annotation('arrow',[0.4537 0.4450],[0.8246 0.8000]);
aa(2) = annotation('arrow',[0.3163 0.3063],[0.6062 0.5908]);
aa(3) = annotation('arrow',[0.4440 0.4338],[0.6062 0.5908]);
aa(4) = annotation('arrow',[0.6453 0.6351],[0.4447 0.4293]);
aa(5) = annotation('arrow',[0.6428 0.6326],[0.2170 0.2016]);

set(aa,'LineWidth',1.5,'HeadLength',7,'HeadWidth',8,'HeadStyle','cback2','color',[0 0 0.8])
set(aa(4:5),'color',[0.95 0.75 0.2])

setpospap([100 40 800 650])

print('-dpng','-r300','bin/Test3b_fig1.png')


%% Now fit the BSA data using the two pool model

IX=2; %<-- corresponds to BSA data

%%% MT requires saturation parameters:
b1sqrdtau = b1sf(IX)^2 * 54.3 * alphasc.^2; % ms uT^2

%%% Now inversion pulse
b1sqrdtauinv = b1sf(IX)^2 * 433.3;
prep=struct;
prep.flip=pi;
prep.t_delay=0; %<-- No delay
prep.B1SqrdTau = b1sqrdtauinv;

%%% Starting value of lineshape
G0=14;

%%% Fit parameters
% [k*1000 f T1a T1b Gsf] Gsf is scaling factor applied to G0, T1a,b in
% seconds

%%% FWD model functions
% 150
sigfun1 = @(x)(EPGX_GRE_MT(b1sf(IX)*alphasc, RF_phase_cycle(npulse,phi_0(1)),...
    b1sqrdtau*b1sf(IX)^2,TR, [x(3)*1e3 x(4)*1e3], T2(IX),x(2),x(1)*1e-3,G0*x(5),'prep',prep,'diff',d{IX},'kmax',30));

% 117
sigfun2 = @(x)(EPGX_GRE_MT(b1sf(IX)*alphasc, RF_phase_cycle(npulse,phi_0(2)),...
    b1sqrdtau*b1sf(IX)^2,TR, [x(3)*1e3 x(4)*1e3], T2(IX),x(2),x(1)*1e-3,G0*x(5),'prep',prep,'diff',d{IX},'kmax',30));

sigfun_tot = @(x)([sigfun1(x);sigfun2(x)]);

% scale signal by 1/(1-f)
sigfun_tot_sc = @(x)(sigfun_tot(x)/(1-x(2)));

% join data
xd_tot = [xd{IX,1};xd{IX,2}];

cf_tot = @(x)(norm(abs(sigfun_tot_sc(x))-xd_tot)^2);


x0 = [3  0.1000 3   0.5  2];
lb = [0   0.05  1.0   0.1  0.5 ];
ub = [8   0.25   5.0   0.99  5 ];


% non linear constraint
nonlcon = @(x)(t1obs_constraint(x,T1(IX),25));

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcns', {  @optimplotx @optimplotfunccount @optimplotfval });

if 0
    x=fmincon(cf_tot,x0,[],[],[],[],lb,ub,nonlcon,options);
else
    % previously found best solution 
    x = [6.2    0.1   1.763    0.363    2.1];
end

R1f = 1e-3/x(3);
kf = x(1)*1e-3;
R1b = 1e-3/x(4);
kb = kf * (1-x(2))/x(2);

R1obs = 0.5*(R1f + kf + R1b + kb)-0.5*sqrt((R1f + kf + R1b + kb).^2 ...
    -4*(R1f*R1b + R1f*kb + R1b*kf));



fprintf(1,'T1bound = %1.1f T1free = %1.1f T1obs =%1.1f kf = %1.3e s^-1 f = %1.3f F=%1.3f G(0) = %1.2f us\n',...
    x(4)*1e3,x(3)*1e3,1/R1obs,x(1),x(2),x(2)/(1-x(2)),x(5)*G0)

fprintf(1,'RMSE before = %1.2f percent, after = %1.2f percent \n',...
100*norm([S{IX,1};S{IX,2}]-xd_tot)/norm(xd_tot),100*norm(abs(sigfun_tot_sc(x))-xd_tot)/norm(xd_tot))



%% Simultaneous fit for single pool version

%%% FWD model functions
% 150
sigfun1s = @(x)(EPG_GRE(b1sf(IX)*alphasc, RF_phase_cycle(npulse,phi_0(1)),...
    TR,x*1e3, T2(IX),'prep',prep,'diff',d{IX},'kmax',30));

% 117
sigfun2s = @(x)(EPG_GRE(b1sf(IX)*alphasc, RF_phase_cycle(npulse,phi_0(2)),...
    TR,x*1e3, T2(IX),'prep',prep,'diff',d{IX},'kmax',30));

sigfun_tots = @(x)([sigfun1s(x);sigfun2s(x)]);

% join data
xd_tot = [xd{IX,1};xd{IX,2}];

cf_tots = @(x)(norm(abs(sigfun_tots(x))-xd_tot)^2);


% only fit T1
lb = 0;
ub = 5;
x0 = 3;

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcns', {  @optimplotx @optimplotfunccount @optimplotfval });
sol=fmincon(cf_tots,x0,[],[],[],[],lb,ub,[],options);

solact_single=sol;

figfp(13)
subplot(2,1,1)
plot(xd_tot)
hold on
plot(abs(sigfun_tots(T1(IX)/1000)))
plot(abs(sigfun_tots(solact_single)))
subplot(2,1,2)
hold on
plot(xd_tot-abs(sigfun_tots(T1(IX)/1000)))
plot(xd_tot-abs(sigfun_tots(solact_single)))

fprintf(1,'T1 fit = %d, error = %1.3f percent \n',fix(solact_single*1e3),...
    100*norm(abs(sigfun_tots(solact_single))-xd_tot)/norm(xd_tot))

%% Final figure showing both sets of residuals

fprintf(1,'T1bound = %1.1f T1free = %1.1f T1obs =%1.1f kf = %1.3e s^-1 f = %1.3f F=%1.3f G(0) = %1.2f us\n',...
    x(4)*1e3,x(3)*1e3,1/R1obs,x(1),x(2),x(2)/(1-x(2)),x(5)*G0)

fprintf(1,'RMSE before = %1.2f percent, after = %1.2f percent \n',...
100*norm([S{IX,1};S{IX,2}]-xd_tot)/norm(xd_tot),100*norm(abs(sigfun_tot_sc(x))-xd_tot)/norm(xd_tot))

cm=colormap(colorcube);
clrs = {cm(48,:),cm(25,:),cm(44,:)};
figfp(11);
lw=1;
fs=15;

subplot(221)
plot(xd{IX,1},'k-*','markersize',3)
hold on
plot(S{IX,1},'linewidth',lw,'color',clrs{1})
plot(abs(sigfun1s(solact_single)),'linewidth',lw,'color',clrs{2})
plot(abs(sigfun1(x))/(1-x(2)),'-','linewidth',lw,'color',clrs{3})

grid on
xlim([1 256])
% ll1=legend('Data','EPG (not fitted)','EPG fitted','EPG-X(MT) fitted','location','north');
ylabel('Signal/M_0')
xlim([0 npulse])
set(gca,'XTick',0:32:256)
xlabel('RF pulse number')
set(gca,'FontSize',fs)
title('Fitted data: \Phi_0=150°')

subplot(222)

plot(xd{IX,2},'k-*','markersize',3,'markerfacecolor',[0 0 0])
hold on
plot(S{IX,2},'linewidth',lw,'color',clrs{1})
plot(abs(sigfun2s(solact_single)),'linewidth',lw,'color',clrs{2})
plot(abs(sigfun2(x))/(1-x(2)),'-','linewidth',lw,'color',clrs{3})

grid on
xlim([1 256])
ll1=legend('Data','EPG (not fitted)','EPG fitted','EPG-X(MT) fitted','location','north');

ylabel('Signal/M_0')
xlim([0 npulse])
set(gca,'XTick',0:32:256)
xlabel('RF pulse number')
set(gca,'FontSize',fs)
title('Fitted data: \Phi_0=117°')


subplot(223)
plot(abs(xd{IX,1}-S{IX,1}),'linewidth',lw,'color',clrs{1})
hold on
plot(abs(xd{IX,1}-abs(sigfun1s(solact_single))),'-','linewidth',lw,'color',clrs{2})
plot(abs(xd{IX,1}-abs(sigfun1(x))/(1-x(2))),'-','linewidth',lw,'color',clrs{3})
grid on
xlim([1 256])
%legend('EPG (not fitted)','EPG fitted','EPG-X(MT) fitted','location','north')
ylabel('error')
xlim([0 npulse])
ylim([0 0.02])
set(gca,'XTick',0:32:256)
xlabel('RF pulse number')
set(gca,'FontSize',fs)
title('Residuals, \Phi_0=150°')


subplot(224)
plot(abs(xd{IX,2}-S{IX,2}),'linewidth',lw,'color',clrs{1})
hold on
plot(abs(xd{IX,2}-abs(sigfun2s(solact_single))),'-','linewidth',lw,'color',clrs{2})
plot(abs(xd{IX,2}-abs(sigfun2(x))/(1-x(2))),'-','linewidth',lw,'color',clrs{3})
grid on
xlim([1 256])
ll2=legend('EPG (not fitted)','EPG fitted','EPG-X(MT) fitted','location','north');
ylabel('error')
xlim([0 npulse])
ylim([0 0.02])
set(gca,'XTick',0:32:256)
xlabel('RF pulse number')
set(gca,'FontSize',fs)
title('Residuals, \Phi_0=117°')

setpospap([100 200 1035 480])

%%% Move
gc=get(gcf,'children');
gc(6).Position = [0.08 0.61 0.33 0.3];
gc(5).Position = [0.48 0.61 0.33 0.3];
gc(4).Position = [0.82 0.7 0.15 0.18];
gc(3).Position = [0.08 0.13 0.33 0.3];
gc(2).Position = [0.48 0.13 0.33 0.3];
gc(1).Position = [0.82 0.2 0.15 0.18];

print('-dpng','-r300','bin/Test3b_fig2.png')
