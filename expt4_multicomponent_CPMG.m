%%% Experiment 4: Multiecho CPMG

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%%% Simplified model. For interpretation, compartment A is the larger
%%% (slow relaxing) and B is the small, fast relaxing compartment (myelin
%%% water?)
f = 0.20; %<--- pool b fraction (this is fast relaxing pool) i.e.  myelin water fraction
fka = 1e-2; % This is ka->b (i.e. Ksf in mcdespot; fb*kb = fa*ka, so kb = (1-fb)*ka/fb -- this is 1/tau
T1 = [1000 500];
T2 = [100 20];


%%% Pulse sequence parameters
ESP=5;
a0 = d2r([90 180*ones(1,50)]);
Npulse = length(a0);
Necho = Npulse-1;




%%% Analyse the echoes using NNLS (Whittall KP, MacKay AL. Quantitative
%%% interpretation of NMR relaxation data. J. Magn. Reson. 1989;84:134?152.)
nt2 = 2001;
t2s = linspace(10,120,nt2);
r2s = 1./t2s;
tt = ESP*(1:Necho);
S =exp(-r2s(:)*tt);



%% 1. Look at estimated T2s/fraction as a function of B1 offset and exchange rate

%%% Leave mt struct as above but modify k term
nk=32;ntx=32;

ka = linspace(0,2.5e-3,nk);
tx = linspace(0.75,1.25,ntx);

t2app = zeros(nk,ntx,2);
fapp = zeros(nk,ntx);
for ii=1:nk
    for jj=1:ntx
        %[s, Fn,Zn,F] = FSE_EPGX_sim(a0*tx(jj),'ESP',ESP,'T1',T1,'T2',T2,'mt',mtt);
        s = EPGX_TSE_BM(a0*tx(jj),ESP,T1,T2,f,ka(ii));
        t2sol = lsqnonneg(S',abs(s(:)));
        [~,pks] = findpeaks(t2sol);
        tmp = t2s(pks);
        t2app(ii,jj,:)=tmp(1:2);%<- in case more than 2 peaks
        % fraction
        fapp(ii,jj) = sum(t2sol(pks(1)+(-5:5))) / sum(t2sol(:)); %<-- look 5 points either side
    end
    disp([ii jj])
end


%% Summary figure

%%% example curve
[s,Fn] = EPGX_TSE_BM(a0*1.1,ESP,T1,T2,f,1e-3);

t2sol = lsqnonneg(S',abs(s(:)));
[~,pks] = findpeaks(t2sol);
t2s(pks)
sum(t2sol(pks(1)+(-5:5)))/ sum(t2sol(:)) %<-- look 5 points either side

figfp(6)
subplot(223)
imagesc(tx,ka*1e3,fapp,[0.1 0.2])
hold
[cc,h]=contour(tx,ka*1e3,fapp,0.1:0.02:0.24);
h.LineColor = [1 1 1];
clabel(cc,h,'Color',[1 1 1],'fontsize',13)

title('Estimated fraction, f')
xlabel('B1 scaling factor')
ylabel('k_a, s^{-1}')
set(gca,'FontSize',12)

colorbar



subplot(224)
imagesc(tx,ka*1e3,t2app(:,:,1),[16 35])
hold
[cc,h]=contour(tx,ka*1e3,t2app(:,:,1),16:4:36);
h.LineColor = [1 1 1];
clabel(cc,h,'Color',[1 1 1],'fontsize',13)

title('Estimated T2b')
xlabel('B1 scaling factor')
ylabel('k_a, s^{-1}')
set(gca,'FontSize',12)

colorbar


        
subplot(2,2,2)
plot(t2s,t2sol)
grid on
hold
xlabel('T2 / ms')
ylabel('au')
title('T2 spectrum from NNLS')
xlim([15 105])
set(gca,'FontSize',12)

subplot(2,2,1)
plot((1:50)*ESP,abs(s))
grid on
hold
plot((1:50)*ESP,squeeze(abs(Fn(52,:,:))))
xlabel('Echo time / ms')
ylabel('Signal (F_0)')
title('Echo amplitudes')
xlim([1 50]*ESP)
set(gca,'FontSize',12)


gg=get(gcf,'children');
%
gg(1).Position = [0.075    0.5838    0.4    0.3412];
gg(2).Position = [0.55    0.5838    0.4    0.3412];
gg(3).Position = [0.91    0.1542    0.0120    0.2334];
gg(4).Position = [0.56    0.1100    0.33    0.3412];
gg(5).Position = [0.43    0.1542    0.0120    0.2334];
gg(6).Position = [0.085    0.1100    0.33    0.3412];

%%% labels
axes(gg(1))
text(-30,-0.1,'(a)','fontsize',16,'fontweight','bold')
text(270,-0.1,'(b)','fontsize',16,'fontweight','bold')
text(-30,-1.5,'(c)','fontsize',16,'fontweight','bold')
text(270,-1.5,'(d)','fontsize',16,'fontweight','bold')

%%% add new axes
axes(gg(2))
ax=axes('Position',[0.6478    0.7002    0.1218    0.1370]);
box on
pp = patch(t2s([-5:5]+pks(1)),t2sol([-5:5]+pks(1)),ones([1 11]));
xlim([20 22])
grid on
pp.FaceColor = [0.75 0.75 0.];
pp.EdgeColor = [0 0.5 0.75];
aa=annotation('arrow',[0.6312 0.5845],[0.6590 0.6123]);

%%% legend
axes(gg(1))
legend('Total echo amplitude','Compartment a','Compartment b')

setpospap([360   231   751   467])
print -r300 -dpng figure9.png
