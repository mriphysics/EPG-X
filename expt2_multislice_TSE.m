%% Experiment 2: Comparison of multislice and single slice TSE

addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%%% White matter model (see Gloor 2008)
f = 0.1166;  %% F=0.132 = f/(1-f) => f=0.1166
kf = 4.3e-3;
kb = kf * (1-f)/f;
R1f = 1/779; % ms^-1
R1b = 1/1000; %<- Gloor fix as T1f, usual literature is 1s
R2f = 1/45;

%%% Run two simulations, each with TSE factor 27
nrefocus = 27;

a0={};

% Case 1 - 180 pulses
aa = FSE_alpha_oneahead_v2(1.0*ones(1,nrefocus));
a0{1} = [pi/2 abs(aa)];
fprintf(1,'RMS flip = %1.2f\n',norm(r2d(aa))/sqrt(length(aa)));

% Case 2 - ~70° pulses 
aa = FSE_alpha_oneahead_v2(0.7*ones(1,nrefocus));
a0{2} = [pi/2 abs(aa)];
fprintf(1,'RMS flip = %1.2f\n',norm(r2d(aa))/sqrt(length(aa)));

%%% Common sequence parameters
ESP = 7.3;
TR = 5000; %<--- use this to work out time gap at the end

slice_order = [1 3 5 2 4];%<- slice of interest is number 3
soi = 3;%<- slice of interest is number 3

Nsl = length(slice_order);
Ntr = 3;
slice_order = repmat(slice_order,[1 Ntr]);
Ntot = Ntr * Nsl;

%%% now bulk out the shot to include dummy pulses
Tshot = TR/Nsl;
Tdelay = Tshot - ESP*(nrefocus+0.5);
ndummy = round(Tdelay/ESP);

for jj=1:2
    if ndummy>0
        a0{jj} = [a0{jj} zeros(1,ndummy)];
    end
end
npulse = length(a0{1}); % same


%%% lineshape & slice offset frequencies
T2b=12e-6;
[ff,G] = SuperLorentzian(T2b);
df = [-2 -1 0 1 2]*1e3;

GG = interp1(ff,G,df); %<-- saturation factor for each slice


z0_ms = {};
mz_ms = {};
ss_ms = {};
z0_ss = {};
mz_ss = {};
ss_ss = {};

for jj=1:2 %<-- loop oevr flip trains
    
    %%% MT pool saturation depends on flips
    b1sqrdtau = 22 * a0{jj}.^2; % ms uT^2
  
    %%% simulate both cases
    z0_ms{jj} = [(1-f) f]; %<-- multi slice
    mz_ms{jj} = zeros([npulse-1 Ntot 2]);
    ss_ms{jj} = zeros([npulse-1 Ntot]);
    z0_ss{jj} = [(1-f) f]; %<-- single slice
    mz_ss{jj} = zeros([npulse-1 Ntot 2]);
    ss_ss{jj} = zeros([npulse-1 Ntot]);
    
    
    
    for ii=1:Ntot
        
        %%% First do multi-slice case
        if slice_order(ii)==soi
            % Local slice
            [s, Fn,Zn] = EPGX_TSE_MT(a0{jj},b1sqrdtau,ESP,[1/R1f 1/R1b],1/R2f,f,kf,GG(1),'zinit',z0_ms{jj});
        else
            % Other slice
            % set flips to zero and use the saturation lineshape for that slice
            [s, Fn,Zn] = EPGX_TSE_MT(a0{jj}*0,b1sqrdtau,ESP,[1/R1f 1/R1b],1/R2f,f,kf,GG(slice_order(ii)),'zinit',z0_ms{jj});
        end
        % update z0
        z0_ms{jj} = squeeze(Zn(1,end,:));
        mz_ms{jj}(:,ii,:) = reshape(Zn(1,:,:),[npulse-1 1 2]);
        ss_ms{jj}(:,ii) = s;
        
        %%% Now single slice case - off resonant slices do nothing, but Mz evolves still
        if slice_order(ii)==soi
           [s, Fn,Zn] = EPGX_TSE_MT(a0{jj},b1sqrdtau,ESP,[1/R1f 1/R1b],1/R2f,f,kf,GG(1),'zinit',z0_ss{jj});
        else
            % set flips to zero but keep saturation
            [s, Fn,Zn] = EPGX_TSE_MT(a0{jj}*0,b1sqrdtau*0,ESP,[1/R1f 1/R1b],1/R2f,f,kf,GG(1),'zinit',z0_ss{jj});
        end
        % update z0
        z0_ss{jj} = squeeze(Zn(1,end,:));
        mz_ss{jj}(:,ii,:) = reshape(Zn(1,:,:),[npulse-1 1 2]);
        ss_ss{jj}(:,ii) = s;
        
    end
end

%% Plots

fs=13;

figfp(2)

for jj=1:2
    subplot(2,2,jj*2-1)
    plot(reshape(mz_ms{jj},[(npulse-1)*Ntot 2]))
    hold on
    plot(reshape(mz_ss{jj},[(npulse-1)*Ntot 2]),'-.')
    if jj==2
        legend({'Multi Slice Z_0^a','Multi Slice {Z}_0^{b}',...
            'Single Slice {Z}_0^{a}','Single Slice {Z}_0^{b}'},...
            'location','Eastoutside','FontSize',11)
    end
    grid on
    xlabel('Excited slice number','fontsize',fs)
    ylabel('$$\tilde{Z}_0^{a,b} / M_0$$','interpreter','latex','fontsize',14)

    switch jj
        case 1
            title('M_z over 3 TR periods, RMS flip = 180')
        case 2
            title('M_z over 3 TR periods, RMS flip = 68')
    end
    set(gca,'xtick',0:(npulse-1):(npulse-1)*Ntot,'xticklabel',slice_order,'fontsize',fs)
    ylim([0 1])
    xlim([0 (npulse-1)*Ntot])
    %%% patch objects
    p1 = [];
    p2=[];
    for ii=1:Ntr
        for kk=1:Nsl
            if slice_order(kk)==soi
                tmp = patch((npulse-1)*Nsl*(ii-1)+(kk-1)*npulse+[0 length(aa) length(aa) 0],...
            [0 0 1 1],[1 1 1]);
                p1 = [p1 tmp];
            else
                tmp = patch((npulse-1)*Nsl*(ii-1)+(kk-1)*npulse+[0 length(aa) length(aa) 0],...
            [0 0 1 1],[1 1 1]);
                p2 = [p2 tmp];
            end
        end
    end
    set(p1,'EdgeAlpha',0,'FaceAlpha',0.05,'FaceColor',[1 0 0])
    set(p2,'EdgeAlpha',0,'FaceAlpha',0.05,'FaceColor',[0 1 0])
end

%%% Echo plots
soi_idx = min(find(slice_order==soi)); %<- this is the order of the slice of interest in time
for jj=1:2
    
    idx = 1:27;
    
    subplot(2,2,2*jj)
    plot(abs(ss_ms{jj}(idx + (npulse-1)*(Nsl*(Ntr-1)+(soi_idx-1)))))
    hold on
    plot(abs(ss_ss{jj}(idx + (npulse-1)*(Nsl*(Ntr-1)+(soi_idx-1)))))
    grid on
    xlim([1 nrefocus])
    switch jj
        case 1
            title('Echo amplitudes, RMS flip = 180','fontsize',fs)
        case 2
            title('Echo amplitudes, RMS flip = 68','fontsize',fs)
    end
    legend('Multi Slice','Single Slice')
    xlabel('Echo number','fontsize',fs)
    ylabel('$$\tilde{F}_0 / M_0$$','interpreter','latex','fontsize',14)
    echo_idx = 11; %<-- TE=80ms
  
    plot(echo_idx,abs(ss_ms{jj}(idx(echo_idx)+ (npulse-1)*(Nsl*(Ntr-1)+(soi_idx-1)))),'*k')
    plot(echo_idx,abs(ss_ss{jj}(idx(echo_idx)+ (npulse-1)*(Nsl*(Ntr-1)+(soi_idx-1)))),'*k')
    text(echo_idx,0.4,sprintf('Ratio = %1.2f',abs(ss_ms{jj}(idx(echo_idx)+ (npulse-1)*(Nsl*(Ntr-1)+(soi_idx-1))))/abs(ss_ss{jj}(idx(echo_idx)+ (npulse-1)*(Nsl*(Ntr-1)+(soi_idx-1))))),'fontweight','bold','fontsize',14)
end

gc = get(gcf,'Children');

%
set(gc(2),'Position',[0.72    0.1100    0.245    0.3412])
set(gc(4),'Position',[0.72    0.5838    0.245    0.3412])
set(gc(5),'Position',[0.5298    0.3961    0.1415    0.1531])
set(gc(6),'Position',[0.05    0.1100    0.608    0.3412])
set(gc(7),'Position',[0.05    0.5838    0.608    0.3412])


setpospap([110         142        1071         540])

print('-dpng','-r300','FSE_multislice_singleslice.png')