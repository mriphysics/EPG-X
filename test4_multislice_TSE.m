%% Test 4: Comparison of multislice and single slice TSE
% 2017-09-13
% version 2 - multiple slices as per in-vivo experiments

addpath(genpath('EPGX-src'))
addpath(genpath('lib'))

%%% Sequence
%%% Number of refocusing pulses
nrefocus = 25;

%%% have 8 different scenarios
nslice = 1:2:15;
nn  = length(nslice);

%%% Common sequence parameters
ESP = 7.7;
TR = 5000; %<--- use this to work out time gap at the end



%%% Generate refocusing pulse series
a0={};
a0{1} = d2r([90 180*ones(1,nrefocus)]);
a0{2} = d2r([90 160 120*ones(1,nrefocus-1)]);

npulse = nrefocus+1;
b1sqrdtau={};
b1sqrdtau{1} = [32.7 213.1*ones(1,nrefocus)]; % uT^2 ms
b1sqrdtau{2} = [36.7 189.4 106.5*ones(1,nrefocus-1)]; % uT^2 ms

%%% slice shifts are different for each expt
df = [10.9 12.27] * 42.57e3 * 6e-3; % G * gamma * dx

%%% lineshape
T2b=12e-6;
[ff,G] = SuperLorentzian(T2b);

%%% Just store signal at echo time for each number of slices in WM and GM
sig={};
% sig=zeros([8 2]);
JX=1;

%%% Loop over WM/GM
for IX = 1:2
    
    switch IX
        case 1
            %%% White matter model (see Gloor 2008)
            f = 0.1166;  %% F=0.132 = f/(1-f) => f=0.1166
            kf = 4.3e-3;
            kb = kf * (1-f)/f;
            R1f = 1/779; % ms^-1
            R1b = 1/779; %<- Gloor fix as T1f, usual literature is 1s
            R2f = 1/45;
            
            
        case 2
            %%% Grey matter, caudate nucleus (also Gloor, 2008)
            f = 0.0610;
            kf = 2.3e-3;
            kb = kf * (1-f)/f;
            R1b = 1/1087;
            R1f = 1/1087;
            R2f = 1/59;

    end
    
    
    %%% loop over the experiments
    for JX = 1:2
        
        z0 = {};
        mz = {};
        ss = {};
        
        for jj = 1:nn %<-- loop over number of slices (each one was a separate experiment)
            
            %%% Set up sequence. Slice order goes odd then even.
            slice_order = [1:2:nslice(jj) 2:2:nslice(jj)];
            soi = ceil(nslice(jj)/2);%<-- slice of interest is the middle slice
            
            %%% generate the lineshape for each slice (depends on frequency
            %%% offset)
            fs = df(JX) * (-floor(nslice(jj)/2):floor(nslice(jj)/2));
            GG = interp1(ff,G,fs);
            
            Nsl = length(slice_order);
            Ntr = 4;
            slice_order = repmat(slice_order,[1 Ntr]);
            Ntot = Ntr * Nsl;
            
            %%% now work out delay to add after each slice
            Tshot = ESP*(nrefocus+0.5);
            Tdelay = TR/nslice(jj) - Tshot;
            
            % how to evolve Z0 between slices
            L = [[-R1f-kf kb];[kf -R1b-kb]];
            C = [R1f*(1-f) R1b*f]';
            Xi = expm(L*Tdelay);
            I=eye(2);
            Zoff = (Xi - I)*inv(L)*C;
            
            %%% Initialise magnetization
            z0 = [(1-f) f];
            ss{jj} = [];
            
            % loop over the slices
            for ii=1:Ntot
                
                
                if slice_order(ii)==soi
                    % Local slice
                    [s, Fn,Zn] = EPGX_TSE_MT(a0{JX},b1sqrdtau{JX},ESP,[1/R1f 1/R1b],1/R2f,f,kf,GG(soi),'zinit',z0);
                    % save signal only in this case
                    ss{jj} = cat(2,ss{jj},s(:));
                else
                    % Other slice
                    % set flips to zero and use the saturation lineshape for that slice
                    [s, Fn,Zn] = EPGX_TSE_MT(a0{JX}*0,b1sqrdtau{JX},ESP,[1/R1f 1/R1b],1/R2f,f,kf,GG(slice_order(ii)),'zinit',z0);
                end
                % update z0
                %%% First take Z0 at end of TSE shot
                z0 = squeeze(Zn(1,end,:));
                
                %%% Now evolve it by amount due to recovery period
                z0 = Xi*z0 + Zoff;
                
                
            end
            sig{JX}(jj,IX) = abs(ss{jj}(13,end));%<- record signal in echo 13 of the last simulated TR period
            disp([JX jj IX])
        end
        
    end
end
sig{1}(:,3)=1; % CSF is defined as having no MT effect, so signals are same
sig{2}(:,3)=1;

%% Load in images and ROI measurements (made in another script)
load bin/test4_imagedata

%% Assume in-vivo data was generated elsewhere - read in above

leg={'White Matter','Gray Matter (Caudate)','Cerebrospinal Fluid'};

figure(21)
clf
sl = [1 4 8];
nr=3;nc=3;
fs=18;
%%% 1st image, mark ROIs
subplot(nr,nc,1)
h1=imagesc(repmat(imrotate((abs(ims{1,1})-200)/1000,-90),[1 1 3]));% grayscale
hold on
for ii=1:3
    roitmp = imrotate(roi{ii},-90);
    roitmp = repmat(roitmp,[1 1 3]);
    roitmp(:,:,setdiff(1:3,ii))=0;
    h2 = imagesc(roitmp);
    set(h2, 'AlphaData', 0.5*imrotate(roi{ii},-90));
end
axis off
title('Single slice TSE','fontsize',fs)

for ii=2:nc
    subplot(nr,nc,ii)
    imsjm(ims{sl(ii),1},[200 1200],'rot',-90,'gray')
    title(sprintf('Multislice (%d slices)',nslice(sl(ii))),'fontsize',fs)
    axis off
end

for kk=1:2
    for ii=1:3
        subplot(nr,nc,ii + nc*(kk-1)+nc)
        hold on
        grid on
        pp=[];
        for jj=1:nn
            pp(jj)=plot(nslice(jj),sig{kk}(jj,ii)/sig{kk}(1,ii),'^R');
        end
        set(pp,'markersize',8,'markerfacecolor',[1 0 0])
        if exist('xm','var')
            %%% compute scaling 
            sf = mean(sig{kk}(:,ii)/sig{kk}(1,ii))/mean(xm{ii,kk});
            %pp2=errorbar(1:2:15,xm{ii,kk}/xm{ii,kk}(1),xs{ii,kk}/xm{ii,kk}(1));% in-vivo data
            pp2=errorbar(1:2:15,xm{ii,kk}*sf,xs{ii,kk}*sf);% in-vivo data
            set(pp2,'marker','.','markersize',8,'color',[0 0 0])
        end
        ylim([0.55 1.15])
        if kk==1
            title(leg{ii})
        else
            xlabel('Number of slices')
        end
        ylabel('signal (au)')
        set(gca,'fontsize',12,'ytick',0.5:0.1:1.1)
        xlim([0 16])
    end
end

%
gc = get(gcf,'children');

axes(gc(6))
text(-5,0.62,'180° pulses','rotation',90,'fontsize',18,'fontweight','bold')
axes(gc(3))
text(-5,0.62,'120° pulses','rotation',90,'fontsize',18,'fontweight','bold')

%
%%% reposition
ww=850;hh=650;
setpospap([100 100 ww*0.9 hh*0.9])

set(gc(9),'Position',[0.025 0.55 0.3*([1 (ww/hh)])])
set(gc(8),'Position',[0.35 0.55 0.3*([1 (ww/hh)])])
set(gc(7),'Position',[0.675 0.55 0.3*([1 (ww/hh)])])

for ii=1:3
    gc(ii).Position = [0.7-0.3*(ii-1) 0.07 0.23 0.18];
    gc(ii+3).Position = [0.7-0.3*(ii-1) 0.3 0.23 0.18];
end

print -dpng -r300 bin/Test4_fig1.png

