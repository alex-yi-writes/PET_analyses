%% modelling

clear; clc; close all
warning('off','all');

% set paths
paths = [];
paths.parent  = '/Volumes/korokdorf/MRPET/new_recon_registered/';
paths.TACs    = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TACs/NewRecon/';
paths.TACs_new= '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TACs/NewRecon/';
paths.ROImask = '/Volumes/korokdorf/MRPET/coreg_roi/';
paths.figure  = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TACs/NewRecon/tacfig/';
paths.rawImg  = '/Volumes/korokdorf/MRPET/new_recon_registered/';
paths.exports = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/';

addpath(genpath('/Users/alex/Dropbox/paperwriting/MRPET/scripts/modelling'))

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

% etc
settledownFrames=0;

%% load TACs for modelling

for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            TACs{i1,d} = [];
        else
            if (i1==8 & d==1) | (i1==21 && d==1)
                TACs{i1,d} = load([paths.TACs_new num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr_noPVC_troubleshoot_newCerC.mat']);
            else
                TACs{i1,d} = load([paths.TACs_new num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr_noPVC_newCerC.mat']);
            end
        end
    end
end

%% run the model: condensed ROIs
%  because 55 ROIs are too many

clc;
set(0, 'DefaultLineLineWidth', 1);
set(gca,'FontSize',10)
Occ_dyn = cell(length(IDs), 2);

for id = 1:length(IDs)
    for d = 1:2
        if days(id,d) == 0      

            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))

        else

            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))

            TACDATA=[];
            for reg={'CerC','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA',...
                    'caudalanteriorcingulate','entorhinal','inferiortemporal','parahippocampal','precuneus',...
                    'rostralanteriorcingulate','superiortemporal','temporalpole','middletemporal','insula'}
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
            end

            Lengths=[60*ones(60,1)]; % time windows                        % build vector of 60 one-minute frame lengths for inflow
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];           % produce [start,end] time pairs (sec) for each inflow frame
            Lengths=60*ones(15,1);                                         % same idea for 15 baseline frames
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];           % start/end pairs for baseline
            Lengths=60*ones(55,1);                                         % and 55 task frames
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];           % start/end pairs for task frames
            %times=[tt1(10:end,:); tt2+95*60; tt3+115*60];
%             if IDs(id)==4032 && d==2
%             times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3(2:11,:)+115*60];
%             else
            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];      % final assembled frame table with correct absolute offsets
            %             end

            % timevecs
            tmid     = mean(times,2);              % seconds               % midtime (sec) for every frame
            tmidMin  = tmid/60;                    % minutes               % midtime in min because our kinetic formulas use minute (gotta be consistent)
            t_points = numel(tmid);                                        % total number of frames
            dtMin    = [tmidMin(1); diff(tmidMin)];  % *** minutes!!!!***  % frame duration (min) computed from midtimes
            break_point = find(times(:,1)>=115*60,1,'first'); %% time of activation start (>=115 min which is task onset for us) %%
            % break_point = break_point+settledownFrames;

            %%%%%%%%%%%%%%%
            PlotStrFit = 0; % plotting can be really annoying if you're trying to do something in parallel
            %%%%%%%%%%%%%%%

            Tthr=2;                                                        % threshold for zscore outlier check? (we don't use it anywhere yet)
            badcases={}; badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1

                for reg={'Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA'} % iterate through 10 ROIs (actually 9 because striatum is caud+put) we analyse now
                    if PlotStrFit
                        figure('Position',[100 100 800 1200],'Visible','on'); hold on; spidx=1;
                    end
                    reftac=TACDATA.CerC;                                   % ref (cerebellar cortex) TAC
                    reftac=reftac(:);

                    mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2]; % mid-frame refs for trapezoidal integration

                    %% construct linear SRTM design matrix
                    %  later lscov(A, ROItac) solves Aθ=Cᵗ(t) to give θ=[R1 k2 k2a]

                    ASRTM = zeros(t_points ,3);                            % preallocate A matrix (rows = frames, cols = regressors)
                    ASRTM(:,1)  = reftac;                                  % reg 1: instantaneous Cʳ(t): radioactivity concentration in ref tissue (CerC) at time t
                    for k = 1:t_points
                        ASRTM(k,2)  = sum(mreftac(1:k).*dtMin(1:k));       % reg 2: ∫₀ᵗ Cʳ(τ)dτ up to frame k: cumulative area under ref curve up to time t
                    end
                    refauc=sum(mreftac.*dtMin);                            % whole-scan AUC of Cʳ(t) (for QC only i think)

                    if PlotStrFit
                        subplot(3,2,[1 2]); plot(tmidMin,reftac,'k-'); hold on; legendtext={'Cerebellum'};
                    end

                    %% select the target TAC depending on ROI
                    %  will add more if i want to

                    switch(reg{1})
                        case 'Striatum'
                            tempTac=[]; tempVol=[];
                            for subReg={'Caud','Put'}
                                tempTac=[tempTac, TACDATA.(subReg{1})];    % collect caud & put TACs
                                tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol]; % corresponding ROI vols
                            end
                            ROItac=sum(tempTac.*tempVol,2)./sum(tempVol);  savestr='Striatum'; % vol-weighted average
                        case 'Caud', ROItac=TACDATA.Caud; savestr='Caud';
                        case 'Put',  ROItac=TACDATA.Put;  savestr='Put';
                        case 'Nac',  ROItac=TACDATA.Nac;  savestr='Nac';
                        case 'ThalProp', ROItac=TACDATA.ThalProp; savestr='ThalProp';
                        case 'Hipp', ROItac=TACDATA.Hipp; savestr='Hipp';
                        case 'Amy',  ROItac=TACDATA.Amy;  savestr='Amy';
                        case 'SN',   ROItac=TACDATA.SN;   savestr='SN';
                        case 'LC',   ROItac=TACDATA.LC;   savestr='LC';
                        case 'VTA',  ROItac=TACDATA.VTA;  savestr='VTA';
                    end
                    ROItac=ROItac(:);

                    snr_early = mean(ROItac(1:20)) / std(ROItac(1:20));
                    snr_late  = mean(ROItac(end-20:end)) / std(ROItac(end-20:end));
                    fprintf('ROI %s: SNR early=%.2f, late=%.2f\n',savestr,snr_early,snr_late);

                    mROItac  = [ROItac(1)/2; (ROItac(2:end)+ROItac(1:end-1))/2]; % mid-frame target ROI vals
                    taskFrames = length(TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac);
                    baselinetaskFrames = length(TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac)+length(TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac);

                    ASRTM(:,3)=zeros(t_points,1);                          % init 3rd regressor
                    for k = 1:t_points
                        % reg 3: -∫₀ᵗ Cᵗ(τ)dτ
                        ASRTM(k,3)  = -sum(mROItac(1:k).*dtMin(1:k));
                        % Cᵗ(t) is the TAC of the target ROI
                        % this is how SRTM is linearised (eq. 4 in gunn et
                        % al. 1997): with it k2a enters with a negative
                        % weight
                    end

                    %% linear estimate (that's MRTM)

                    [parest,se_srtm,mse_srtm] = lscov(ASRTM,ROItac);       % solve for [R1,k2,k2a]
                    fittac = ASRTM*parest;                                 % reconstructed TAC from linear coefficients
                    BP = parest(2)/parest(3) - 1;                          % BPND MRTM
                    k2p = parest(2)/parest(1);                             % efflux of the ref tissue: parest(1)=R1 (delivery ratio) & parest(2)=k2 (efflux constant for target)
                    % we gotta know k2p whenever we later fix it in the
                    % ESRTM or in lp-ntPET so that fewer parameters vary
                    % during optimisation

                    %% non-linear SRTM fit (whole scan)
                    options = optimset('MaxFunEvals',3e5,'MaxIter',3e5);% config optimisation budget for fminsearch

                    %%%%%%%%%%%

                    % APPROACH 1
                    % trying poisson-variance weights because the occupancy
                    % is still quite noisy in both ESRTM and lp-ntPET.
                    % (Petibon et al 2020, NI; Normandin, Koeppe, & Morris 2012 PMB)
                    % counts = max( ROItac .* dtMin , 1 );                   % expected counts  (activity * framelength in min)
                    % weighs = sqrt(1 ./ counts);                            % inverse-variance weights  (var≈counts)
                    
                    % APPROACH 2 (this one works)

                    if IDs(id)==4011 && d==2
                    weighs = [0*ones(20,1); 0.25*ones(5,1); 0.5*ones(5,1); ones(t_points-30,1)];
                    else
                    weighs = [0.25*ones(10,1); ones(t_points-10,1)];
                    end
                    % weighs = [0.25*ones(30,1); ones(t_points-30,1)];%[0.25*ones(30,1); ones(t_points-30,1)];       % weights (down-weight early inflow frames that are noisy)
                    
                    % APPROACH 3 
                    % find the frame where reference TAC slope first falls below 1 % change/frame
                    % if IDs(id)==4011 && d==2
                    % weighs = [0.25*ones(32,1); ones(t_points-32,1)];
                    % else
                    % refSlope = diff(reftac) ./ reftac(1:end-1);
                    % inflowEnd = find(abs(refSlope) < 0.01, 1, 'first');
                    % if isempty(inflowEnd), inflowEnd = 30; end  % fallback if no clear plateau
                    % weighs = [0.25*ones(inflowEnd,1); ones(t_points-inflowEnd,1)];
                    % end

                    %%%%%%%%%%%%%

                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(...
                        tmidMin,reftac,t_points,x(1),x(2)/x(1),x(3)*ones(t_points,1))-ROItac).*weighs);
                        % tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1))-ROItac).*weighs);
                    % objective: RMS of weighted residuals given R1, k2, BP

                    [parest_srtm,~] = fminsearch(@(x)fobj(x),[1 .3 2],options); % nelder-mead search is used here (searches for the parameter triplet that minimises the weighted root-mean-square residuals)
                    R1__=parest_srtm(1); k2__=parest_srtm(2); BP__=parest_srtm(3);

                    
                    modfit_esrtm = simESRTMfixk2p_1_0_0(...
                        tmidMin,reftac,t_points,R1__,k2__,BP__*ones(t_points,1)); % nonlinear fit TAC

                    %% nonlinear SRTM fit (baseline only)
                    fobj = @(x) norm((simESRTMfixk2p_1_0_0(...
                        tmidMin(1:end-taskFrames),reftac(1:end-taskFrames),t_points-taskFrames,x(1),x(2)/x(1),x(3)*ones(t_points-taskFrames,1)) ...
                        -ROItac(1:end-taskFrames)).*weighs(1:end-taskFrames));

                    [parest_srtm,~] = fminsearch(@(x)fobj(x),[R1__,k2__,BP__],options); % maybe this is better informed?? 
                    % [parest_srtm,~] = fminsearch(@(x)fobj(x),[R1__,k2p,BP__],options); 
                    % [parest_srtm,~] = fminsearch(@(x)fobj(x),[1 .3 2],options);
                    R1_bl=parest_srtm(1); k2_bl=parest_srtm(2); BP_bl=parest_srtm(3);

                    modfit_esrtm_bl = simESRTMfixk2p_1_0_0(...
                        tmidMin,reftac,t_points,R1_bl,k2_bl,BP_bl*ones(t_points,1));

                    %% ESRTM step-change search (noisy)
                    BP_base = BP_bl;                                       % BP baseline
                    k2p_esrtm = k2_bl/R1_bl;                               % fix k2p (this is done to reduce free parameters from 4 to 3, giving the model enough stability for frame-by-frame fitting)

                    fESRTM = @(bp_task) norm( ROItac - simESRTMfixk2p_1_0_0(...
                        tmidMin,reftac,t_points,R1_bl,k2p_esrtm,...
                        [BP_base*ones(break_point-1,1); bp_task*ones(t_points-break_point+1,1)] ) );

                    BP_task   = fminbnd(fESRTM,0,BP_base);                 % by constraining the interval to 0 we forbid any neg occupancy solution and avoid weird optimisation and strange values
                    DeltaBP_ESRTM = (BP_base - BP_task)/BP_base;                     % ΔBP if positive implies DA release
                    Occ_ESRTM   = 100*DeltaBP_ESRTM;               % occupancy in percent

                    %%  lp-ntPET (noisy)
                    Alpntpet=zeros(t_points,4); Alpntpet(:,1:3)=ASRTM;
                    best_mse=10^20; best_parest=[];                        % init trackers

                    for point_rise=break_point:length(tmid)-1              % point_rise is the first frame allowed to rise
                        t2p_index = find(tmid>tmid(point_rise));
                        if length(t2p_index)>1
                            for t_ind=t2p_index(1):t2p_index(end-1)        % try multiple peaks after rise
                                for alpha=[0.25 1 4]                       % try three shape factors (smaller the value spikier it is)
                                    t_peak  = tmid(t_ind)-tmid(point_rise);% when the peak occurs after rise
                                    p       = [1 alpha tmid(point_rise)+t_peak tmid(point_rise)];
                                    actfun  = zeros(size(tmid));
                                    actfun(point_rise:t_points)=gamma_variate_madsen(p,tmid(point_rise:t_points)); % the classical gamma variate fx (A(t))
                                    actfun  = actfun./max(actfun);         % normalise so the peak is 1... now gamma means the peak increase in k2a (al-bachari 2023)

                                    roitac_gamma = ROItac.*actfun;         % product for integral (Cᵗ(t)·A(t)) (normandin et al. 2012)
                                    mroitac_gamma=[roitac_gamma(1)/2; (roitac_gamma(2:end)+roitac_gamma(1:end-1))/2];

                                    Alpntpet(:,4)=0;                       % prep dynamic 4th column
                                    for k=break_point:t_points
                                        Alpntpet(k,4)=-sum(mroitac_gamma(break_point:k).*dtMin(break_point:k));
                                    end

                                    k2a = k2__/(1+BP_bl);                  % fix baseline k2a estimate
                                    
                                    [parest,modelfit]=lpntpetBound(Alpntpet,ROItac,k2a); % constrained LSQ on γ≥0.
                                    mse = sum((ROItac-modelfit).^2);       % residual sume of squares
                                    

                                    if best_mse > mse                      % keep the winner
                                        best_mse=mse;
                                        best_parest=parest;                % [R1 k2 k2a gamma]
                                        best_modelfit=modelfit;
                                        best_actfun=actfun;
                                    end
                                end
                            end
                        end
                    end

                    rssLP{id,d}=sum((ROItac-modelfit).^2); dfLP=4;
                    rssSRTM{id,d}=sum((ROItac-fittac).^2); dfSRTM=3;

                    %%% extra F-test: try accepting gamma only if significant %%%
                    % (otherwise they're noise because it's physiologically
                    % very improbable to be negative)
                    if ~modelCompare(rssSRTM{id,d},rssLP{id,d},dfSRTM,dfLP,t_points)   
                        gamma=NaN; k2a=parest(3);
                    else
                        gamma=parest(4);
                    end

                    y=-ASRTM(:,3)./ROItac; x=ASRTM(:,2)./ROItac;           % logan coordinates
                    [pp,~] = polyfit(x(end-taskFrames:end),y(end-taskFrames:end),1);       % task slope -> BP_task_Logan
                    [pp2,~]= polyfit(x(end-baselinetaskFrames:end-taskFrames),y(end-baselinetaskFrames:end-taskFrames),1); % baseline slope -> BP_base_Logan

                    %% b-ntPET (work in progress, experimenting with other approaches because lt-ntPET is too noisy)
                    % doi: 10.3389/fphys.2020.00498
                    % instead of finding one best fit, it uses markov chain
                    % monte carlo to sample from the posterior distribution
                    % of all param + explicit noise variance
                    
                    % % priors and markov-chain monte-carlo settings
                    % %           R1    k2    k2a   gamma     tD(release time)    dt alpha
                    % theta_min = [1,    0,    0,   0,   tmidMin(break_point),  5, 10]; % min of 7 model param
                    % theta_max = [2,  0.5,  0.1,   k2a,         tmidMin(end), 25, 20]; % max of 7 model param
                    % nIter      = 60000; % num of MCMC iterations (hope this is enough)
                    % burnIn     = 5000;  % how many initial draws to discard (so we only use the well mixed sample)
                    % step       = [0.05, 0.02, 0.005, 0.05, 0.5, 1, 1]; % proposal s.d. for each param before changing
                    % a0 = 1e-3; b0 = 1e-3; % weak inverse-gamma hyperpriors (very loose hyperpriors.. they control how uncertain we are about the overall noise lvl)
                    % 
                    % % init the chain
                    % theta  = [R1__, k2__, k2a, 0, tmidMin(break_point), 10, 15]; % we start at reasonable values from SRTM and set gamma to 0
                    % omega2 = var(ROItac); % initial estimate of the variance of measurement noise (how much ROItac fluctuates)
                    % samples = zeros(nIter,7); % will store every sampled theta vector for later
                    % 
                    % % metropolis-within-gibbs sampling
                    % for j = 1:nIter % start building up a 'cloud' of theta samples
                    % 
                    %     % force gamma to be between 0 and the current k2a
                    %     % (otherwise it doesn't make sense)
                    %     theta_min(4) = 0;
                    %     theata_max(4) = theta(3);
                    % 
                    %     % update each model param theta_k
                    %     for k = 1:7
                    %         prop = theta;
                    %         prop(k) = theta(k) + step(k)*randn; % propose a tiny random change
                    %         if prop(k)>=theta_min(k) && prop(k)<=theta_max(k) % reject if it's outside the minmax
                    %             lp0     = logPosterior(theta,omega2,tmidMin,reftac,ROItac); % otherwise calc log-probability of the data given the old and proposed theta
                    %             lp_prop = logPosterior(prop, omega2,tmidMin,reftac,ROItac);
                    %             if log(rand) < lp_prop - lp0
                    %                 theta = prop;
                    %             end
                    %         end
                    %     end
                    %     % sample omega sqared from its inverse-gamma
                    %     % conditional distribution
                    %     fit    = simulateBnTPET(theta,tmidMin,reftac,ROItac);
                    %     resid  = ROItac - fit;
                    %     a_post = a0 + numel(ROItac)/2;
                    %     b_post = b0 + sum(resid.^2)/2;
                    %     omega2 = 1./gamrnd(a_post,1/b_post);
                    %     samples(j,:) = theta; % store the new theta
                    % end
                    % 
                    % % posterior inference (HPD+MMSE)
                    % post      = samples(burnIn+1:end, :); % discard the first 5k draws and keep only the stable part of the chain
                    % theta_est = hpdm(post, 0.1);    % returns 7*1 vector of MMSE within 10% HPD
                    % [R1_b, k2_b, k2a_b, gamma_b, tD_b, dt_b, alpha_b] = deal(theta_est{:});
                    % 
                    % % reconstruct the fit and calc metrics
                    % [actfun_b, fit_b] = simulateBnTPET(cell2mat(theta_est), tmidMin, reftac, ROItac);
                    % BP_bnt      = k2_b/k2a_b - 1;
                    % BP_peak_b   = k2_b/(k2a_b + gamma_b) - 1;
                    % DeltaBP_bnt = BP_bl - BP_peak_b;
                    % Occ_bnt     = 100 * DeltaBP_bnt / BP_bl;                    
                   
                    %% compute final metrics

                    k2       = best_parest(2); k2a=best_parest(3);
                    BP_lp    = k2/k2a - 1;
                    BP_srtm  = BP__;
                    BP_srtm_bsl=BP_bl;
                    G        = best_parest(4);

                    %%%%%%%%%%%%%%% added on 13-02-2026 %%%%%%%%%%%%%%%%%%%
                    % F-test: lp-ntPET (4 params) vs baseline SRTM (3 params)
                    rss_bl   = sum((ROItac - modfit_esrtm_bl).^2);
                    rss_lp   = sum((ROItac - best_modelfit).^2);
                    F_val    = ((rss_bl - rss_lp) / (4 - 3)) / (rss_lp / (t_points - 4));
                    pF       = 1 - fcdf(F_val, 4-3, t_points-4);
                    if pF >= 0.05
                        G = 0;
                    end
                    % based on what's described in Kim et al. (2014) HBM & 
                    % Wang et al. (2017) NI
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    BPND     = (R1__*k2p)/k2a - 1;
                    DBP      = k2./(k2a+G*best_actfun) - 1;
                    OCC      = 100*(1-DBP/BP_lp);

                    if any(abs([BP_lp  Occ_ESRTM  BP_task]) > 30)
                        fprintf('%s  BP0 %.2f  BP_lp %.2f  Occ_es %.1f  γ %.2e\n',...
                            savestr, BP_bl, BP_lp, Occ_ESRTM, best_parest(4));
                    end

                    %%%%%%%%%%%% 2nd attempt (2025) %%%%%%%%%%%%
                    % lp-ntPET fits one constant baseline k2a (efflux from
                    % target to plasma) & a time-varying extra term
                    % gamma*A(t) that starts at 0 and rises only after
                    % the activation point (break_point/task onset)
                    % because gamma*A(t) is forced to be non-negative and
                    % to be zero before the break_point, the time before
                    % activation is pure baseline.
                    % BP0 = k2 / k2a_baseline – 1 uses only that baseline
                    % k2a, so numerically it is a baseline-only BPND
                    BP0     = k2/k2a-1;
                    BP_peak = k2/(k2a+G)-1;
                    DeltaBP = BP0-BP_peak;
                    Occ_peak= 100*DeltaBP/BP0;
                    disp(['Occ_peak=' num2str(min(max(Occ_peak, 0), 30))]);
                    BP_dyn  = k2./(k2a+G*best_actfun)-1;
                    Occ_dyn{id,d}.(reg{1}) = 100*(1-BP_dyn./BP0);
                    AUC_BP_trap = trapz( tmidMin, BP_dyn );

                    % store
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.'     savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.'    savestr '=OCC;' ])
                    eval(['BPND_save{id,d}.'    savestr '=BPND;' ])
                    eval(['BP_srtm_save{id,d}.' savestr '=BP_srtm;' ])
                    eval(['BP_srtm_Bsl_save{id,d}.' savestr '=BP_srtm_bsl;' ])
                    eval(['BestActivationFunction{id,d}.' savestr '=best_actfun;'])
                    eval(['BestModelFit{id,d}.'          savestr '=best_modelfit;'])
                    eval(['gamma_lp_save{id,d}.'         savestr '=best_parest(4);' ])
                    eval(['Residuals{id,d}.'             savestr '.BaselineFit=ROItac-fittac;' ])
                    eval(['Residuals{id,d}.'             savestr '.AllFit_ESRTM=ROItac-modfit_esrtm;' ])
                    eval(['Residuals{id,d}.'             savestr '.lpntPETFit=ROItac-best_modelfit;' ])
                    eval(['CompensatoryFunction{id,d}.'  savestr '=best_parest(4)*best_actfun;' ])
                   
                    BPdataSave{id,d}.(reg{1}).BP_mrtm=BP;
                    BPdataSave{id,d}.(reg{1}).BP_srtm=BP__;
                    BPdataSave{id,d}.(reg{1}).BP_srtm_bl=BP_bl;
                    BPdataSave{id,d}.(reg{1}).BP_lpnt=BP_lp;
                    BPdataSave{id,d}.(reg{1}).BP_logan=pp(1)-1;
                    BPdataSave{id,d}.(reg{1}).BPND=BPND;
                    BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
                    BPdataSave{id,d}.(reg{1}).DBP=DBP;
                    BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);
                    BPdataSave{id,d}.(reg{1}).DeltaBP_gamma=DeltaBP;
                    BPdataSave{id,d}.(reg{1}).Occ_peak=Occ_peak;
                    BPdataSave{id,d}.(reg{1}).dynBP_AUC=AUC_BP_trap;
                    BPdataSave{id,d}.(reg{1}).DeltaBP_ESRTM = DeltaBP_ESRTM;
                    BPdataSave{id,d}.(reg{1}).Occ_ESRTM     = Occ_ESRTM;
                    % BPdataSave{id,d}.(savestr).BP_bnt    = BP_bnt;
                    % BPdataSave{id,d}.(savestr).DeltaBP   = DeltaBP_bnt;
                    % BPdataSave{id,d}.(savestr).Occ_bnt   = Occ_bnt;

                    % BestActFun{id,d}.(savestr)           = actfun_b;
                    % BestFit_bnt{id,d}.(savestr)          = fit_b;


                    BPdata.BP_mrtm(Subj{r})   = BP;
                    BPdata.BP_srtm(Subj{r})   = BP__;
                    BPdata.BP_srtm_bl(Subj{r})= BP_bl;
                    BPdata.BP_lpnt(Subj{r})   = BP_lp;
                    BPdata.BP_logan(Subj{r})  = pp(1)-1;

                    if PlotStrFit
                        legendtext{end+1}=[reg{1} ' raw'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{Baseline}'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{All}'];
                        legendtext{end+1}=[reg{1} ' fit lpnt-PET_{All}'];
                        legendtext{end+1}=['Activation function'];
                        SE=abs(BP)*sqrt((se_srtm(2)/parest(2))^2+(se_srtm(3)/parest(3))^2);
                        subplot(3,2,[1 2]);
                        yyaxis left;
                        plot(tmidMin,ROItac,'ro',tmidMin,modfit_esrtm_bl,'r-',tmidMin,modfit_esrtm,'b-',tmidMin,best_modelfit,'k--'); hold on;
%                         xlim([0 60]);
                        xlabel('Time (min)');
                        ylabel('Radioactivity concentration');
                        legend(legendtext,'Location','best');
                        yyaxis right;
                        plot(tmidMin,best_parest(4)*best_actfun,'c-'); % plot gamma
                        ylabel('Compensatory function','Color','c');
                        ylim([0 4*10^(-4)]);
                        title([Subj{r} ' compartmental fits: BP_{Baseline}=' num2str(BP_bl,'%1.2f') ', BP_{All}=' num2str(BP__,'%1.2f') ', BP_{lp-nt}=' num2str(BP_lp,'%1.2f')]);
                        subplot(3,2,3);
                        plot(tmidMin,ROItac-fittac,'bo',tmidMin,ROItac-modfit_esrtm,'go',tmidMin,ROItac-modfit_esrtm_bl,'co',tmidMin,ROItac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        plot(tmidMin,ROItac-fittac,'bo',tmidMin,ROItac-modfit_esrtm,'bo',tmidMin,ROItac-modfit_esrtm_bl,'ro',tmidMin,ROItac-best_modelfit,'ko',[0 180],[0 0],'k--');
                        [h p]=runstest(ROItac-fittac);
                        [h1 p1]=runstest(ROItac-modfit_esrtm);
                        [h2 p2]=runstest(ROItac-best_modelfit);
                        if p1<0.05
                            badcases{end+1}=Subj{r};
                        end
                        if p2<0.05
                            badcases2{end+1}=Subj{r};
                        end
                        title(['Residuals (runstest p=' num2str(p,'%1.2f') ', p=' num2str(p1,'%1.2f')  ', p=' num2str(p2,'%1.2f') ')']);
                        xlabel('Time (min)');

                        subplot(3,2,4);
                        y=-ASRTM(:,3)./ROItac;
                        x=ASRTM(:,2)./ROItac;
                        [pp ss]=polyfit(x(end-taskFrames:end),y(end-taskFrames:end),1);
                        [pp2 ss2]=polyfit(x(end-baselinetaskFrames:end-taskFrames),y(end-baselinetaskFrames:end-taskFrames),1);
                        plot(x,y,'ko',x(end-baselinetaskFrames:end),polyval(pp,x(end-baselinetaskFrames:end)),'k-',x(end-baselinetaskFrames:end),polyval(pp2,x(end-baselinetaskFrames:end)),'k--');
                        title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        xlabel(['\int REF/ROI']);
                        ylabel('\int ROI/ROI')
                        subplot(3,2,[5 6]);
                        plot(tmidMin,ROItac./reftac,'ko',[0 180],[0 0],'k--');
                        xlim([0 190]); %ylim([0.5 30]);
                        title('Target to reference ratio');
                        xlabel('Time (min)');

                        set(findall(gcf,'-property','FontSize'),'FontSize',12);

                        % print('-dpsc2','-append','-bestfit',fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_TAC_Fit_lpntpet_logan_' date '_s3.ps']));
                        close(gcf)

                    end

                end   % close ROI loop
            end

            close all
            keep IDs days paths id d TACs BPdata BP_lp_save BPdataSave ...
                DBP_save Occupancy BPND_save BP_srtm_save BP_srtm_Bsl_save ...
                BestActivationFunction BestModelFit gamma_lp_save Residuals ...
                CompensatoryFunction PlotStrFit Occ_dyn rssLP rssSRTM

        end    % close session conditional
    end        % close session loop
end            % closse subject loop

disp('done')

save(['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPpackage_condensed_' ...
    date '_newRecon_noPVC.mat'],...
    'BPdataSave','BP_lp_save','DBP_save','Occupancy','BPND_save',...
    'BP_srtm_save','BP_srtm_Bsl_save','BestActivationFunction',...  
    'BestModelFit','gamma_lp_save','Residuals','CompensatoryFunction')

%% try plotting occupancy?
% 
% close all
% for reg={'Striatum','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA'}
%     figure;
%     for id = 1:length(IDs)
%         for d = 1:2
%             if  days(id,d) == 0
%                 disp('skipped')
%             else
% 
%                 plot(Occ_dyn{id,d}.(reg{1})); hold on
% 
%             end
%         end
%     end
%     title(reg{1});ylim([0 100])
% end

for id = 1:length(IDs)
    for d = 1:2
        if  days(id,d) == 0
            inspectOccPeak_caud(id,d)=NaN; 
            inspectOccPeak_hpc(id,d)=NaN;
            % inspectOccPeak_caud_ES(id,d)=NaN; 
            % inspectOccPeak_hpc_ES(id,d)=NaN;
            inspectdeltaBP_caud(id,d)=NaN;
            inspectdeltaBP_hpc(id,d)=NaN;
            % inspectdeltaBP_caud_ES(id,d)=NaN;
            % inspectdeltaBP_hpc_ES(id,d)=NaN;
        else
             inspectOccPeak_caud(id,d)=BPdataSave{id,d}.Caud.Occ_peak;
             inspectOccPeak_hpc(id,d)=BPdataSave{id,d}.Hipp.Occ_peak;
             % inspectOccPeak_caud_ES(id,d)=BPdataSave{id,d}.Caud.Occ_ESRTM;
             % inspectOccPeak_hpc_ES(id,d)=BPdataSave{id,d}.Hipp.Occ_ESRTM;

            inspectdeltaBP_caud(id,d)=BPdataSave{id,d}.Caud.DeltaBP_gamma;
            inspectdeltaBP_hpc(id,d)=BPdataSave{id,d}.Hipp.DeltaBP_gamma;
            % inspectdeltaBP_caud_ES(id,d)=BPdataSave{id,d}.Caud.DeltaBP_ESRTM;
            % inspectdeltaBP_hpc_ES(id,d)=BPdataSave{id,d}.Hipp.DeltaBP_ESRTM;
             

        end
    end
end

%% save as a excel file

load(['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPpackage_condensed_10-Jul-2025_newRecon_noPVC.mat'], ...
     'BPdataSave','IDs','days');       % bring only what we need

ROIs = {'Amy','Caud','Hipp','LC','Nac','Put','SN','ThalProp', ...
        'VTA'};

metricMap = { ...
%   ROIfield        column stub
   'BP_srtm_bl'     ,'BPbsl_SRTM'      ; ...
   'BP_srtm'        ,'BPtask_SRTM'     ; ...
   'BP_lpnt'        ,'BPbsl_lpnt'      ; ...
   'Occ_ESRTM'      ,'Occupancy_ESRTM' ; ...
   'DeltaBP_ESRTM'  ,'DeltaBP_ESRTM'   ; ...
   'gamma'          ,'Gamma'           ; ...
   'dynBP_AUC'      ,'dynBP_AUC'       ; ...
   'Occ_peak'       ,'Occupancy_lpnt'  };

sessionTag = {'highDA','lowDA'};   % day=1 high / day=2 low

colNames = {};
for m = 1:size(metricMap,1)
    stub = metricMap{m,2};
    for r = 1:numel(ROIs)
        for d = 1:2
            colNames{end+1} = sprintf('%s_%s_%s',stub,ROIs{r},sessionTag{d});
        end
    end
end
nCol   = numel(colNames);
nSubj  = numel(IDs);
values = NaN(nSubj,nCol);           

for s = 1:nSubj                                 % loop subjects
    for d = 1:2                                 % loop sessions
        if days(s,d)==0, continue, end          % skip if session missing

        baseOffset = (d-1);   % 0 for high, 1 for low — used in column index

        for r = 1:numel(ROIs)
            ROI = ROIs{r};
            ROIstruct = BPdataSave{s,d}.(ROI);

            for m = 1:size(metricMap,1)
                fld  = metricMap{m,1};
                stub = metricMap{m,2};

                if isfield(ROIstruct,fld)
                    col = (m-1)*numel(ROIs)*2 + (r-1)*2 + baseOffset + 1;
                    values(s,col) = ROIstruct.(fld);
                end
            end
        end
    end
end

T = array2table(values, ...
        'ID',     cellstr(num2str(IDs(:))), ...
        'VariableNames',colNames);

outFile = ['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPtable_condensedROIs_noPVC_' date '.xlsx'];
writetable(T,outFile,'WriteRowNames',true);
fprintf('excel sheet written:  %s\n',outFile);


%% save: baseline only

% load(['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPpackage_condensed_16-May-2025_newRecon_noPVC.mat'], ...
%      'BPdataSave','IDs','days');       % bring only what we need

ROIs = {'Amy','Caud','Hipp','LC','Nac','Put','SN','ThalProp', ...
        'VTA'};

metricMap = { ...
%   ROIfield        column stub
   'BP_srtm_bl'     ,'BPbsl_SRTM'      ; ...
   % 'BP_srtm'        ,'BPtask_SRTM'     ; ...
   'BP_lpnt'        ,'BPbsl_lpnt'      }; ...
   % 'gamma'          ,'Gamma'           ; ...
   % 'dynBP_AUC'      ,'dynBP_AUC'       ; ...
   % 'Occ_peak'       ,'Occupancy_lpnt'  };

sessionTag = {'highDA','lowDA'};   % day=1 high / day=2 low

colNames = {};
for m = 1:size(metricMap,1)
    stub = metricMap{m,2};
    for r = 1:numel(ROIs)
        for d = 1:2
            colNames{end+1} = sprintf('%s_%s_%s',stub,ROIs{r},sessionTag{d});
        end
    end
end
nCol   = numel(colNames);
nSubj  = numel(IDs);
values = NaN(nSubj,nCol);           

for s = 1:nSubj                                 % loop subjects
    for d = 1:2                                 % loop sessions
        if days(s,d)==0, continue, end          % skip if session missing

        baseOffset = (d-1);   % 0 for high, 1 for low — used in column index

        for r = 1:numel(ROIs)
            ROI = ROIs{r};
            ROIstruct = BPdataSave{s,d}.(ROI);

            for m = 1:size(metricMap,1)
                fld  = metricMap{m,1};
                stub = metricMap{m,2};

                if isfield(ROIstruct,fld)
                    col = (m-1)*numel(ROIs)*2 + (r-1)*2 + baseOffset + 1;
                    values(s,col) = ROIstruct.(fld);
                end
            end
        end
    end
end

T = array2table(values, ...
        'ID',     cellstr(num2str(IDs(:))), ...
        'VariableNames',colNames);

outFile = ['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPbsl_condensedROIs_noPVC_' date '.xlsx'];
writetable(T,outFile,'WriteRowNames',true);
fprintf('excel sheet written:  %s\n',outFile);

%% save: ∆BP_ESRTM only

% load(['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPpackage_condensed_16-May-2025_newRecon_noPVC.mat'], ...
%      'BPdataSave','IDs','days');       % bring only what we need

ROIs = {'Amy','Caud','Hipp','LC','Nac','Put','SN','ThalProp', ...
        'VTA'};

metricMap = { ...
%   ROIfield        column stub
   % 'BP_srtm_bl'     ,'BPbsl_SRTM'      ; ...
   % 'BP_srtm'        ,'BPtask_SRTM'     ; ...
   % 'BP_lpnt'        ,'BPbsl_lpnt'      }; ...
   % 'Occ_peak'       ,'Occ_lpnt'  ; ...  
   % 'gamma'          ,'Gamma'           ; ...
   % 'dynBP_AUC'      ,'dynBP_AUC'       ; ...
   'DeltaBP_ESRTM'  ,'BPchange_ESRTM'    };

sessionTag = {'highDA','lowDA'};   % day=1 high / day=2 low

colNames = {};
for m = 1:size(metricMap,1)
    stub = metricMap{m,2};
    for r = 1:numel(ROIs)
        for d = 1:2
            colNames{end+1} = sprintf('%s_%s_%s',stub,ROIs{r},sessionTag{d});
        end
    end
end
nCol   = numel(colNames);
nSubj  = numel(IDs);
values = NaN(nSubj,nCol);           

for s = 1:nSubj                                 % loop subjects
    for d = 1:2                                 % loop sessions
        if days(s,d)==0, continue, end          % skip if session missing

        baseOffset = (d-1);   % 0 for high, 1 for low — used in column index

        for r = 1:numel(ROIs)
            ROI = ROIs{r};
            ROIstruct = s{s,d}.(ROI);

            for m = 1:size(metricMap,1)
                fld  = metricMap{m,1};
                stub = metricMap{m,2};

                if isfield(ROIstruct,fld)
                    col = (m-1)*numel(ROIs)*2 + (r-1)*2 + baseOffset + 1;
                    values(s,col) = ROIstruct.(fld);
                end
            end
        end
    end
end

% T = array2table(values, ...
%         'ID',cellstr(num2str(IDs(:))), ...
%         'VariableNames',colNames);
T = array2table([IDs(:), values], ...
    'VariableNames', ['ID', colNames]);

outFile = ['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPchange_ESRTM_condensedROIs_noPVC_' date '.xlsx'];
writetable(T,outFile,'WriteRowNames',true);
fprintf('excel sheet written:  %s\n',outFile);

%% remove outliers


for foldcomments=1

clear dataTable
clc;

filename='MRPET_BPchange_ESRTM_condensedROIs_noPVC_13-Jun-2025';

dataTable = readtable(['/Users/alex/Dropbox/paperwriting/MRPET/data/' filename '.xlsx' ]);

% Assuming dataTable is your table
[numRows, numCols] = size(dataTable);

% Loop over each variable (column) in the table
for col = 2:numCols
    dataColumn = dataTable{:, col};
    
    % Variable to count the number of outliers
    numOutliers = 0;

    % Check if the column contains numeric data
    if isnumeric(dataColumn)
        
        clear k1 k2 k3
        
        % Calculate Q1 and Q3
        Q1 = quantile(dataColumn, 0.25);
        Q3 = quantile(dataColumn, 0.75);

        % Calculate IQR
        IQR = Q3 - Q1;

        % Calculate lower and upper bounds
        lowerBound = Q1 - 1.5 * IQR;
        upperBound = Q3 + 1.5 * IQR;

        % Count the number of outliers below the lower bound
        numOutliers = numOutliers + sum(dataColumn < lowerBound);
        k1=find(dataColumn < lowerBound);
        
        % Replace outliers below the lower bound with NaN
        dataColumn(dataColumn < lowerBound) = NaN;

        % Count the number of outliers above the upper bound
        numOutliers = numOutliers + sum(dataColumn > upperBound);
        k2=find(dataColumn > upperBound);
        
        % Replace outliers above the upper bound with NaN
        dataColumn(dataColumn > upperBound) = NaN;

        % Replace the original column with the updated column
        dataTable{:, col} = dataColumn;
        
        % Display the number of outliers for this variable
        if numOutliers==0
        else
            k3=[k1;k2];
            disp('***')
            fprintf('Variable %s had %d outliers.\n', dataTable.Properties.VariableNames{col}, numOutliers);
            
%             try
                disp(strcat('subject number: ', cellstr(num2str(k3))))
%             catch
%                 disp(['subject number: ' num2str(k3(1)) ])
%             end
            
        end
    end
end

writetable(dataTable, ['/Users/alex/Dropbox/paperwriting/MRPET/data/' filename '_outlierRemoved.xlsx'], 'Sheet', 'Cleaned Data');

end
% data_original=data;
% data=dataTable;

disp('done')

%% 

filename = 'MRPET_BPchange_ESRTM_condensedROIs_noPVC_12-Jun-2025';
in  = ['/Users/alex/Dropbox/paperwriting/MRPET/data/' filename '.xlsx'];
out = ['/Users/alex/Dropbox/paperwriting/MRPET/data/' filename '_MADclean.xlsx'];

clc

thr = 4;                                 % raise to 4-5 if this is still too strict

T           = readtable(in);
totalFlagged = 0;

for c = 2:width(T)                       % skip ID column
    if isnumeric(T{:,c})
        colmask = isoutlier(T{:,c},'median','ThresholdFactor',thr);
        idx     = find(colmask);
        if ~isempty(idx)
            totalFlagged = totalFlagged + numel(idx);
            fprintf('\n%s – %d outlier(s)\n',T.Properties.VariableNames{c},numel(idx));
            disp(table(idx,T{idx,1},T{idx,c}, ...
                'VariableNames',{'RowIdx','SubjectID','Value'}))
            T{idx,c} = NaN;              % safe assignment
        end
    end
end

fprintf('\nTotal outliers flagged: %d\n',totalFlagged);
writetable(T,out,'Sheet','Cleaned Data');
disp done



%% fx

% ---------------------------------------------------------
% F-test/ΔAIC (keep lp-ntPET only if it beats SRTM)
% ---------------------------------------------------------
function useLP = modelCompare(rssSRTM,rssLP,dfSRTM,dfLP,n)
F  = ((rssSRTM-rssLP)/(dfLP-dfSRTM)) / (rssLP/(n-dfLP));
pF = 1-fcdf(F,dfLP-dfSRTM,n-dfLP);      % F-test
AIC_S = n*log(rssSRTM/n)+2*dfSRTM;
AIC_L = n*log(rssLP  /n)+2*dfLP;
useLP = (pF<0.05) && (AIC_L < AIC_S);
end

% ---------------------------------------------------------
% lp-ntPET with bound:  0 < gamma ≤ gamma_max (=0.3 k2a)
%    (uses lsqlin nearly same to lscov)
% ---------------------------------------------------------
function [pars,fit] = lpntpetBound(A,y,k2a)
% why is it 0.1765 cap here?: from Riccardi et al. (2006) "Greatest 
% displacements were seen in striatal subdivisions - 5.6% in caudate, 11.2% 
% in putamen, 7.2% in ventral striatum, and 6.6% in substantia nigra. 
% Lesser decrements were seen in amygdala - 4.4%, temporal cortex - 3.7%, 
% and thalamus - 2.8%." ([Riccardi et al., 2006, p. 1016]).
% if it is: gamma=k2a*Occ/(1-Occ), then for 15% occupancy it's
% gamma=(0.15/0.85)*k2a=0.1765*k2a. 
gammaMax = 0.1765*k2a;
gammaMax = max(gammaMax,1e-6);
lb = [-Inf;-Inf;-Inf;   1e-6];
ub = [ Inf; Inf; Inf; gammaMax];
opts = optimoptions('lsqlin','Display','off');
pars = lsqlin(A,y,[],[],[],[],lb,ub,[],opts);  % [R1 k2 k2a gamma]
fit  = A*pars;
end

% ---------------------------------------------------------
% ESRTM  (step model: baseline -> task)
% ---------------------------------------------------------
function [BPbase,BPtask] = esrtm(T,CR,CT,breakIdx,R1,k2p)
% baseline fit
A1 = [CR, cumtrapz(T,CR), -cumtrapz(T,CT)];
p1 = lscov(A1(1:breakIdx,:),CT(1:breakIdx));
k2a1 = p1(3);
% task fit with R1 & k2p fixed
BPbase = k2p/k2a1 - 1;
k2a2   = fminsearch(@(k) norm(CT(breakIdx+1:end) - ...
    esrtmForward(T(breakIdx+1:end),CR(breakIdx+1:end),...
    CT(breakIdx+1:end),R1,k2p,k)), k2a1);
BPtask = k2p/k2a2 - 1;
end

function CT = esrtmForward(T,CR,~,R1,k2p,k2a)
A = [CR, cumtrapz(T,CR), -cumtrapz(T,R1*CR)];   % uses fixed R1,k2p
CT = A*[R1; k2p; k2a];                          % forward TAC
end

% ---------------------------------------------------------
% bayesian ntPET
% ---------------------------------------------------------
function lp = logPosterior(theta, omega2, tmidMin, reftac, ROItac)
    try
        % simulate fit (may error or return nans)
        [~, fit] = simulateBnTPET(theta, tmidMin, reftac, ROItac);
        % if fit contains any non-finite, treat as zero posterior
        if any(~isfinite(fit))
            lp = -Inf;
            return;
        end
        % compute gaussian log‐likelihood
        resid = ROItac - fit;
        ll = -0.5 * sum((resid.^2)./omega2 + log(2*pi*omega2));
        lp = ll;
    catch
        % any error (incl. lsqlin interior-point failure) -> reject
        lp = -Inf;
    end
end


function [actfun, fit] = simulateBnTPET(theta, tmidMin, reftac, ROItac)
R1     = theta(1); 
k2     = theta(2); 
k2a    = theta(3);
gamma  = theta(4); 
tD     = theta(5); 
dt     = theta(6);
alpha  = theta(7); 
t_peak = tD + dt;
% build activation function h(t)
actfun = zeros(size(tmidMin));
idx    = find(tmidMin >= tD);
p      = [1 alpha dt tD];
actfun(idx) = gamma_variate_madsen(p, tmidMin(idx));
actfun = actfun / max(abs(actfun));
% build design matrix Alp
N      = numel(tmidMin);
dtMin  = [tmidMin(1); diff(tmidMin)];
Alp    = zeros(N,4);
Alp(:,1:3) = [reftac, cumtrapz(tmidMin,reftac), -cumtrapz(tmidMin,ROItac)];
for k = idx(1):N
    Alp(k,4) = -sum((ROItac(idx(1):k).*actfun(idx(1):k)) .* dtMin(idx(1):k));
end
% fix k2a_baseline = k2/(1+BP_bl)
k2a_fix = k2 / (1 + (k2/k2a - 1));
try
        [~, fit] = lpntpetBound(Alp, ROItac, k2a_fix);
    catch
        fit = nan(size(ROItac));
end
end

function theta_est = hpdm(samples, prop)
%  Compute for each column the smallest interval containing 'prop'*100% of
%  samples, then return the column-wise mean of samples within that interval.
D = size(samples,2);
theta_est = cell(1,D);
for k = 1:D
    s = sort(samples(:,k));
    L = numel(s);
    win = floor(prop*L);
    minWidth = Inf; idx0 = 1;
    for i = 1:(L-win)
        w = s(i+win) - s(i);
        if w < minWidth
            minWidth = w; idx0 = i;
        end
    end
    theta_est{k} = mean(s(idx0:idx0+win));
end
end
