%% modelling

clear; clc; close all
warning('off','all');

% set paths
paths = [];
paths.parent  = '/Volumes/korokdorf/MRPET/new_recon_registered/';
paths.TACs    = '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/NewRecon/noPVC/';
paths.TACs_new= '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/NewRecon/noPVC/';
paths.ROImask = '/Volumes/korokdorf/MRPET/coreg_roi/';
paths.figure  = '/Users/alex/Dropbox/paperwriting/MRPET/data/TACs/NewRecon/tacfig/';
paths.rawImg  = '/Volumes/korokdorf/MRPET/new_recon_registered/';
paths.exports = '/Users/alex/Dropbox/paperwriting/MRPET/data/';

addpath(genpath('/Users/alex/Dropbox/paperwriting/MRPET/scripts/modelling'))

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

%% load TACs for modelling

for i1 = 1:length(IDs)
    for d = 1:2
        if days(i1,d) == 0
            TACs{i1,d} = [];
        else
            if (i1==8 & d==1) | (i1==21 && d==1) % people had to come out of the scanner earlier during inflow scan because they wanted to go to the bathroom
            TACs{i1,d} = load([paths.TACs_new num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr_noPVC_troubleshoot.mat']);
            else
            TACs{i1,d} = load([paths.TACs_new num2str(IDs(i1)) num2str(d) '/' num2str(IDs(i1)) num2str(d) '_TACs_deccorr_noPVC.mat']);
            end
        end
    end
end


%% run the model: condensed ROIs

set(0, 'DefaultLineLineWidth', 1);

for id = 1%:length(IDs)
    for d = 1%:2 % which session?: 1 = high reward, 2 = low reward
        if days(id,d) == 0
            
            fprintf(['\n *************\n no session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
        else
            
            fprintf(['\n *************\n analysing \n session %1.d data for ID %4.d\n *************\n'],d,IDs(id))
            
            TACDATA=[];
            
            for reg={'CerC','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA'}
            % for reg={'CerC','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA',...
            %             'caudalanteriorcingulate','entorhinal','inferiortemporal','parahippocampal','precuneus',...
            %             'rostralanteriorcingulate','superiortemporal','temporalpole','middletemporal','insula'}
                TACDATA.(reg{1})=[TACs{id,d}.TACDATA_InFlow.(reg{1}).Bilateral.tac; ...
                    TACs{id,d}.TACDATA_Baseline.(reg{1}).Bilateral.tac;...
                    TACs{id,d}.TACDATA_Task.(reg{1}).Bilateral.tac];
            end
            
            % please adjust time bins according to your PET data - in my
            % case, every frame was 60 seconds
            Lengths=[60*ones(60,1)]; % inflow, 60 frames & 60 seconds each
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];

            Lengths=60*ones(15,1); % baseline, 15 frames & 60 seconds each
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];

            Lengths=60*ones(55,1); % task, 55 frames & 60 seconds each
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)];

            times=[tt1(1:length(TACs{id,d}.TACDATA_InFlow.CerC.Bilateral.tac),:); tt2+95*60; tt3+115*60];
            
            tmid=mean(times,2);
            tmidMin=tmid/60;
            t_points    = length(tmid);
            dt      = [tmid(1); tmid(2:length(tmid))-tmid(1:length(tmid)-1)];
            break_point=find(times(:,1)>=115*60,1,'first'); %% time of activation start
            
            %%%%%%%%%%%
            PlotStrFit=1; % show modelling plots or not (the plots are saved in a folder regardless)
            %%%%%%%%%%%
            
            Tthr=2;
            badcases={};
            badcases2={};
            BPdata=array2table(NaN*zeros(1,5));
            Subj={[num2str(IDs(id)) num2str(d)]};
            BPdata.Properties.RowNames=Subj;
            BPdata.Properties.VariableNames={'BP_mrtm','BP_srtm','BP_srtm_bl','BP_lpnt','BP_logan'};
            for r=1
                for reg={'Striatum','Caud','Put','Nac','ThalProp','Hipp','Amy','SN','LC','VTA','CerC'} % loop through dopaminergic ROIs: please change this according to your ROI list
                    if PlotStrFit
                        figure('Position',[100 100 800 1200],'Visible','on'); hold on;
                        spidx=1;
                    end
                    reftac=TACDATA.CerC;
                    mreftac  = [reftac(1)/2; (reftac(2:end)+reftac(1:end-1))/2];
                    %% set the SRTM part of A

                    ASRTM = zeros(t_points ,3);
                    ASRTM(:,1)  = reftac;
                    for k = 1:t_points
                        ASRTM(k,2)  = sum(mreftac(1:k).*dt(1:k));
                    end
                    refauc=sum(mreftac.*dt);
                    
                    if PlotStrFit
                        subplot(3,2,[1 2]);
                        plot(tmidMin,reftac,'k-'); hold on;
                        legendtext={'Cerebellum'};
                    end
                    switch(reg{1})
                        case 'Striatum'
                            tempTac=[];
                            tempVol=[];
                            for subReg={'Caud','Put'}
                                tempTac=[tempTac, TACDATA.(subReg{1})];
                                tempVol=[tempVol, TACs{id,d}.TACDATA_Baseline.(subReg{1}).Bilateral.vol];
                            end
                            ROItac=sum(tempTac.*tempVol,2)./sum(tempVol);
                            savestr = 'Striatum';
                        case 'Caud'
                            ROItac=TACDATA.Caud;
                            savestr = 'Caud';
                        case 'Put'
                            ROItac=TACDATA.Put;
                            savestr = 'Put';
                        case 'Nac'
                            ROItac=TACDATA.Nac;
                            savestr = 'Nac';
                        case 'ThalProp'
                            ROItac=TACDATA.ThalProp;
                            savestr = 'ThalProp';
                        case 'Hipp'
                            ROItac=TACDATA.Hipp;
                            savestr = 'Hipp';
                        case 'Amy'
                            ROItac=TACDATA.Amy;
                            savestr = 'Amy';
                        case 'Pall'
                            ROItac=TACDATA.Pall;
                        case 'SN'
                            ROItac=TACDATA.SN;
                            savestr = 'SN';
                        case 'LC'
                            ROItac=TACDATA.LC;
                            savestr = 'LC';
                        case 'VTA'
                            ROItac=TACDATA.VTA;
                            savestr = 'VTA';    
                        % case 'caudalanteriorcingulate'
                        %     ROItac=TACDATA.caudalanteriorcingulate;
                        %     savestr = 'caudalanteriorcingulate';
                        % case 'entorhinal'
                        %     ROItac=TACDATA.entorhinal;
                        %     savestr = 'entorhinal';
                        % case 'inferiortemporal'
                        %     ROItac=TACDATA.inferiortemporal;
                        %     savestr = 'inferiortemporal';
                        % case 'parahippocampal'
                        %     ROItac=TACDATA.parahippocampal;
                        %     savestr = 'parahippocampal';
                        % case 'precuneus'
                        %     ROItac=TACDATA.precuneus;
                        %     savestr = 'precuneus';
                        % case 'rostralanteriorcingulate'
                        %     ROItac=TACDATA.rostralanteriorcingulate;
                        %     savestr = 'rostralanteriorcingulate';
                        % case 'superiortemporal'
                        %     ROItac=TACDATA.superiortemporal;
                        %     savestr = 'superiortemporal';
                        % case 'temporalpole'
                        %     ROItac=TACDATA.temporalpole;
                        %     savestr = 'temporalpole';
                        % case 'middletemporal'
                        %     ROItac=TACDATA.middletemporal;
                        %     savestr = 'middletemporal';
                        % case 'insula'
                        %     ROItac=TACDATA.insula;
                        %     savestr = 'insula';
                        case 'CerC'
                            ROItac=TACDATA.CerC;
                            savestr = 'CerC';    
                    end
                    
                    mROItac  = [ROItac(1)/2; (ROItac(2:end)+ROItac(1:end-1))/2];
                    ASRTM(:,3)=zeros(t_points,1);
                    for k = 1:t_points
                        ASRTM(k,3)  = -sum(mROItac(1:k).*dt(1:k));
                    end

                    % ====== multilinear SRTM (MRTM) using lscov ======
                    [parest,se_srtm,mse_srtm] = lscov(ASRTM,ROItac);
                    %  parest(1) => R1=K1/K1'
                    %  parest(2) => V_T
                    %  parest(3) => V_ND
                    fittac = ASRTM * parest;
                    BP = parest(2)/parest(3) - 1;
                    k2p = parest(2)/parest(1);  

                    % %LSQ-estimation using lscov
                    % [parest se_srtm mse_srtm]   = lscov(ASRTM,ROItac); % parest(1)->R1=K1/K1', parest(2)->V_T, parest(3)->V_ND
                    % fittac=ASRTM*parest;
                    % BP=parest(2)/parest(3)-1;
                    % k2p=parest(2)/parest(1); % k2 of the reference region
                    

                    % ============ jarkko's original 'real SRTM' (ESRTM) ============= %
                    options = optimset('MaxFunEvals',6000, 'MaxIter',6000);
                    weighs  = [0.25*ones(30,1); ones(t_points-30,1)];
                    fobj_srtm = @(x) norm(( ...
                        simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points,x(1),x(2),x(3)*ones(t_points,1)) ...
                        - ROItac).*weighs );
                    [parest_srtm, minnorm] = fminsearch(@(x) fobj_srtm(x),[1 .3 2],options);
                    R1__=parest_srtm(1);
                    k2__=parest_srtm(2);
                    BP__=parest_srtm(3);
                    modfit_esrtm=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points, ...
                        R1__, k2__, BP__*ones(t_points,1));
                    
                    % SRTM up to end of baseline
                    fobj_bsl = @(x) norm(( ...
                        simESRTMfixk2p_1_0_0(tmidMin(1:end-55),reftac(1:end-55),t_points-55,x(1),x(2),x(3)*ones(t_points-55,1)) ...
                        - ROItac(1:end-55)).*weighs(1:end-55));
                    [parest_srtm_bsl, minnorm_bsl] = fminsearch(@(x) fobj_bsl(x),[1 .3 2],options);
                    R1_bl=parest_srtm_bsl(1);
                    k2_bl=parest_srtm_bsl(2);
                    BP_bl=parest_srtm_bsl(3);
                    modfit_esrtm_bl=simESRTMfixk2p_1_0_0(tmidMin,reftac,t_points, ...
                        R1_bl, k2_bl, BP_bl*ones(t_points,1));


                    % ============================================================ % 
                    % %%%%% SRTM2 %%%%%
                    %
                    % now fit only 2 parameters (R1, BP) using the new simSRTM2_1_0_0

                    fobj_ref = @(x) norm( sim1TC_ref(tmidMin, t_points, x(1), x(2)) - reftac ); % x(1)->K1, x(2)->k2
                    [xbest, ~] = fminsearch(fobj_ref, [1, 0.1]);

                    K1_ref  = xbest(1);
                    k2_ref  = xbest(2); % best estimate of the efflux from the reference region (our k2p)
                    fprintf('fitted single-tissue ref: K1=%.3f, k2=%.3f\n',K1_ref,k2_ref); % check if they're reasonable

                    % let's experiment with weighting:

                    % start with everything = 1:
                    weighs2 = ones(t_points, 1); % initialise

                    % (1) inflow frames 1..60 -> down-weight them, for example 0.25
                    weighs2(1:9) = 0; % exclude the first few noisy frames - the earliest time points might be affected by vascular or bolus effects
                    weighs2(10:60) = 1;

                    % (2) Baseline frames 61..75 -> normal weight (already 1), so do nothing
                    % weighs2(61:75) = 1.0;

                    % (3) task frames 76..130 ->
                    %       give normal weight (1) to the first 76..112
                    %       extra weight (2) to the last 18 frames 113..130
                    % “Specifically bound [18F]fallypride was rapidly displaced with haloperidol 
                    % (1 mg/kg, intravenous) with a dissociation rate of 0.0385 min^-1 and thus 
                    % a halftime, t_1/2 = 18 min." (Mukherjee et al., 1995, p. 13)
                    weighs2(113:130) = 2.0;

                    fobj_srtm2 = @(p) norm(( ...
                       simSRTM2_1_0_0(tmidMin, reftac, t_points, k2_ref, p(1), p(2)) ...
                       - ROItac ).*weighs2);
                   
                    [parest_srtm2, minnorm_srtm2] = fminsearch(@(x) fobj_srtm2(x), [1 2], options);
                    R1_srtm2 = parest_srtm2(1);
                    BP_srtm2 = parest_srtm2(2);
                    
                    % for plotting or saving, generate the SRTM2 fit:
                    modfit_srtm2 = simSRTM2_1_0_0(tmidMin, reftac, t_points, ...
                                                  k2_ref, R1_srtm2, BP_srtm2);

                    % save the SRTM2 outcome in BPdataSave:
                    eval(['BPdataSave{id,d}.' savestr '.BP_srtm2 = BP_srtm2;']);
                    eval(['BPdataSave{id,d}.' savestr '.R1_srtm2 = R1_srtm2;']);
                    eval(['Residuals{id,d}.' savestr '.AllFit_SRTM2 = ROItac - modfit_srtm2;']);
                    %  ============================================================ %


                    % ============ jarkko's original lp-ntPET =========== %

                    %%% do lp-ntPET !!!!
                    Alpntpet=zeros(t_points,4);
                    Alpntpet(:,1:3)=ASRTM;
                    best_mse=10^20;
                    best_parest=[];
                    best_se=[];
                    
                    % for estimating the optimal peak and alpha parameters
                    for point_rise=break_point:length(tmid)-1
                        t2p_index = find(tmid>tmid(point_rise));
                        if length(t2p_index)>1
                            for t_ind=t2p_index(1):t2p_index(end-1) 
                                for alpha=[0.25 1 4]
                                    t_peak=tmid(t_ind)-tmid(point_rise);
                                    p = [1 alpha tmid(point_rise)+t_peak tmid(point_rise)];
                                    actfun = zeros(size(tmid));
                                    actfun(point_rise:t_points) = gamma_variate_madsen(p,tmid(point_rise:t_points));
                                    roitac_gamma = ROItac.*actfun;
                                    mroitac_gamma  = [roitac_gamma(1)/2; (roitac_gamma(2:length(ROItac))+roitac_gamma(1:length(ROItac)-1))/2];
                                    
                                    Alpntpet(:,4)=0;
                                    for k = break_point:t_points
                                        Alpntpet(k,4)  = -sum(mroitac_gamma(break_point:k).*dt(break_point:k));
                                    end
                                    
                                    %LSQ-estimation using lscov
                                    [parest se mse]   = lscov(Alpntpet,ROItac);
                                    
                                    %estimated model TAC
                                    modelfit = Alpntpet*parest;
                                    
                                    %%%%%% magic (finding the best model fit) %%%%%
                                    if (best_mse > mse)
                                        best_mse = mse;
                                        best_parest=parest; % R1 - k2 - k2a - gamma
                                        best_modelfit=modelfit;
                                        best_actfun=actfun;
                                        best_se=se;
                                    end
                                end
                            end
        %                        breakpoint2=breakpoint;
                        end
                    end
                    
                    y=-ASRTM(:,3)./ROItac;
                    x=ASRTM(:,2)./ROItac;
                    [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                    [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                    
                    % save these results
                    k2=best_parest(2);
                    k2a=best_parest(3);
                    BP_lp=k2/k2a-1; % baseline binding potential of lpntPET
                    BP_srtm=BP__; % binding potential until the end of the task
                    BP_srtm_bsl=BP_bl; % binding potential until the end of baseline
                    G=best_parest(4); % gamma
                    BPND = ((R1__*k2p)/k2a)-1;
%                     DBP =  k2./(k2a + best_actfun) - 1; % dynamic binding potential
                    DBP =  k2./(k2a + G*best_actfun) - 1; % is it this way...? confused
                    OCC = 100*(1-DBP/BP_lp);% receptor occupancy
                    eval(['BP_lp_save{id,d}.' savestr '=BP_lp;' ])
                    eval(['DBP_save{id,d}.' savestr '=DBP;' ])
                    eval(['Occupancy{id,d}.' savestr '=OCC;' ])
                    eval(['BPND_save{id,d}.' savestr '=BPND;' ])
                    eval(['BP_srtm_save{id,d}.' savestr '=BP_srtm;' ])
                    eval(['BP_srtm2_save{id,d}.' savestr '=BP_srtm2;' ])
                    eval(['BP_srtm_Bsl_save{id,d}.' savestr '=BP_srtm_bsl;' ])
                    eval(['BestActivationFunction{id,d}.' savestr '=best_actfun;'])
                    eval(['BestModelFit{id,d}.' savestr '=best_modelfit;'])
                    eval(['gamma_lp_save{id,d}.' savestr '=best_parest(4);' ])
                    eval(['Residuals{id,d}.' savestr '.BaselineFit=ROItac-fittac;' ])
                    eval(['Residuals{id,d}.' savestr '.AllFit_ESRTM=ROItac-modfit_esrtm;' ])
                    eval(['Residuals{id,d}.' savestr '.lpntPETFit=ROItac-best_modelfit;' ])
                    eval(['CompensatoryFunction{id,d}.' savestr '=best_parest(4)*best_actfun;' ])
                    
                  
                    BPdata.BP_mrtm(Subj{r})=BP;
                    BPdata.BP_srtm(Subj{r})=BP__;
                    BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                    BPdata.BP_lpnt(Subj{r})=BP_lp;
                    BPdata.BP_logan(Subj{r})=pp(1)-1;
                    
                    % save
                    BPdataSave{id,d}.(reg{1}).BP_mrtm=BP;
                    BPdataSave{id,d}.(reg{1}).BP_srtm=BP__;
                    BPdataSave{id,d}.(reg{1}).BP_srtm_bl=BP_bl;
                    BPdataSave{id,d}.(reg{1}).BP_lpnt=BP_lp;
                    BPdataSave{id,d}.(reg{1}).BP_logan=pp(1)-1;
                    BPdataSave{id,d}.(reg{1}).BPND=BPND;
                    BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
                    BPdataSave{id,d}.(reg{1}).DBP=DBP;
                    BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);


                    if PlotStrFit
                        legendtext{end+1}=[reg{1} ' raw'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{Baseline}'];
                        legendtext{end+1}=[reg{1} ' fit SRTM_{All}'];
                        legendtext{end+1}=[reg{1} ' fit lpnt-PET_{All}'];
                        legendtext{end+1}=[reg{1} ' fit SRTM2'];
                        legendtext{end+1}=['Activation function'];
                        SE=abs(BP)*sqrt((se_srtm(2)/parest(2))^2+(se_srtm(3)/parest(3))^2);
                        subplot(3,2,[1 2]);
                        yyaxis left;
                        plot(tmidMin,ROItac,'ro', ...
                            tmidMin,modfit_esrtm_bl,'r-', ...
                            tmidMin,modfit_esrtm,'b-',...
                            tmidMin,best_modelfit,'k--',...
                            tmidMin, modfit_srtm2, 'g-'); hold on;
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
                        plot(tmidMin,ROItac-fittac,'yo',...
                            tmidMin,ROItac-modfit_esrtm,'bo',...
                            tmidMin,ROItac-modfit_esrtm_bl,'ro',...
                            tmidMin,ROItac-best_modelfit,'ko',...
                            [0 180],[0 0],'k--'); hold on;
                        plot(tmidMin, ROItac - modfit_srtm2, 'go');

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
                        [pp ss]=polyfit(x(end-11:end),y(end-11:end),1);
                        [pp2 ss2]=polyfit(x(end-25:end-11),y(end-25:end-11),1);
                        plot(x,y,'ko',x(end-25:end),polyval(pp,x(end-25:end)),'k-',x(end-25:end),polyval(pp2,x(end-25:end)),'k--');
                        title(['Logan fit: BP(Baseline)=' num2str(pp2(1)-1,'%1.2f') ' BP(Task)=' num2str(pp(1)-1,'%1.2f') ]);
                        xlabel(['\int REF/ROI']);
                        ylabel('\int ROI/ROI')
                        subplot(3,2,[5 6]);
                        plot(tmidMin,ROItac./reftac,'ko',[0 180],[0 0],'k--');
                        xlim([0 190]); %ylim([0.5 30]);
                        title('Target to reference ratio');
                        xlabel('Time (min)');

                        print('-dpsc2','-append','-bestfit',fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_TAC_Fit_SRTM2_logan_' date '_s3.ps']));
                        savefig(fullfile(paths.figure, [ num2str(IDs(id)) num2str(d) '_' reg{1} '_TAC_Fit_SRTM2_logan_' date '_s3.fig']))
                        close(gcf)
                        BPdata.BP_mrtm(Subj{r})=BP;
                        BPdata.BP_srtm(Subj{r})=BP__;
                        BPdata.BP_srtm_bl(Subj{r})=BP_bl;
                        BPdata.BP_lpnt(Subj{r})=BP_lp;
                        BPdata.BP_logan(Subj{r})=pp(1)-1;
                        
                        % save
                        BPdataSave{id,d}.(reg{1}).BP_mrtm=BP;
                        BPdataSave{id,d}.(reg{1}).BP_srtm=BP__;
                        BPdataSave{id,d}.(reg{1}).BP_srtm_bl=BP_bl;
                        BPdataSave{id,d}.(reg{1}).BP_lpnt=BP_lp;
                        BPdataSave{id,d}.(reg{1}).BP_logan=pp(1)-1;
                        BPdataSave{id,d}.(reg{1}).BPND=BPND;
                        BPdataSave{id,d}.(reg{1}).Occupancy=OCC;
                        BPdataSave{id,d}.(reg{1}).DBP=DBP;
                        BPdataSave{id,d}.(reg{1}).gamma=best_parest(4);
                        
                        
                        continue;
                    end
                end
                
                
            end
            close all
            keep IDs days paths id d TACs BPdata Occupancy gamma_lp_save BP_lp_save DBP_save ROI BPdataSave BP_lp_save DBP_save Occupancy BPND_save BP_srtm_save BP_srtm2_save BP_srtm_Bsl_save BestActivationFunction...
    BestModelFit gamma_lp_save Residuals CompensatoryFunction
        end
    end
end

disp('done')

save(['/Users/alex/Dropbox/paperwriting/MRPET/data/MRPET_BPpackage_' date '_newRecon_noPVC_SRTM2.mat'],...
    'BPdataSave', 'BP_lp_save', 'DBP_save', 'Occupancy', 'BPND_save', 'BP_srtm_save', 'BP_srtm_Bsl_save', 'BestActivationFunction',...
    'BestModelFit', 'gamma_lp_save', 'Residuals', 'CompensatoryFunction')