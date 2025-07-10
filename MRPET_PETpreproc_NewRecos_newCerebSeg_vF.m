%% new PET registrations

clc;clear;

path_parent = '/Volumes/korokdorf/MRPET/PVEc_new/';
path_segs   = '/Volumes/korokdorf/MRPET/coreg_roi/';
path_tacs   = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TACs/NewRecon/noPVC/';
path_fig    = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TACs/NewRecon/tacfig/';

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];


%% copy new segs

for id = 11:length(IDs)

    for d = 1:2

        if days(id,d) == 0
            warning('Skipped')
        else
            copyfile(['/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/newseg/' num2str(IDs(id)) num2str(d) '1/aparc+aseg_pt1_nat_labelled_updated.nii'],...
                [path_segs  num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled_updated.nii'])
            copyfile(['/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/newseg/' num2str(IDs(id)) num2str(d) '2/aparc+aseg_pt2_nat_labelled_updated.nii'],...
                [path_segs  num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled_updated.nii'])

        end
    end
end

%% set env

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin']);
setenv('ANTSPATH','/usr/local/antsbin/bin')

%% write out TACs

opengl hardwarebasic

% percentage of radioactivity concentrations trimmed-out when calculated
% ROI-average

TrimPerc=15;

clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab

            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask1   =[ path_segs num2str(IDs(id)) num2str(d) '/aparc+aseg_pt1_nat_labelled_updated.nii'];
            
            % 4-D PET file
            PETflow = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/inflow_on_t1.nii.gz'];
            
            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask1);
            ROIMask=round(ROI.img);
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            % Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data

            [A B ROIDef]=xlsread('/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/scripts/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
            maskidx=[];
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,2}).LongName=strtrim(ROIDef{ROITabidx,3});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1}+1000; % Add 1000 to the index
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,1}+2000; % Add 2000 to the index
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}+1000, ROIDef{ROITabidx,1}+2000];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,2}).(Hemi{1})=IndList;
                end
            end
            
            
            [A B ROIDef]=xlsread('/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/scripts/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,3}).LongName=strtrim(ROIDef{ROITabidx,4});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1};
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,2};
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}, ROIDef{ROITabidx,2}];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,3}).(Hemi{1})=IndList;
                end
            end
            
            %% Read in 4D-PET data and extract ROI-averages in each frame
            
            % inFlow
            DynPET=load_nii(PETflow);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                disp(maskidx.(ROI{1}).LongName)
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'inFlow';
            
            mkdir([ path_tacs num2str(IDs(id)) num2str(d) ])
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_noPVC_newCerC.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData            
        end
    end
end

clear id d

for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            % Freesurfer segmentation, if .mgh use mri_read from FreeSurfer/Matlab
            clear Mask CurPET_task CurPET_flow CurPET_BSL
            Mask2   =[ path_segs num2str(IDs(id)) num2str(days(id,d)) '/aparc+aseg_pt2_nat_labelled_updated.nii']; % task and the baseline
            
            % 4-D PET file
            if (id==8 & d==1) | (id==21 && d==1)
                PETtask = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/rtask_on_t1.nii.gz'];
                PETbsl  = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/rbsl_on_t1.nii.gz'];
            elseif (id==13 && d==2) | (id==15 && d==1) | (id==16 && d==1) | (id==28 && d==2)
                PETtask = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/rtask_on_t1.nii'];
                PETbsl  = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/bsl_on_t1.nii.gz'];
            else
                PETtask = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/task_on_t1.nii.gz'];
                PETbsl  = [path_parent num2str(IDs(id)) num2str(days(id,d)) '/bsl_on_t1.nii.gz'];
            end

            %% Read in FreeSurfer mask data
            ROI=load_nii(Mask2);
            ROIMask=round(ROI.img); 
            
            %% Read in Dictionary for Desikan-Killiany atlas and find voxel coordinates from this individual mask
            %  Voxel indices are first collected in a structure (maskidx) and later applied to extract time-activity data
            [A B ROIDef]=xlsread('/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/scripts/ExtractPETTACs/Dictionary_FSaparc2004_Desikan_ROIs.xls');
            maskidx=[];
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,2}).LongName=strtrim(ROIDef{ROITabidx,3});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1}+1000; % Add 1000 to the index
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,1}+2000; % Add 2000 to the index
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}+1000, ROIDef{ROITabidx,1}+2000];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,2}).(Hemi{1})=IndList;
                end
            end
            
            
            [A B ROIDef]=xlsread('/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/scripts/ExtractPETTACs/Subcortical_Dictionary_2.xls');
            
            for ROITabidx=2:length(ROIDef)
                maskidx.(ROIDef{ROITabidx,3}).LongName=strtrim(ROIDef{ROITabidx,4});
                for Hemi={'Left','Right','Bilateral'}
                    switch Hemi{1}
                        case 'Left'
                            CurIdx=ROIDef{ROITabidx,1};
                        case 'Right'
                            CurIdx=ROIDef{ROITabidx,2};
                        case 'Bilateral'
                            CurIdx=[ROIDef{ROITabidx,1}, ROIDef{ROITabidx,2}];
                    end
                    IndList=[];
                    for ROIInd=CurIdx
                        if ~isnan(ROIInd)
                            IndList=unique(union(IndList,find(ROIMask == ROIInd)));
                        end
                    end
                    maskidx.(ROIDef{ROITabidx,3}).(Hemi{1})=IndList;
                end
            end
            
            %% Read in 4D-PET data and extract ROI-averages in each frame
            %% task
            DynPET=load_nii(PETtask);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'RewardTask';
            if (id==8 & d==1) | (id==21 && d==1)
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_noPVC_corrected_newCerC.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            else
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_noPVC_newCerC.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            end
            
            %% baseline
            DynPET=load_nii(PETbsl);
            temp=size(DynPET.img);
            ImgData=reshape(DynPET.img,prod(temp(1:3)),temp(4));
            for ROI=fieldnames(maskidx)'
                TACDATA.(ROI{1}).LongName=maskidx.(ROI{1}).LongName;
                for Hemi={'Left','Right','Bilateral'}
                    maskidxcur = maskidx.(ROI{1}).(Hemi{1});
                    TACDATA.(ROI{1}).(Hemi{1}).vol=length(maskidxcur)*(1-TrimPerc/100);
                    TACDATA.(ROI{1}).(Hemi{1}).tac=trimmean(ImgData(maskidxcur,:),TrimPerc)';
                end
                
            end
            TACDATA.info = 'Baseline';
            if (id==8 & d==1) | (id==21 && d==1)
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_noPVC_corrected_newCerC.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            else
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_noPVC_newCerC.mat'],'TACDATA'); clear TACDATA DynPET temp ImgData
            end
        end
    end
end

%% originals

close all

t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)

    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_noPVC_newCerC.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            
            if (id==8 & d==1) | (id==21 && d==1)
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_noPVC_corrected_newCerC.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_noPVC_corrected_newCerC.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            else
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_noPVC_newCerC.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_noPVC_newCerC.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            end
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    % bilaterals
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac);
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac);
                    
                    % lefts
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac);
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac);
                    
                    % rights
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac);
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac);
                end
            end
            
            if (id==8 & d==1) | (id==21 && d==1)
            save([path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_noPVC_corrected_newCerC.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            else
            save([path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_noPVC_newCerC.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            end
            
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
            try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                if (id==8 & d==1) | (id==21 && d==1)
                print('-dpdf','-bestfit',[ path_fig num2str(IDs(id)) num2str(d) '_decayUncorr_noPVC_corrected_newCerC.pdf']);
                else
                print('-dpdf','-bestfit',[ path_fig num2str(IDs(id)) num2str(d) '_decayUncorr_noPVC_newCerC.pdf']);
                end
                
            catch
                
                Cer=[vertcat(TACDATA_InFlow.CerC.Bilateral.tac,nan(14,1)); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[vertcat(TACDATA_InFlow.Put.Bilateral.tac,nan(14,1)); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[vertcat(TACDATA_InFlow.Caud.Bilateral.tac,nan(14,1)); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                
                if (id==8 & d==1) | (id==21 && d==1)
                print('-dpdf','-bestfit',[ path_fig num2str(IDs(id)) num2str(d) '_decayUncorr_noPVC_corrected_newCerC.pdf']);
                else
                print('-dpdf','-bestfit',[ path_fig num2str(IDs(id)) num2str(d) '_decayUncorr_noPVC_newCerC.pdf']);
                end

                
            end
            
            
        end
    end
end


%% retro correct the decay

close all

t0_frame_bsl    = 95;
t0_frame_task   = 115;
 
for id = 1:length(IDs)
    
    for d = 1:2
        
        if days(id,d) == 0
            warning('Skipped')
        else
            
            if (id==8 & d==1) | (id==21 && d==1)
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_noPVC_newCerC.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_noPVC_corrected_newCerC.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_noPVC_corrected_newCerC.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            else
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_InFlow_noPVC_newCerC.mat']);
            TACDATA_InFlow=TACDATA; clear TACDATA
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Baseline_noPVC_newCerC.mat']);
            TACDATA_Baseline=TACDATA; clear TACDATA
            load([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACDATA_Task_noPVC_newCerC.mat']);
            TACDATA_Task=TACDATA; clear TACDATA
            end
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60];
            
            
            clear fields
            fields = fieldnames(TACDATA_InFlow);
            for f1 = 1:length(fields)
                if strcmp((fields{f1}),'info')
                    disp('info field skipped')
                else
                    
                    % bilaterals
                    if IDs(id)==4008 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_bsl+4)/109));
                    elseif IDs(id)==4020 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs) 
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_bsl-17)/109));
                    elseif (IDs(id)==4026 && d==1) | (IDs(id)==4026 && d==2) % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_bsl-10)/109));
                    else
                    TACDATA_Baseline.(fields{f1}).Bilateral.tac=(TACDATA_Baseline.(fields{f1}).Bilateral.tac).*(2^(t0_frame_bsl/109));
                    end
                    
                    if IDs(id)==4020 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_task+24)/109 ));
                    elseif IDs(id)==4020 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_task-23)/109 ));
                    elseif (IDs(id)==4026 && d==1) | (IDs(id)==4026 && d==2) % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_task-10)/109 ));
                    elseif IDs(id)==4027 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_task+6)/109 ));
                    elseif IDs(id)==4008 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_task+3)/109 ));
                    elseif IDs(id)==4032 && d==2 % first task frame dropped
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^( (t0_frame_task+5)/109 )); 
                    else
                    TACDATA_Task.(fields{f1}).Bilateral.tac=(TACDATA_Task.(fields{f1}).Bilateral.tac).*(2^(t0_frame_task/109));
                    end
                    
                    % lefts
                    if IDs(id)==4008 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^( (t0_frame_bsl+4)/109));
                    elseif IDs(id)==4020 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^( (t0_frame_bsl-17)/109));
                    elseif (IDs(id)==4026 && d==1) | (IDs(id)==4026 && d==2) % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^( (t0_frame_bsl-10)/109));
                    else
                    TACDATA_Baseline.(fields{f1}).Left.tac=(TACDATA_Baseline.(fields{f1}).Left.tac).*(2^(t0_frame_bsl/109));
                    end
                    
                    if IDs(id)==4020 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^( (t0_frame_task+24)/109 ));
                    elseif IDs(id)==4020 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^( (t0_frame_task-23)/109 ));
                    elseif (IDs(id)==4026 && d==1) | (IDs(id)==4026 && d==2) % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^( (t0_frame_task-10)/109 ));
                    elseif IDs(id)==4027 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^( (t0_frame_task+6)/109 ));
                    elseif IDs(id)==4008 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^( (t0_frame_task+3)/109 ));
                    elseif IDs(id)==4032 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^( (t0_frame_task+5)/109 ));    
                    else
                    TACDATA_Task.(fields{f1}).Left.tac=(TACDATA_Task.(fields{f1}).Left.tac).*(2^(t0_frame_task/109));
                    end
                    
                    % rights
                    if IDs(id)==4008 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^( (t0_frame_bsl+4)/109));
                    elseif IDs(id)==4020 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^( (t0_frame_bsl-17)/109));
                    elseif (IDs(id)==4026 && d==1) | (IDs(id)==4026 && d==2) % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^( (t0_frame_bsl-10)/109));
                    else
                    TACDATA_Baseline.(fields{f1}).Right.tac=(TACDATA_Baseline.(fields{f1}).Right.tac).*(2^(t0_frame_bsl/109));
                    end
                    
                    if IDs(id)==4020 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^( (t0_frame_task+24)/109 ));
                    elseif IDs(id)==4020 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^( (t0_frame_task-23)/109 ));
                    elseif (IDs(id)==4026 && d==1) | (IDs(id)==4026 && d==2) % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^( (t0_frame_task-10)/109 ));
                    elseif IDs(id)==4027 && d==1 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^( (t0_frame_task+6)/109 ));
                    elseif IDs(id)==4008 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^( (t0_frame_task+3)/109 ));
                    elseif IDs(id)==4032 && d==2 % started at different time, correcting with the actual timing (refer the experimental logs)
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^( (t0_frame_task+5)/109 ));    
                    else
                    TACDATA_Task.(fields{f1}).Right.tac=(TACDATA_Task.(fields{f1}).Right.tac).*(2^(t0_frame_task/109));
                    end
                    
                end
            end
            
            if (id==8 & d==1) | (id==21 && d==1)
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_noPVC_AdjustedActualStart_newCerC.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            else
            save([ path_tacs num2str(IDs(id)) num2str(d) '/' num2str(IDs(id)) num2str(d) '_TACs_deccorr_noPVC_newCerC.mat'],'TACDATA_Baseline','TACDATA_InFlow','TACDATA_Task');
            end
            
            Lengths=[60*ones(60,1)];
            tt1=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(15,1);
            tt2=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Lengths=60*ones(55,1);
            tt3=[[0;cumsum(Lengths(1:end-1))], cumsum(Lengths)]; clear Lengths
            Times=[tt1; tt2+95*60; tt3+115*60]
            try
                
                Cer=[TACDATA_InFlow.CerC.Bilateral.tac; TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[TACDATA_InFlow.Put.Bilateral.tac; TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[TACDATA_InFlow.Caud.Bilateral.tac; TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                if (id==8 & d==1) | (id==21 && d==1)
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_noPVC_AdjustedActualStart_newCerC.pdf']);
                else
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_noPVC_newCerC.pdf']);
                end
                
            catch
                
                Cer=[vertcat(TACDATA_InFlow.CerC.Bilateral.tac,nan(14,1)); TACDATA_Baseline.CerC.Bilateral.tac; TACDATA_Task.CerC.Bilateral.tac];
                Put=[vertcat(TACDATA_InFlow.Put.Bilateral.tac,nan(14,1)); TACDATA_Baseline.Put.Bilateral.tac; TACDATA_Task.Put.Bilateral.tac];
                Caud=[vertcat(TACDATA_InFlow.Caud.Bilateral.tac,nan(14,1)); TACDATA_Baseline.Caud.Bilateral.tac; TACDATA_Task.Caud.Bilateral.tac];
                tmid=mean(Times,2)/60;
                
                % now draw
                figure('Renderer', 'painters ')
                plot(tmid,Cer,'ko-',tmid,Put,'ro-',tmid,Caud,'bo-');
                xlabel('Time (min)'); ylabel('Radioactivity (Bq/mL)');
                legend('Cerebellum','Putamen','Caudate');
                title([ num2str(IDs(id)) '_' num2str(d) ]);
                ax = gca; ax.YAxis.Exponent = 0;
                if (id==8 & d==1) | (id==21 && d==1)
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_noPVC_AdjustedActualStart_newCerC.pdf']);
                else
                print('-dpdf','-bestfit',[ num2str(IDs(id)) num2str(d) '_decaycorrected_noPVC_newCerC.pdf']);
                end

                
        end
            
            
        end
    end
end
