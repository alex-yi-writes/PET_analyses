%% new PET registrations

clc;clear;

path_parent = '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/';
path_segs   = '/Users/alex/Dropbox/paperwriting/MRPET/data/new_recon/segs/';

% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 0; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];

%% copy files to the workspace

for id=1:length(IDs)
    
    for d=1:2
        
%         mkdir([path_parent 'segs/' num2str(IDs(id)) num2str(d)])
        
        copyfile(['/Volumes/korokdorf/MRPET/coreg_roi/' num2str(IDs(id)) num2str(d) '/data/T1pt1.nii'],...
            [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii'])
        copyfile(['/Volumes/korokdorf/MRPET/coreg_roi/' num2str(IDs(id)) num2str(d) '/data/T1pt2.nii'],...
            [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii'])
        
        
%         % first baseline absolute
%         clear tmp
%         tmp=dir([path_parent 'raw/' num2str(IDs(id)) num2str(d) '/*/*PET_B_Sa_*.nii']);
%         copyfile(fullfile(tmp.folder,tmp.name),...
%             [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl_absolute.nii'])
%         % baseline relative
%         clear tmp
%         tmp=dir([path_parent 'raw/' num2str(IDs(id)) num2str(d) '/*/*PET_B_Sr_*.nii']);
%         copyfile(fullfile(tmp.folder,tmp.name),...
%             [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl_relative.nii'])
%         
%         % task absolute
%         clear tmp
%         tmp=dir([path_parent 'raw/' num2str(IDs(id)) num2str(d) '/*/*PET_T_Sa_*.nii']);
%         copyfile(fullfile(tmp.folder,tmp.name),...
%             [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task_absolute.nii'])
%         % task absolute
%         clear tmp
%         tmp=dir([path_parent 'raw/' num2str(IDs(id)) num2str(d) '/*/*PET_T_Sr_*.nii']);
%         copyfile(fullfile(tmp.folder,tmp.name),...
%             [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task_relative.nii'])


    end
    
end


%% set env

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin']);
setenv('ANTSPATH','/usr/local/antsbin/bin')

%% convert seg files

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            

            movefile([path_segs num2str(IDs(id)) num2str(d) '1/T1WB.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'])
            movefile([path_segs num2str(IDs(id)) num2str(d) '2/T1WB.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'])
            movefile([path_segs num2str(IDs(id)) num2str(d) '1/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'])
            movefile([path_segs num2str(IDs(id)) num2str(d) '2/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'])
        end
        
    end
end

%% extra

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            
            mkdir(['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d)])
            copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt2.nii'],...
                ['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d) '/T1pt2.nii'])
            copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2_nat_labelled.nii'],...
                ['/Volumes/ALEX3/MRPET/segs/' num2str(IDs(id)) num2str(d) '/aparc+aseg_pt2_nat_labelled.nii'])
        end
        
    end
end

%% move files for ANTs coreg

% for id=1:length(IDs)
%     
%     for d=1:2
%         if days(id,d)==0
%             warning('skipped')
%         else
%             
%             mkdir(['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/'])
% 
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/T1pt1.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/T1pt2.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/aparc+aseg_pt1.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/aparc+aseg_pt2.nii'])
%             copyfile([path_parent num2str(IDs(id)) '_' num2str(d) '/meana' num2str(IDs(id)) '_MRI_4D_MT' num2str(d) '.nii'], ['/Volumes/MY PASSPORT/MRPET/' num2str(IDs(id)) num2str(d) '/data/meanEPI.nii'])
%         end
%         
%     end
% end


%% convert seg files

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            

            movefile([path_segs num2str(IDs(id)) num2str(d) '1/T1.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt1.nii'])
            movefile([path_segs num2str(IDs(id)) num2str(d) '2/T1.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/T1pt2.nii'])
            movefile([path_segs num2str(IDs(id)) num2str(d) '1/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt1.nii'])
            movefile([path_segs num2str(IDs(id)) num2str(d) '2/aparc+aseg.nii'], [path_parent num2str(IDs(id)) '_' num2str(d) '/aparc+aseg_pt2.nii'])
            
        end
        
    end
end



%% first, realign and create

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            %% inflow

            list_inflow=[]; clear numvol
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
            for v1=10:numvol
                list_inflow{v1-9,1} = [path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii,' num2str(v1)];
            end

            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_inflow};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            spm_jobman('run',matlabbatch)


            mkdir([path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D'])

            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D']};
            spm_jobman('run',matlabbatch)

            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent num2str(IDs(id)) '_' num2str(d) '/mean' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii -f '...
                path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlowtoT1pt1_' ])
            % apply transformations to the images
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlowtoT1pt1_0GenericAffine.mat -i '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '_' sprintf('%05i',v2) '.nii -r '...
                    path_parent num2str(IDs(id)) '_' num2str(d) '/' num2str(IDs(id)) '_MRI_4D_MPRAGE' num2str(d) '_pt1.nii -o ' path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1_' sprintf('%05i',v2) '.nii'])
            end

            % assemble again
            list_inflow=[];
            numvol = length(spm_vol([path_parent num2str(IDs(id)) '_' num2str(d) '/r' num2str(IDs(id)) '_PET_4D_InFlow' num2str(d) '.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1_' sprintf('%05i',v1) '.nii,1'];
            end
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_inflow;
            matlabbatch{1}.spm.util.cat.name = ['coreg_InFlow' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)

            movefile([path_parent num2str(IDs(id)) '_' num2str(d)  '/inflow3D/coreg_InFlow' num2str(d) '_on_T1.nii'],...
                [path_parent num2str(IDs(id)) '_' num2str(d)  '/coreg_InFlow' num2str(d) '_on_T1.nii'])
            rmdir([path_parent num2str(IDs(id)) '_' num2str(d) '/inflow3D'],'s')
            
            %% realign
            
            list_inflow=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow_absolute.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow_absolute.nii,'  num2str(v1)];
            end
            
            % list_bsl=[];
            % numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl_absolute.nii']));
            % for v1=1:numvol
            %     list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl_absolute.nii,'  num2str(v1)];
            % end
            % 
            % list_task=[];
            % numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task_absolute.nii']));
            % for v1=1:numvol
            %     list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task_absolute.nii,'  num2str(v1)];
            % end
            % 
            clear matlabbatch
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {list_inflow};
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
            matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
%             matlabbatch{2}.spm.spatial.realign.estwrite.data = {list_task};
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 1;
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
%             matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
%             matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
%             matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
%             matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
%             matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
%             matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
            
            spm_jobman('run',matlabbatch)

            %% inflow
            
            list_inflow=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rinflow_absolute.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rinflow_absolute.nii,'  num2str(v1)];
            end
                        
            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D'])
            
            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rinflow_absolute.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D']};
            spm_jobman('run',matlabbatch)
            
            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/meaninflow_absolute.nii -f '...
                path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_InflowtoT1pt1_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_InflowtoT1pt1_0GenericAffine.mat -i '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D/rinflow_absolute_' sprintf('%05i',v2) '.nii -r '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt1.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/inflow3D/crinflow_absolute_' sprintf('%05i',v2) '.nii'])
            end
            
            % assemble again
            list_inflow=[];
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/inflow3D/crinflow_absolute_' sprintf('%05i',v1) '.nii,1'];
            end
            
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_inflow;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Inflow' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)
            
            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/inflow3D/coreg_Inflow' num2str(d) '_on_T1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Inflow' num2str(d) '_on_T1.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/inflow3D'],'s')

            %% baseline
            list_bsl=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rbsl_absolute.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rbsl_absolute.nii,'  num2str(v1)];
            end

            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D'])

            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rbsl_absolute.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D']};
            spm_jobman('run',matlabbatch)

            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/meanbsl_absolute.nii -f '...
                path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_BaselinetoT1pt2_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_BaselinetoT1pt2_0GenericAffine.mat -i '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D/rbsl_absolute_' sprintf('%05i',v2) '.nii -r '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/bsl3D/crbsl_absolute_' sprintf('%05i',v2) '.nii'])
            end

            % assemble again
            list_bsl=[];
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/bsl3D/crbsl_absolute_' sprintf('%05i',v1) '.nii,1'];
            end

            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_bsl;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Baseline' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)

            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/bsl3D/coreg_Baseline' num2str(d) '_on_T1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Baseline' num2str(d) '_on_T1.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/bsl3D'],'s')


            %% task
            list_task=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rtask_absolute.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rtask_absolute.nii,'  num2str(v1)];
            end

            mkdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D'])

            % unpack 4d for registration
            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.split.vol = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/rtask_absolute.nii']};
            matlabbatch{1}.spm.util.split.outdir = {[path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D']};
            spm_jobman('run',matlabbatch)

            % register meanPET to T1
            eval(['!antsRegistrationSyN.sh -d 3 -t r -m ' path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/meantask_absolute.nii -f '...
                path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_TasktoT1pt2_' ])
            % apply transformations to the images
            for v2=1:numvol
                eval(['!antsApplyTransforms -d 3 -v 0 -n Linear -t ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_TasktoT1pt2_0GenericAffine.mat -i '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D/rtask_absolute_' sprintf('%05i',v2) '.nii -r '...
                    path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/T1pt2.nii -o ' path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/task3D/crtask_absolute_' sprintf('%05i',v2) '.nii'])
            end

            % assemble again
            list_task=[];
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/task3D/crtask_absolute_' sprintf('%05i',v1) '.nii,1'];
            end

            clear matlabbatch
            spm_jobman('initcfg')
            matlabbatch{1}.spm.util.cat.vols = list_task;
            matlabbatch{1}.spm.util.cat.name = ['coreg_Task' num2str(d) '_on_T1.nii'];
            matlabbatch{1}.spm.util.cat.dtype = 4;
            spm_jobman('run',matlabbatch)

            movefile([path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/task3D/coreg_Task' num2str(d) '_on_T1.nii'],...
                [path_parent 'preproc/' num2str(IDs(id)) num2str(d)  '/coreg_Task' num2str(d) '_on_T1.nii'])
            rmdir([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/task3D'],'s')
            
            
        end
    end
    
end

%% smoothe

for id=1:length(IDs)
    
    for d=1:2
        if days(id,d)==0
            warning('skipped')
        else
            list_inflow=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Inflow' num2str(d) '_on_T1.nii']));
            for v1=1:numvol
                list_inflow{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Inflow' num2str(d) '_on_T1.nii,'  num2str(v1)];
            end
            list_bsl=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Baseline' num2str(d) '_on_T1.nii']));
            for v1=1:numvol
                list_bsl{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Baseline' num2str(d) '_on_T1.nii,'  num2str(v1)];
            end
            list_task=[];
            numvol = length(spm_vol([path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Task' num2str(d) '_on_T1.nii']));
            for v1=1:numvol
                list_task{v1,1} = [path_parent 'preproc/' num2str(IDs(id)) num2str(d) '/coreg_Task' num2str(d) '_on_T1.nii,'  num2str(v1)];
            end
            
            clear matlabbatch
            matlabbatch{1}.spm.spatial.smooth.data = list_inflow;
            matlabbatch{1}.spm.spatial.smooth.fwhm = [3 3 3];
            matlabbatch{1}.spm.spatial.smooth.dtype = 0;
            matlabbatch{1}.spm.spatial.smooth.im = 0;
            matlabbatch{1}.spm.spatial.smooth.prefix = 's3_';
            matlabbatch{2}.spm.spatial.smooth.data = list_bsl;
            matlabbatch{2}.spm.spatial.smooth.fwhm = [3 3 3];
            matlabbatch{2}.spm.spatial.smooth.dtype = 0;
            matlabbatch{2}.spm.spatial.smooth.im = 0;
            matlabbatch{2}.spm.spatial.smooth.prefix = 's3_';
            matlabbatch{3}.spm.spatial.smooth.data = list_task;
            matlabbatch{3}.spm.spatial.smooth.fwhm = [3 3 3];
            matlabbatch{3}.spm.spatial.smooth.dtype = 0;
            matlabbatch{3}.spm.spatial.smooth.im = 0;
            matlabbatch{3}.spm.spatial.smooth.prefix = 's3_';
            spm_jobman('run',matlabbatch)
            
        end
    end
end



