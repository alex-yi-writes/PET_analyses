%% calc TIV

%% prep

clc; clear

% paths

path_analysis   = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TIV/t1/';
path_source     = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TIV/t1/';


% IDs
IDs  = [4001 4002 4003 4004 4005 4006 4007 4008 4009 4010 4011 4012 4013 4014 4015 4016 4017 4018 4019 4020 4021 4022 4023 4024 4025 4026 4027 4028 4029 4030 4031 4032 4033];
days = [1 2; 1 2; 1 0; 1 2; 1 2; 0 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 1 0; 1 2; 0 2; 1 2; 1 2; 1 2; 1 2];


cnt=0;
for i1=1:length(IDs)
    
    mkdir([path_analysis num2str(IDs(i1))])

    for d = 1:2
        if days(i1,d) == 0
        else
            cnt=cnt+1;
            ID{cnt}=[num2str(IDs(i1)) num2str(d)];
            
            
            copyfile([path_source num2str(IDs(i1)) '_' num2str(d) '/T1pt1.nii'],...
                [path_analysis num2str(IDs(i1)) '/T1pt1_' num2str(d) '.nii'])
            copyfile([path_source num2str(IDs(i1)) '_' num2str(d) '/T1pt1.nii'],...
                [path_analysis num2str(IDs(i1)) '/T1pt2_' num2str(d) '.nii'])
            
        end
    end
end


setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/usr/local/antsbin/bin']);
setenv('ANTSPATH','/usr/local/antsbin/bin')

%%
for id=1:length(IDs)
    
    clear tmp
    tmp=dir([path_analysis num2str(IDs(id)) '/T1*.nii']);
    
%     for tt=2:length(tmp)
%     eval(['!antsRegistrationSyNQuick.sh -d 3 -n 4 -t r -o ' path_analysis num2str(IDs(id)) '/T1' num2str(tt) '_ -f ' path_analysis num2str(IDs(id)) '/' tmp(1).name ' -m ' path_analysis num2str(IDs(id)) '/' tmp(tt).name]);
%     end
%     
    eval(['!AverageImages 3 ' path_analysis num2str(IDs(id)) '/T1avg.nii 0 ' path_analysis num2str(IDs(id)) '/T1*_Warped.nii.gz ' path_analysis num2str(IDs(id)) '/' tmp(1).name]);
    
%     eval(['!rm ' path_analysis num2str(IDs(id)) '/T1*Warped.nii.gz'])
%     eval(['!rm ' path_analysis num2str(IDs(id)) '/T1*.mat'])
%     eval(['!rm ' path_analysis num2str(IDs(id)) '/T1pt*.nii'])

end


%% segment

addpath(genpath('/Users/alex/Documents/MATLAB/biasCorrection/'))

for id=1:length(IDs)
    
%     try
    
    clear matlabbatch
    
%     matlabbatch{1}.spm.spatial.preproc.channel.vols = {[path_analysis num2str(IDs(id)) '/T1avg.nii,1']};
%     matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
%     matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
%     matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,1'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,2'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,3'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,4'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,5'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,6'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
%     matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
%     matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
%     matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
%     matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
%     matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
%     matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
%     matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
%         NaN NaN NaN];
    matlabbatch{1}.spm.util.imcalc.input = {
                                        [path_analysis num2str(IDs(id)) '/c1T1avg.nii,1']
                                        [path_analysis num2str(IDs(id)) '/c2T1avg.nii,1']
                                        [path_analysis num2str(IDs(id)) '/c3T1avg.nii,1']
                                        };
    matlabbatch{1}.spm.util.imcalc.output = 'brainmask';
    matlabbatch{1}.spm.util.imcalc.outdir = {[path_analysis num2str(IDs(id))]};
    matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)>0.7';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    spm_jobman('run',matlabbatch)
    
%     catch
%     clear matlabbatch
%     
%     matlabbatch{1}.spm.spatial.preproc.channel.vols = {[path_analysis num2str(IDs(id)) '/T1pt1_' num2str(sum(days(id,:))) '.nii,1']};
%     matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
%     matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
%     matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,1'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,2'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,3'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,4'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,5'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/Users/alex/Documents/MATLAB/spm12/tpm/TPM.nii,6'};
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
%     matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
%     matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
%     matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
%     matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
%     matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
%     matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
%     matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
%     matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
%         NaN NaN NaN];
%     matlabbatch{2}.spm.util.imcalc.input(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
%     matlabbatch{2}.spm.util.imcalc.input(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
%     matlabbatch{2}.spm.util.imcalc.input(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
%     matlabbatch{2}.spm.util.imcalc.output = 'brainmask';
%     matlabbatch{2}.spm.util.imcalc.outdir = {[path_analysis num2str(IDs(id))]};
%     matlabbatch{2}.spm.util.imcalc.expression = '(i1+i2+i3)>0.2';
%     matlabbatch{2}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%     matlabbatch{2}.spm.util.imcalc.options.dmtx = 0;
%     matlabbatch{2}.spm.util.imcalc.options.mask = 0;
%     matlabbatch{2}.spm.util.imcalc.options.interp = 1;
%     matlabbatch{2}.spm.util.imcalc.options.dtype = 4;
%     
%     spm_jobman('run',matlabbatch)    
%     end
%     
    
end


%% calculate intracranial volume

TIV=[];
for id=1:length(IDs)
    clear img
    img = spm_read_vols(spm_vol([path_analysis num2str(IDs(id)) '/brainmask.nii']));
    TIV(id,1)=IDs(id);
    TIV(id,2)=numel(find(img>0));
end


%% calculate brainstem volume

clc

brainstemVolumes = [];
brainstemVolumes = cell(length(IDs), 1);

for i = 1:length(IDs)
    
    fileData = fileread(['/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TIV/brainstem/' num2str(IDs(i)) '/brainstemSsVolumes.v13.txt']);
   
    lines = strsplit(fileData, '\n');
    tempStruct = struct();
    
    for j = 1:length(lines)
        if isempty(lines{j})
            continue;
        end
        
        % split line into brain area and volume
        splitLine = strsplit(lines{j});
        area = splitLine{1};
        volume = str2double(splitLine{2});
        
        tempStruct.(area) = volume;
    end
    
    brainstemVolumes{i} = tempStruct;
    
    copyfile(['/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TIV/brainstem/' num2str(IDs(i)) '/reg_brainstemSsLabels.v13.FSvoxelSpace.nii.gz'],...
        ['/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/TIV/brainstemSegs/' num2str(IDs(i)) '_BrainstemMask.nii.gz'])
end


%% threshold and calculate grey matter volume

GM_nat=[];
for id=1:length(IDs)
    clear img
    img = spm_read_vols(spm_vol([path_analysis num2str(IDs(id)) '/c1T1avg.nii']));
    GM_nat(id,1)=IDs(id);
    GM_nat(id,2)=numel(find(img>0.8));
end


%% normalise the GM volume

GMnorm(:,1)=TIV(:,1);
GMnorm(:,2)=(GM_nat(:,2)./TIV(:,2));
