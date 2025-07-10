%% extract LC contrast ratio

clc;clear

setenv('PATH', [getenv('PATH') ':/Applications/freesurfer/mni/bin:/Users/alex/ants/bin']);
setenv('ANTSPATH','/Users/alex/ants/bin')

% load('/Users/alex/Documents/LCaverage/IDs.mat')
ID={'4001','4002','4003','4004','4005','4006','4007','4008','4009','4010','4011','4012','4013','4014',...
    '4015','4016','4017','4018','4019','4020','4021','4022','4023','4024','4025','4026','4027','4028','4029',...
    '4030','4031','4032','4033'}

path_parent='/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/LCaverage_MRPET/';
% pons='/Users/yeojin/Desktop/E_data/EE_atlases_templates/mni_icbm152_ponsmask_25vox.nii';
pons_native = '/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/LCaverage_pons_and_templatemask/';

%% registration T1->LCslab averaged

for id=1:length(ID)
    
    clear FixedImage MovingImage OutputImage
    MovingImage     = [path_parent ID{id} '/T1WB_corrected.nii'];
    FixedImage      = [path_parent ID{id} '/LCslab_averaged.nii'];
    OutputPrefix    = [path_parent ID{id} '/T1toLCslab_'];
    
    % register
%      eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -f ' FixedImage ' -m ' MovingImage ' -o ' OutputPrefix])
    % transform
    eval(['!antsApplyTransforms -d 3 -v 0 -n NearestNeighbor  -t ' path_parent ID{id}...
        '/T1toLCslab_0GenericAffine.mat -t [' path_parent ID{id} '/NLreg_T1WB_to_template_0GenericAffine.mat,1] -t ' ...
        path_parent 'NLreg_template_to_MNI_1InverseWarp.nii.gz -t [' path_parent '/NLreg_template_to_MNI_0GenericAffine.mat,1] -i '...
        pons ' -r ' FixedImage ' -o ' path_parent ID{id} '/pons_native.nii'])
    
end

%% LC template mask to native

LCtemplate = '/Users/alex/Qumulo/IKND/Alex/Masks/Masks/mni_icbm152/mni_icbm152_LCmetaMask_final.nii';

for id=1%:length(ID)

    %     try
    %         copyfile([path_ponsmask ID{id} '/pons_native.nii'],...
    %             [path_LCslab ID{id} '/pons_native.nii'])
    %     catch
    clear FixedImage MovingImage OutputImage
    MovingImage     = ['/Volumes/korokdorf/MRPET/coreg_mri/bothsessions/' ID{id} '/data/T1WB.nii'];
    FixedImage      = [path_parent ID{id} '/LCslab_averaged_noNorm.nii'];
    OutputPrefix    = [pons_native ID{id} '/T1toLCslabNoNorm_'];

    % register
    eval(['!antsRegistrationSyNQuick.sh -d 3 -t r -f ' FixedImage ' -m ' MovingImage ' -o ' OutputPrefix])
    % transform
    eval(['!antsApplyTransforms -d 3 -v 0 -n NearestNeighbor -t ' OutputPrefix '0GenericAffine.mat -t ' pons_native ID{id} '/NLreg_T1WB_to_template_1InverseWarp.nii.gz'...
        ' -t [' pons_native ID{id} '/NLreg_T1WB_to_template_0GenericAffine.mat,1] -t ' ...
        pons_native 'NLreg_template_to_MNI_1InverseWarp.nii.gz -t ['...
        pons_native 'NLreg_template_to_MNI_0GenericAffine.mat,1] -i '...
        LCtemplate ' -r ' FixedImage ' -o ' pons_native ID{id} '/LCtemplate_on_native.nii'])
    %     end

end

%% extract LC CR (in native?)

for id=1:length(ID)
    
    clear LC S_LC_mean S_pons_mean S_LC_max S_pons_max
    
    % load images
    LC          = spm_read_vols(spm_vol([path_parent ID{id} '/LCmask_YY.nii']));
    pons_nat    = spm_read_vols(spm_vol([pons_native ID{id} '/pons_native.nii']));
    LCslab      = spm_read_vols(spm_vol([path_parent ID{id} '/LCslab_averaged_noNorm.nii']));
    
    % define masks
    [x_lc,y_lc,z_lc]= ind2sub(size(LC),find(LC~=0));
    [x_p,y_p,z_p]   = ind2sub(size(pons_nat),find(pons_nat~=0));
    coord_LC        = [x_lc,y_lc,z_lc];
    coord_pons      = [x_p,y_p,z_p];
    slices          = unique(z_lc); % identify slices in LC masks
    
    % extract voxel values from each mask, in each slice
    for sl=1:length(slices)
        
        clear tmp_lc voxvals_lc tmp_p voxvals_p
        
        % LC
        tmp_lc  = coord_LC(z_lc==slices(sl),:);
        tmp_p   = coord_pons(z_p==slices(sl),:);
        for vox=1:size(tmp_lc,1)
            voxvals_lc(vox,1)=LCslab(tmp_lc(vox,1),tmp_lc(vox,2),tmp_lc(vox,3));
        end
        S_LC_mean(sl,1) = mean(voxvals_lc(:));
        S_LC_max(sl,1)  = max(voxvals_lc(:));
        
        % pons
        for vox=1:size(tmp_p,1)
            voxvals_pons(vox,1)=LCslab(tmp_p(vox,1),tmp_p(vox,2),tmp_p(vox,3));
        end
        S_pons_mean(sl,1) = mean(voxvals_pons(:));
        S_pons_max(sl,1)  = max(voxvals_pons(:));
        
    end
    
    LC_CR_mean(id,1) = mean((S_LC_mean - S_pons_mean)./ S_pons_mean);
    LC_CR_max(id,1)  = max((S_LC_max - S_pons_max)./ S_pons_max);
    
end


%% extract LC CR (in LC template space)

LC_template          = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EE_atlases_templates/KerenLCmask_on_LCslabtemplate.nii']));
pons_template        = spm_read_vols(spm_vol(['/Users/yeojin/Desktop/E_data/EE_atlases_templates/pons25vox_on_LCslabtemplate.nii']));


% define masks
[x_lc,y_lc,z_lc]= ind2sub(size(LC_template),find(LC_template~=0));
[x_p,y_p,z_p]   = ind2sub(size(pons_template),find(pons_template~=0));
coord_LC        = [x_lc,y_lc,z_lc];
coord_pons      = [x_p,y_p,z_p];
slices          = unique(z_lc); % identify slices in LC masks


for id=1:length(ID)
        
    % load images

    LCslab      = spm_read_vols(spm_vol([path_parent 'LCslab_on_Template_' ID{id} '.nii']));

    % extract voxel values from each mask, in each slice
    for sl=1:length(slices)

        clear tmp_lc voxvals_lc tmp_p voxvals_p

        % LC
        tmp_lc  = coord_LC(z_lc==slices(sl),:);
        tmp_p   = coord_pons(z_p==slices(sl),:);
        for vox=1:size(tmp_lc,1)
            voxvals_lc(vox,1)=LCslab(tmp_lc(vox,1),tmp_lc(vox,2),tmp_lc(vox,3));
        end
        S_LC_mean(sl,1) = mean(voxvals_lc(:));
        S_LC_max(sl,1)  = max(voxvals_lc(:));

        % pons
        for vox=1:size(tmp_p,1)
            voxvals_pons(vox,1)=LCslab(tmp_p(vox,1),tmp_p(vox,2),tmp_p(vox,3));
        end
        S_pons_mean(sl,1) = mean(voxvals_pons(:));
        S_pons_max(sl,1)  = max(voxvals_pons(:));

    end

    CR_mean(id,1) = mean((S_LC_mean - S_pons_mean)./ S_pons_mean);
    CR_max(id,1)  = max((S_LC_max - S_pons_max)./ S_pons_max);
    
end

%% LC volume from the manual segmentation

load('/Users/alex/Qumulo/IKND/Alex/paperwriting/MRPET/data/MRPET_TIVandGMnorm.mat')

for id=1:length(ID)
    
    clear LC 
    
    % load images
    LC          = spm_read_vols(spm_vol([path_parent ID{id} '/LCmask_YY.nii']));
    
    LC_vol_orig(id,1) = sum(LC(:));
    LC_vol_TIVnorm(id,1)  = LC_vol_orig(id,1)/TIV(id,2);
    
end

