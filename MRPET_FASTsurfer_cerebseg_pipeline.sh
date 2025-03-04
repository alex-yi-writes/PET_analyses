#!/bin/bash

for ID in 40021 #40041 40042 40051 40052 40062 40071 40081 40082 40091 40092 40101 40102 40111 40112 40121 40122 40131 40132 40141 40142 40151 40152 40161 40162 40171 40172 40181 40182 40191 40192 40201 40202 40211 40212 40221 40222 40231 40232 40241 40242 40251 40252 40261 40262 40271 40281 40282 40292 40301 40302 40311 40312 40321 40322 40331 40332 
do

	mkdir /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1
	mkdir /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2

	cd /Users/alex/FastSurfer
	./run_fastsurfer.sh --sid ${ID}1 --sd /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1 --t1 /Volumes/korokdorf/MRPET/coreg_roi/${ID}/data/T1pt1.nii --vox_size 1 --seg_only --no_hypothal 
	mri_convert -it mgz -ot nii /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/${ID}1/mri/orig.mgz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/${ID}1/mri/orig.nii.gz
	./run_fastsurfer.sh --sid ${ID}2 --sd /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2 --t1 /Volumes/korokdorf/MRPET/coreg_roi/${ID}/data/T1pt2.nii --vox_size 1 --seg_only --no_hypothal 
	mri_convert -it mgz -ot nii /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/${ID}2/mri/orig.mgz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/${ID}2/mri/orig.nii.gz

	antsRegistrationSyNQuick.sh -d 3 -f /Volumes/korokdorf/MRPET/coreg_roi/${ID}/data/T1pt1.nii -m /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/${ID}1/mri/orig.nii.gz -t r -o /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/reg_FStoNAT_
	antsApplyTransforms -d 3 -i /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/${ID}1/mri/cerebellum.CerebNet.nii.gz -r /Volumes/korokdorf/MRPET/coreg_roi/${ID}/data/T1pt1.nii -o /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cerebseg_on_pt1.nii.gz -n NearestNeighbor -t /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/reg_FStoNAT_0GenericAffine.mat

	antsRegistrationSyNQuick.sh -d 3 -f /Volumes/korokdorf/MRPET/coreg_roi/${ID}/data/T1pt2.nii -m /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/${ID}2/mri/orig.nii.gz -t r -o /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/reg_FStoNAT_
	antsApplyTransforms -d 3 -i /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/${ID}2/mri/cerebellum.CerebNet.nii.gz -r /Volumes/korokdorf/MRPET/coreg_roi/${ID}/data/T1pt2.nii -o /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cerebseg_on_pt2.nii.gz -n NearestNeighbor -t /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/reg_FStoNAT_0GenericAffine.mat



	##################################################
	######## segmentation part 1, inflow segs ########
	##################################################

	# init an empty image (using fslmaths; create a 0 image from cerebseg_on_pt1)
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cerebseg_on_pt1.nii.gz -mul 0 /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz

	# left cerebellum grey: map subregion labels to label 8 of recon-all label
	for lab in 601 603 605 608 611 614 617 620 623 626; do
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cerebseg_on_pt1.nii.gz -thr $lab -uthr $lab -bin -mul 8 tmp.nii.gz
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz
	done

	# left cerebellum white: label 7
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cerebseg_on_pt1.nii.gz -thr 7 -uthr 7 -bin -mul 7 tmp.nii.gz
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz

	# right cerebellum grey: map subregion labels to label 47 of recon-all label
	for lab in 602 604 607 610 613 616 619 622 625 628; do
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cerebseg_on_pt1.nii.gz -thr $lab -uthr $lab -bin -mul 47 tmp.nii.gz
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz
	done

	# right cerebellum white: label 46
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cerebseg_on_pt1.nii.gz -thr 46 -uthr 46 -bin -mul 46 tmp.nii.gz
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz

	# remove any remaining unwanted labels by thresholding; set labels not equal to 7,8,46,47 to 0
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz -uthr 47 -thr 7 /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz

	rm tmp.nii.gz

	echo "new cerebellar segmentation created"

	# create a binary mask for all cerebellar voxels in the original parcellation
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -uthr 8 -thr 8 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/left_grey_mask_pt1.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -uthr 7 -thr 7 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/left_white_mask_pt1.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -uthr 47 -thr 47 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/right_grey_mask_pt1.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -uthr 46 -thr 46 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/right_white_mask_pt1.nii.gz

	# combine masks
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/left_grey_mask_pt1.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/left_white_mask_pt1.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/right_grey_mask_pt1.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/right_white_mask_pt1.nii.gz -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cereb_mask_pt1.nii.gz

	# replace those voxels with 0:
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cereb_mask_pt1.nii.gz -mul -1 -add 1 /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/inverted_cereb_mask_pt1.nii.gz && fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -mul /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/inverted_cereb_mask_pt1.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/aparc_no_cereb_pt1.nii.gz

	echo "cerebellar labels removed from whole-brain parcellation"

	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -thr 7 -uthr 8 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/left_cereb_mask_pt1.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt1_nat_labelled.nii -thr 46 -uthr 47 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/right_cereb_mask_pt1.nii.gz

	# combine the left and right cerebellum masks into one
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/left_cereb_mask_pt1.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/right_cereb_mask_pt1.nii.gz -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cereb_mask_pt1.nii.gz

	# multiply your new segmentation with the cerebellum mask to make sure that only voxels inside the cerebellum remain
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_on_pt1.nii.gz -mul /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/cereb_mask_pt1.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_pt1_labelled.nii.gz

	# now, when merging, only the cerebellum area gets labelled anew
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/aparc_no_cereb_pt1.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/new_cerebseg_pt1_labelled.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}1/aparc+aseg_pt1_nat_labelled_upated.nii.gz




	##################################################
	####### segmentation part 2, bsl+task segs #######
	##################################################

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cerebseg_on_pt2.nii.gz -mul 0 /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz

	for lab in 601 603 605 608 611 614 617 620 623 626; do
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cerebseg_on_pt2.nii.gz -thr $lab -uthr $lab -bin -mul 8 tmp.nii.gz
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz
	done

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cerebseg_on_pt2.nii.gz -thr 7 -uthr 7 -bin -mul 7 tmp.nii.gz
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz

	for lab in 602 604 607 610 613 616 619 622 625 628; do
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cerebseg_on_pt2.nii.gz -thr $lab -uthr $lab -bin -mul 47 tmp.nii.gz
	  fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz
	done

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cerebseg_on_pt2.nii.gz -thr 46 -uthr 46 -bin -mul 46 tmp.nii.gz
	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz -add tmp.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz


	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz -uthr 47 -thr 7 /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz

	rm tmp.nii.gz

	echo "new cerebellar segmentation created"

	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -uthr 8 -thr 8 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/left_grey_mask_pt2.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -uthr 7 -thr 7 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/left_white_mask_pt2.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -uthr 47 -thr 47 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/right_grey_mask_pt2.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -uthr 46 -thr 46 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/right_white_mask_pt2.nii.gz

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/left_grey_mask_pt2.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/left_white_mask_pt2.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/right_grey_mask_pt2.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/right_white_mask_pt2.nii.gz -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cereb_mask_pt2.nii.gz

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cereb_mask_pt2.nii.gz -mul -1 -add 1 /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/inverted_cereb_mask_pt2.nii.gz && fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -mul /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/inverted_cereb_mask_pt2.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/aparc_no_cereb_pt2.nii.gz

	echo "cerebellar labels removed from whole-brain parcellation"

	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -thr 7 -uthr 8 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/left_cereb_mask_pt2.nii.gz
	fslmaths /Volumes/korokdorf/MRPET/coreg_roi/${ID}/aparc+aseg_pt2_nat_labelled.nii -thr 46 -uthr 47 -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/right_cereb_mask_pt2.nii.gz

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/left_cereb_mask_pt2.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/right_cereb_mask_pt2.nii.gz -bin /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cereb_mask_pt2.nii.gz

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_on_pt2.nii.gz -mul /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/cereb_mask_pt2.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_pt2_labelled.nii.gz

	fslmaths /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/aparc_no_cereb_pt2.nii.gz -add /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/new_cerebseg_pt2_labelled.nii.gz /Users/alex/Dropbox/paperwriting/MRPET/data/newseg/${ID}2/aparc+aseg_pt2_nat_labelled_upated.nii.gz

	echo "all done"
done
