function [ aap,resp ] = aamod_extrema(aap, task, subj)

% extract the whole-brain min and max voxel value in a nifti file, and, optionally,
% for a set of ROI specified either as TPM niftis or a list of coordinate/radius
% (if seed radius > 0, the seed ROI are extracted using marsbar, which must be
% installed).
%
% We also extract the mean and weighted mean because why not.
%

resp='';

switch task
	
    case 'report'
    case 'doit'
        			
		%---------------------------------------------------------------------------------
		% process (renamed) instream
		%---------------------------------------------------------------------------------
		
		% aside: we might also call this at the study level (e.g., second-level maps)
		      
		input_streamname = aap.tasklist.currenttask.inputstreams(1).stream{1};
		nifti_list = aas_getfiles_bystream(aap, subj, input_streamname);
	 				
		%---------------------------------------------------------------------------------
		% process
		%---------------------------------------------------------------------------------
					
		WHOLEBRAIN_max = [];
		weighted_WHOLEBRAIN_max = [];
		WHOLEBRAIN_min = [];
		weighted_WHOLEBRAIN_min = [];
		WHOLEBRAIN_mean = [];
		weighted_WHOLEBRAIN_mean = [];
		
		% whole brain (always computed)
		% "weighted" wholebrain value is just so nifti_max and weighted_nifti_max will have same length at end
			
		for index = 1:size(nifti_list,1)
			[ this_max, this_weighted_max, this_min, this_weighted_min, this_mean, this_weighted_mean ] = nifti_wholebrain_extrema(nifti_list(index,:), aap.tasklist.currenttask.settings.weight);
			WHOLEBRAIN_max = horzcat(WHOLEBRAIN_max, this_max);
			weighted_WHOLEBRAIN_max = horzcat(weighted_WHOLEBRAIN_max, this_weighted_max);
			WHOLEBRAIN_min = horzcat(WHOLEBRAIN_min, this_min);
			weighted_WHOLEBRAIN_min = horzcat(weighted_WHOLEBRAIN_min, this_weighted_min);
			WHOLEBRAIN_mean = horzcat(WHOLEBRAIN_mean, this_mean);
			weighted_WHOLEBRAIN_mean = horzcat(weighted_WHOLEBRAIN_mean, this_weighted_mean);
		end
		
		% process TPM (if defined)
		
		if (~isempty(aap.tasklist.currenttask.settings.TPMdir))
			
			TPMdir = aap.tasklist.currenttask.settings.TPMdir;
			
			% make any errors fatal so as to not potentially silently generate wrong results
			
			if (~exist(TPMdir,'dir'))
				aap = aas_log(aap,true,sprintf('%s: TPMdir (%s) does not exist.', mfilename, TPMdir));
			else
				% need full paths for TPM fies
				TPM_list = cellstr(spm_select('FPList', TPMdir, '^*.nii'));
				if (isempty(TPM_list)) 
					aap = aas_log(aap,true,sprintf('%s: No .nii files in TPMdir (%s).', mfilename, TPMdir));
				end
				% also need the filenames for ROI naming
				TPM_names = cellstr(spm_select('List', TPMdir, '^*.nii'));
			end
				
			TPM_max = [];
			weighted_TPM_max = [];
			TPM_min = [];
			weighted_TPM_min = [];
			TPM_mean = [];
			weighted_TPM_mean = [];
			
			for index = 1:size(nifti_list,1)
				[ this_max, this_weighted_max, this_min, this_weighted_min, this_mean, this_weighted_mean ] = nifti_extrema_in_TPM(nifti_list(index,:), TPM_list);
				% append results from this pass (note results returned as rowvecs)
				TPM_max = horzcat(TPM_max, this_max);
				weighted_TPM_max = horzcat(weighted_TPM_max, this_weighted_max);
				TPM_min = horzcat(TPM_min, this_min);
				weighted_TPM_min = horzcat(weighted_TPM_min, this_weighted_min);
				TPM_mean = horzcat(TPM_mean, this_mean);
				weighted_TPM_mean = horzcat(weighted_TPM_mean, this_weighted_mean);
			end
			
		end
		
		% process seed (if defined)
		
		if (~isempty(aap.tasklist.currenttask.settings.SEEDS))
		
			SEED_list = aap.tasklist.currenttask.settings.SEEDS;
			SEED_names = extractfield(SEED_list.seed,'name');
			SEED_names = SEED_names';
							
			SEED_max = [];
			weighted_SEED_max = [];
			SEED_min = [];
			weighted_SEED_min = [];
			SEED_mean = [];
			weighted_SEED_mean = [];
			
			for index = 1:size(nifti_list,1)
				[ this_max, this_weighted_max, this_min, this_weighted_min, this_mean, this_weighted_mean ] = nifti_extrema_in_seed_ROI(aap, nifti_list(index,:), SEED_list.seed);			
				SEED_max = horzcat(SEED_max, this_max');
				weighted_SEED_max = horzcat(weighted_SEED_max, this_weighted_max');	
				SEED_min = horzcat(SEED_min, this_min');
				weighted_SEED_min = horzcat(weighted_SEED_min, this_weighted_min');	
				SEED_mean = horzcat(SEED_mean, this_mean');
				weighted_SEED_mean = horzcat(weighted_SEED_mean, this_weighted_mean');	
			end
			
		end
		
		%---------------------------------------------------------------------------------
		% Desc outputs
		%---------------------------------------------------------------------------------
		
		% combine whole-brain/TPM/seed results for output
		% results should go [ #ROI x #files-processed (e.g. contrasts) ]
		% note: vertcat cool with empties
		
		nifti_max = vertcat(WHOLEBRAIN_max, TPM_max, SEED_max);
		weighted_nifti_max = vertcat(weighted_WHOLEBRAIN_max, weighted_TPM_max, weighted_SEED_max);
		nifti_min = vertcat(WHOLEBRAIN_min, TPM_min, SEED_min);
		weighted_nifti_min = vertcat(weighted_WHOLEBRAIN_min, weighted_TPM_min, weighted_SEED_min);
		nifti_mean = vertcat(WHOLEBRAIN_mean, TPM_mean, SEED_mean);
		weighted_nifti_mean = vertcat(weighted_WHOLEBRAIN_mean, weighted_TPM_mean, weighted_SEED_mean);

		ROI_names = vertcat('wholebrain', TPM_names, SEED_names);
	
		savedir = aas_getsubjpath(aap,subj);
		
		fname=fullfile(savedir, 'nifti_max.mat');
 		save(fname, 'nifti_max');
		aap=aas_desc_outputs(aap, subj, 'nifti_max', fname);
		
		fname=fullfile(savedir, 'weighted_nifti_max.mat');
		save(fname, 'weighted_nifti_max');
		aap=aas_desc_outputs(aap, subj, 'weighted_nifti_max', fname);
		
		fname=fullfile(savedir, 'nifti_min.mat');
 		save(fname, 'nifti_min');
		aap=aas_desc_outputs(aap, subj, 'nifti_min', fname);
		
		fname=fullfile(savedir, 'weighted_nifti_min.mat');
		save(fname, 'weighted_nifti_min');
		aap=aas_desc_outputs(aap, subj, 'weighted_nifti_min', fname);
		
		fname=fullfile(savedir, 'nifti_mean.mat');
 		save(fname, 'nifti_mean');
		aap=aas_desc_outputs(aap, subj, 'nifti_mean', fname);
		
		fname=fullfile(savedir, 'weighted_nifti_mean.mat');
		save(fname, 'weighted_nifti_mean');
		aap=aas_desc_outputs(aap, subj, 'weighted_nifti_mean', fname);
		
		fname=fullfile(savedir, 'ROI_names.mat');
		save(fname, 'ROI_names');
		aap=aas_desc_outputs(aap, subj, 'ROI_names', fname);

 			
	case 'checkrequirements'
		% add marsbar check here
		
    otherwise
        aas_log(aap,1,sprintf('%s: Unknown task %s',mfilename, task));
		
	end
	
end


%---------------------------------------------------------------------------------------------------------------------
function [ vmax, weighted_vmax, vmin, weighted_vmin, vmean, weighted_vmean ] = nifti_wholebrain_extrema(fname, weight)
%---------------------------------------------------------------------------------------------------------------------

handle = spm_vol(fname);
data = spm_read_vols(handle);
data(isnan(data)) = 0;

vmax = max(max(max(data)));
weighted_vmax = weight * vmax;

vmin = min(min(min(data)));
weighted_vmin = weight * vmin;

vmean = mean(mean(mean(data)));
weighted_vmean = weight * vmean;

end


%-------------------------------------------------------------------------------------------------------------------
function [ vmax, weighted_vmax, vmin, weighted_vmin, vmean, weighted_vmean ] = nifti_extrema_in_TPM(fname, TPMlist)
%-------------------------------------------------------------------------------------------------------------------

if (isempty(TPMlist))
	vmax = [];
	weighted_vmax = [];
	return;
end

vmax = zeros(numel(TPMlist),1);
weighted_vmax = zeros(numel(TPMlist),1);

vmin = zeros(numel(TPMlist),1);
weighted_vmin = zeros(numel(TPMlist),1);

vmean = zeros(numel(TPMlist),1);
weighted_vmean = zeros(numel(TPMlist),1);

% reslice data to match ROI templates

resliceOpts = [];
resliceOpts.mask = false;		% no masking
resliceOpts.mean = false;		% don't write a mean image
resliceOpts.interp = 1;			% default interp
resliceOpts.which = 1;			% don't reslice the first image
resliceOpts.wrap = [1 1 0];		% fMRI wrap around
resliceOpts.prefix = 'r';		% filename really doesn't matter - we delete it

spm_reslice({TPMlist{1} fname}, resliceOpts);

[p,n,e] = fileparts(fname);
tempfile = fullfile(p,['r' n e]);
handle = spm_vol(tempfile);
data = spm_read_vols(handle);
data(isnan(data)) = 0;

for index = 1:numel(TPMlist)
	handle = spm_vol(TPMlist{index});
	roi = spm_read_vols(handle);
	vmax(index) = max(max(max((~~roi).*data)));
	weighted_vmax(index) = max(max(max(roi.*data)));
	vmin(index) = min(min(min((~~roi).*data)));
	weighted_vmin(index) = min(min(min(roi.*data)));
	vmean(index) = mean(mean(mean((~~roi).*data)));
	weighted_vmean(index) = mean(mean(mean(roi.*data)));
end

delete(tempfile);

end



%---------------------------------------------------------------------------------------------------------------------------
function [ vmax, weighted_vmax, vmin, weighted_vmin, vmean, weighted_vmean ] = nifti_extrema_in_seed_ROI(aap, fname, SEEDS)
%---------------------------------------------------------------------------------------------------------------------------

handle = spm_vol(fname);
imgMat = handle.mat;
imgDim = handle.dim;
  
vmax = [];
weighted_vmax = [];

vmin = [];
weighted_vmin = [];

vmean = [];
weighted_vmean = [];

for index = 1:numel(SEEDS)

	seed_name = SEEDS(index).name;
	seed_center = SEEDS(index).center;
	radius = SEEDS(index).radius;
	weight = SEEDS(index).weight;
	
	if (isempty(weight)); weight=1.0; end;

	% marsbar ROI extraction
	%
	% "seed" is one voxel; "ROI" is all voxels within radius
	%
	% maroi_sphere doesn't like radius = 0 (voxpts returns empty ROI) so
	% for 1-voxel seed (i.e., radius=0) just do mat\[seedmm 1]'
	% (note maroi_sphere sometimes fails if radius = 1)

	if (radius> 0) 

		ROI = maroi_sphere(struct('centre', seed_center, 'radius', radius));
		roiIJK = voxpts(ROI, struct('mat', imgMat, 'dim', imgDim));
		if (isempty(roiIJK))
			aas_log(aap, false, sprintf('%s: Seed %s resulted in empty ROI. Skipping...', mfilename, seed_name));
			voxels = 0;
		else
			roiIJK = [ roiIJK; ones(1, size(roiIJK,2)) ];
			voxels = spm_get_data(handle, roiIJK);
		end

	else
		
		roiIJK =  round(imgMat\[seed_center 1]');
		voxels = spm_get_data(handle, roiIJK);

	end

	vmax = [ vmax max(max(max(voxels))) ];
	weighted_vmax = [ weighted_vmax weight*max(max(max(voxels)))];
	
	vmin = [ vmin min(min(min(voxels))) ];
	weighted_vmin = [ weighted_vmin weight*min(min(min(voxels)))];
	
	vmean = [ vmean mean(mean(mean(voxels))) ];
	weighted_vmean = [ weighted_vmean weight*mean(mean(mean(voxels)))];

end

end

