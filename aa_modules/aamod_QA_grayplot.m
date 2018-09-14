function [aap,resp] = aamod_QA_grayplot(aap,task,subj,sess)
%
% generate epi grayplot and QA metrics DVARS, FD, and slice variance
% 
% Reference:
%
% 	Jonathan D. Power
% 	A simple but useful way to assess fMRI scan qualities
% 	NeuroImage 154, 150-158, (2017) 
% 	https://doi.org/10.1016/j.neuroimage.2016.08.009
%

resp='';

switch task
	
    case 'domain'
		
        resp = 'session';
		
    case 'description'
		
        resp = 'generate epi grayplot and QA metrics';

    case 'report'
		
		sesspath = aas_getsesspath(aap, subj, sess);
		grayplot = fullfile(sesspath,'QA_grayplot.jpg');
		aap = aas_report_add(aap, subj, '<table><tr><td>');
		aap = aas_report_addimage(aap, subj, grayplot);
		aap = aas_report_add(aap, subj, '</td></tr></table>');
    
    case 'doit'

		fname.bold = aas_getfiles_bystream(aap, subj, sess, 'epi');
		fname.epimask = aas_getfiles_bystream(aap, subj, 'native_brainmask');
				
		for ind=1:size(fname.bold,1)
			bold_headers = spm_vol(sprintf('%s',fname.bold(ind,:)));
			img.bold = spm_read_vols(bold_headers);
		end
			
		% sanity check -- the mask must match the epi. It won't if the
		% user didn't bother to reslice it in the tasklist (masks come
		% from segmentation and segmentation comes from the structural)
				
		mask_header = spm_vol(fname.epimask);
		representative_epi_header = bold_headers(1);
		
		if ( ~isequal(mask_header.dim, representative_epi_header.dim) || ...
				norm(mask_header.mat-representative_epi_header.mat)>0.01 )

			% fix or bail?
			
			if (aap.tasklist.currenttask.settings.resliceGMMaskIfNecessary)

				aas_log(aap, false, sprintf('\n%s: Reslicing mask to match epi...\n', mfilename));
				
				resliceOpts = [];
				resliceOpts.mask = false;
				resliceOpts.mean = false;
				resliceOpts.interp = 0;		% NN (because it's a MASK)
				resliceOpts.which = 1;		% don't reslice the first image
				resliceOpts.wrap = [0 0 0];
				resliceOpts.prefix = 'r';
				
				spm_reslice({representative_epi_header.fname mask_header.fname}, resliceOpts);
				
				[ p,n,e ] = fileparts(fname.epimask);
				fname.epimask = fullfile(p,[resliceOpts.prefix n e]);
				
			else
				aas_log(aap, true, sprintf('\nMask incompatible with epi. Reslice in tasklist prior to using %s\n', mfilename));
			end
			
		end
		
		img.epimask = spm_read_vols(spm_vol(fname.epimask));

		% ----------------------------------------------------------------
		% calculate SV (slice variance)
		% ----------------------------------------------------------------
		
		tsfn = aas_getfiles_bystream(aap, subj, sess, 'tsdiffana');
% 		load(tsfn);

		% this is for current tsdiffana

		tsdiffana = load(tsfn);
		slicediff = tsdiffana.qa.slice.diff;
		globals = tsdiffana.qa.global.mean;
		td = tsdiffana.qa.global.diff;

		imgno = size(slicediff,1)+1;
		mom = mean(globals);
        
		sv = td/mom;
		SV = [ 0; sv ];

		% ----------------------------------------------------------------
		% calculate DVARS
		% ----------------------------------------------------------------
			
		dd = size(img.bold);
		im = reshape(img.bold,[dd(1)*dd(2)*dd(3)],dd(4));

		% detrend across volumes (note transpose) and extract GM voxels
		% ~~ is a clever trick to make sure mask is 0/1 (not 0/something)
					
		detrend_img = detrend(im')'; 
		mask_img = detrend_img(~~img.epimask(:),:);
		
		% we often get a transient at the edges of the epi
		% define a tukey window to smooth them off (we have
		% to roll our own because Matlab keeps finding the
		% broken version of tukey in the fieldtrip toolbox).
				
		n = size(mask_img,2);
		r = 0.1; % 10%
		t = linspace(0,1,n)';
		per = r/2; 
		tl = floor(per*(n-1))+1;
		th = n-tl+1;
		w = [((1+cos(pi/per*(t(1:tl)-per)))/2); ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end)-1+per)))/2)];
		% mask_img = mask_img .* w' doesn't work on older versions of matlab - unroll it
		for index=1:size(mask_img,1)
			mask_img(index,:) = mask_img(index,:) .* w';
		end
					
		dv = diff(mask_img,1,2);
		[~,dimc] = size(dv);
		dv = sqrt(mean(dv(:,1:dimc).^2));
		DVARS = [0;dv'];
		
		% ----------------------------------------------------------------
		% calculuate FD
		% ----------------------------------------------------------------
		
		M = spm_load(aas_getfiles_bystream(aap, subj, sess, 'realignment_parameter'));
		v = [1 1 1 50 50 50]; 
		V = diag(v);
		M = M*V; % convert to equivalent 50mm radius circle
		dM = diff(M);
		[~,colstr] = size(dM);
		d0 = zeros(1,colstr);
		dM = [d0;dM];
		dM = abs(dM);
		mat = ones(colstr,1);
		FD = dM*mat;

		% plot
		
		variables = [ globals SV DVARS FD ];
		var_names = {'Global mean' 'Scaled Var' 'DVARS' 'FD'};
		threshold_SV = nanmean(SV)+3*mad(SV);
		threshold_DVARS = nanmean(dv)+3*mad(dv);
		threshold_FD = nanmean(FD)+3*mad(FD);
		dthres = [threshold_SV threshold_DVARS threshold_FD];
		
		% the movegui trick here centers the window, which looks cool
		% but it breaks plotting if running headless (i.e. cluster)
		
		if (strcmp(aap.options.wheretoprocess, 'localsingle'))
			h = figure('Position',[0 0 800 960],'Visible', 'off');
			movegui(h,'center');
			set(h,'Visible', 'on');
		else
			h = figure('Position',[0 0 800 960]);
		end
				
		for i=1:4
			
			m=3*(i-1)+1;n=3*(i-1)+2;
			subplot(6,3,[m,n]);
			
			if i==1
				% global mean plot is odd duck
				plot(1:dd(4),variables(:,i));
				axis([0.5 imgno+0.5 -Inf Inf]);
			else 
				plot(2:dd(4),variables(2:end,i));
				axis([0.5 imgno+0.5 0 Inf]);
				hline = refline([0 dthres(i-1)]);
				hline.Color='r';
			end
				
			ylabel(var_names{i});
			ax=gca;p=ax.Position;
			pvector=[p(1)+p(3) p(2) 0.2 p(4)];
			[counts,centers]=hist(variables(:,i),50);
			subplot('position', pvector);
			barh(centers,counts,'FaceColor',[0.5 0.5 0.5]);
			axis([0 Inf -Inf Inf]);
			axis off;
			
		end
		
		subplot(6,3,[1,2]);
		title(strrep([ aas_getsubjname(aap,subj) '/' aas_getsessname(aap,sess) ],'_','-'));

		% create and save diagnostic gray plot
		
		points_to_plot = aap.tasklist.currenttask.settings.numberOfTimePointsToPlot;

		[row,~] = size(mask_img);		
		sp = round(linspace(1,row, points_to_plot));
		smooth_img = imfilter(mask_img(sp,:), fspecial('Gaussian'));
		
		subplot(6,3,[13 14 16 17]);
		axis([0.5 imgno+0.5 -Inf Inf]);
		imagesc(smooth_img);
		xlabel('volume');ylabel('voxel');
		colormap(gray);
		ax=gca;p=ax.Position;
		pvectors=[p(1)+p(3)+0.02 p(2) 0.02 p(4)];
		hb=colorbar;
		set(hb,'Units', 'normalized', 'position', pvectors);
		outname = fullfile(fileparts(fname.bold),'QA_grayplot.jpg');		
		set(h,'Renderer','opengl');
		set(findall(h,'Type','text'),'FontUnits','normalized');
		print(h, '-djpeg', '-r150', outname);
		close(h);
							
		% desc the outputs
		
		outfile = fullfile(fileparts(fname.bold), 'QA.mat'); 
		save(outfile,'SV','DVARS','FD');
		aap = aas_desc_outputs(aap, subj, sess, 'QA', outfile);
        
		close all;

end