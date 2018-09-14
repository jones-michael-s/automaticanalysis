function [aap,resp] = aamod_scrub_epis(aap, task, subj, sess)
%
% aamod_scrub_epis
%
%	scrub epi volumes based on QA metrics and/or explicit criteria
%
%	NB: this module generates a list of delta regressors to "scrub" data
%	from the GLM -- if you want to literally delete the frames, use
%	aamod_delete_epi_frames instead
% 
% 	the scrub list is constructed as the union of an input scrublist,
%	motion frames derived from a QA metric, and frames specified using 
%	an explicit criteria (see "additional_criteria")
%
% Revision History
%
% spring 2018 -- new [MSJ]
%

resp='';

switch task
	
    case 'report'
        
    case 'doit'
	
		% initialization stuff
		
		input_scrublist = [ ];
		motion_scrublist = [ ];
		additional_scrublist = [ ];
		master_scrublist = [ ];
		
		session_path = aas_getsesspath(aap, subj, sess);
		
		% get the number of volumes (frames) in the epi
		
		epi_file_list = aas_getimages_bystream(aap, subj, sess, 'epi');
		numvol = size(epi_file_list, 1);
		
		% if 4D epi (i.e. one file) we need to extract numvol from header
		
		if (numvol == 1)
			numvol = numel(spm_vol(epi_file_list));
		end
		
		% 1) begin with input scrublist if the stream exists

		if (aas_stream_has_contents(aap, 'scrublist'))
			temp = aas_getfiles_bystream(aap, subj, sess, 'scrublist');
			input_scrublist = load(temp);
		end
		

		% 2) add a motion-based scrublist using a QA metric if the stream exists

        if aas_stream_has_contents(aap,'QA')
			
			sanity_check = true;
			
			if (isempty(aap.tasklist.currenttask.settings.QA_metric))
				aas_log(aap, false, 'Found QA stream but no QA metric (DVARS|FD|SV) selected. Ignoring QA stream...\n');
				sanity_check = false;
			end
			
			if (isempty(aap.tasklist.currenttask.settings.QA_threshold))
				aas_log(aap, false, 'Found QA stream but no QA threshold defined. Ignoring QA stream...\n');
				sanity_check = false;
			end
			
			if sanity_check

				temp = aas_getfiles_bystream(aap, subj, sess, 'QA');
				GP_Parameter = load(temp);

				QAfieldname = aap.tasklist.currenttask.settings.QA_metric;
				dat = getfield(GP_Parameter, QAfieldname);

				threshold = aap.tasklist.currenttask.settings.QA_threshold;
				motion_scrublist = find(dat>threshold);

			end
		
		end
     
		% 3) add frames from additional criteria parameter if one is defined
	
		if (~isempty(aap.tasklist.currenttask.settings.additional_criteria))
			additional_criteria = aap.tasklist.currenttask.settings.additional_criteria;
			additional_scrublist = eval(additional_criteria);
		end
		
		% matlab "union" will sort and remove duplicates
		
		master_scrublist = union(input_scrublist, motion_scrublist);
		master_scrublist = union(master_scrublist, additional_scrublist);
		
		% we may extend scrubbing a fixed number of frames before and
		% after any identified frames targeted for scrubbing - the
		% range is specified by the prekillzone and postkillzone 
		% parameters (0 = don't extend)
		
		% easiest if we convert to indicator variables to do this:
		
		master_scrub_indicator = zeros(numvol,1);
		master_scrub_indicator(master_scrublist) = 1;

		if (~isempty(aap.tasklist.currenttask.settings.prekillzone))
			prekillzone = aap.tasklist.currenttask.settings.prekillzone;
		else
			prekillzone = 0;
		end
		
		if (prekillzone > 0)
			for index=1:length(master_scrub_indicator)
				if (master_scrub_indicator(index))
					istart = index - prekillzone;
					istart = max(istart,1);
					master_scrub_indicator(istart:index) = 1;
				end
			end
		end
				
		if (~isempty(aap.tasklist.currenttask.settings.postkillzone))
			postkillzone = aap.tasklist.currenttask.settings.postkillzone;
		else
			postkillzone = 0;
		end
		
		if (postkillzone > 0)
			for index=length(master_scrub_indicator):-1:1
				if (master_scrub_indicator(index))
					iend = index + postkillzone;
					iend = min(iend,length(master_scrub_indicator));
					master_scrub_indicator(index:iend) = 1;
				end
			end
		end		
		
		% process indicator list into scrub and keeper lists
							
		scrublist = []; keeperlist = [];
		for index = 1:length(master_scrub_indicator)
			if (master_scrub_indicator(index))
				scrublist = [ scrublist index ];
			else
				keeperlist = [ keeperlist index ];
			end
		end
		
		% make a simple diagnostic figure of the scrub for the report
		
		if (strcmp(aap.options.wheretoprocess, 'localsingle'))
			% centering trick only works if not cluster
			h = figure('Position',[0 0 800 100], 'Visible', 'off', 'MenuBar', 'none');
			movegui(h, 'center');
			set(h, 'Visible', 'on');
		else
			h = figure('Position',[0 0 800 100],'MenuBar','none');
		end
		master_scrub_indicator = [ master_scrub_indicator ; 0 ]; % for pcolor weirdness
		master_scrub_indicator = [ master_scrub_indicator' ; master_scrub_indicator' ];
		pcolor(master_scrub_indicator);
		axis off; colormap('flag');
		title(strrep([ aas_getsubjname(aap,subj) '/' aas_getsessname(aap,sess) ],'_','-'));
		set(h,'Renderer','opengl');
		set(findall(h,'Type','text'),'FontUnits','normalized');
		fname = fullfile(session_path, 'keeplist.jpg');
		print(h, '-djpeg', '-r150', fname);
		close(h);

			
        %  desc -----------------------------------------------------------
		
		% save the scrub and keeper lists. These are just a list of integers.
		% save( ) converts these to double, so we use dlmwrite( ) here which
		% doesn't. Either works with load( ) if you need the data later...
		      
        scrublist_fname = fullfile(session_path, 'scrublist.txt');
		dlmwrite(scrublist_fname, scrublist);
		aap = aas_desc_outputs(aap, subj, sess, 'scrublist', scrublist_fname);
				
		keeperlist_fname = fullfile(session_path, 'keeperlist.txt');
		dlmwrite(keeperlist_fname, keeperlist);
		aap = aas_desc_outputs(aap, subj, sess, 'keeperlist', keeperlist_fname);
        
        % aamod_firstlevel_model uses the stream name "listspikes" (not
        % "scrublist") so we have to save a copy of the scrub info to this 
		% streamname if we expect to use it in the GLM
		
		% also, have to kludge our scrublist into TSspikes and Mspikes --
		% aas_firstlevel_model_nuisance expects these to be column vectors
		% and does a union on them, so here we set both equal to transpose
		% of the scrublist (which gets generated above as a row vec)
		
		% also be advised aamod_firstlevel_model won't use this unless
		% you truthy the includespikes option
					
		TSspikes=scrublist'; Mspikes=scrublist';
		listspikes_fname = fullfile(session_path, 'listspikes.mat');
		save(listspikes_fname, 'TSspikes', 'Mspikes');
        aap = aas_desc_outputs(aap, subj, sess, 'listspikes', listspikes_fname);

		
    case 'checkrequirements'
        
    otherwise
        aas_log(aap, 1, sprintf('%s: Unknown task %s', mfilename, task));
		
		
end