function [aap,resp]=aamod_QA_threshplot(aap,task)

resp='';

switch task
	
    case 'report'
		
        fdiag = dir(fullfile(aap.acq_details.root,'*.jpg'));
        for d = 1:numel(fdiag)
            aap = aas_report_add(aap,[],'<table><tr><td>');
            aap=aas_report_addimage(aap,[],fullfile(aap.acq_details.root,fdiag(d).name));
            aap = aas_report_add(aap,[],'</td></tr></table>');
		end 
		
    case 'doit'
				
        nsub = length(aap.acq_details.subjects);
        nsess = length(aap.acq_details.sessions);
		
		QA_thresholds = [];  % for output stream
						
		% this needs good coverage of transition zone (ergo 0.01)

		quantiles_to_plot = [0:0.01:1];
		
		metric_names = { 'SV' 'DVARS' 'FD' };
			
		for mindex = 1:length(metric_names)
			
			metric_name = metric_names{mindex};
			
			h = figure(	'Units','points',...
						'Position',[0 0 800 800],...
						'Visible', 'off',...
						'Color', [1 1 1],...
						'Toolbar','none',...
						'NumberTitle','off',...
						'MenuBar','none' );

			movegui(h, 'center');
			set(h, 'Visible', 'on');
			
			metric_subject_values = zeros(nsub,length(quantiles_to_plot));
			
			% gather data from each subject across all sessions for this metric
			
			for m = 1:nsub
				
				subject_data = [ ];
				
				for n = 1:nsess
					QA_data = load(aas_getfiles_bystream(aap, m, n, 'QA'));
					subject_data = [ subject_data ; QA_data.(metric_name) ];
				end
			
				metric_values = quantile(subject_data, quantiles_to_plot);
				plot(metric_values, 100*(1-quantiles_to_plot), 'b', 'LineWidth', 1);
				hold on;

				metric_subject_values(m,:) = metric_values;
			
			end
			
			metric_mean_values = mean(metric_subject_values,1);
			
 			plot(metric_mean_values, 100*(1-quantiles_to_plot),'k', 'LineWidth', 4);
					
			% reference lines at 2.5%, 5%, etc. thresholds
						
			thresholds = [ 2.5 5.0 10 15 20 25 30 40 ];
			threshold_values = interp1(100*(1-quantiles_to_plot), metric_mean_values, thresholds);

			% trim long tail
			a = axis; axis([a(1) 2*threshold_values(1) 0 100]);

			os = [ 1 1.5 1.5 1.5 1.5 1.5 1.5 1.5 ];	% text offsets
			
			for tindex = 1:length(thresholds)
				
				a = axis;
				text_x = 0.9 * a(2);					
				text_y = thresholds(tindex);

				hline = refline([0 text_y]);
				hline.Color = [0.8 0.8 0.8]; hline.LineWidth = 1;
				
				text(text_x, text_y+os(tindex), num2str(threshold_values(tindex),'%.3g'),'FontName', 'Helvetica', 'FontSize', 14, 'Color', [0.8 0.8 0.8]);
				
			end
					
			ylabel('data loss (%)');
			xlabel(metric_name);
			title([ metric_name ' vs. data loss (n = ' num2str(nsub) ')']);
			
			set(gca,'FontName','Helvetica','FontWeight','bold','FontSize',16,'LineWidth',2);
			
			
			if (nsub > 1)

				% boxplot the data losses across all subjects @ the mean threshold values

				boxdata = zeros(nsub,length(thresholds));

				for m = 1:nsub
					boxdata(m,:) = interp1(metric_subject_values(m,:), 100*(1-quantiles_to_plot), threshold_values);
				end

				axes('Position',[0.4 0.55 0.4 0.3]);
				bp = boxplot(fliplr(boxdata),'Labels',num2str(fliplr(threshold_values)','%0.3g'),'LabelOrientation','inline');
				a = axis; axis([a(1) a(2) 0 100 ]);
				set(bp,{'linew'},{2});
				set(findobj(gca,'Type','text'), 'FontSize', 14, 'HorizontalAlignment', 'center');
				ylabel('data loss (%)');
				title({['Data loss across all subjects at indicated threshold'];[' ']});
				set(gca,'FontName','Helvetica','FontWeight','bold','FontSize',12,'LineWidth',2);
				txt = findobj(gca,'Type','text');
				set(txt, 'FontWeight', 'bold');
				
			end

			fname = fullfile(aap.acq_details.root,[ 'threshold_' metric_name '.jpg']);
			set(h,'Renderer','opengl');
			set(findall(h,'Type','text'),'FontUnits','normalized');
			print(h, '-djpeg', '-r150', fname);
			close(h);

			% save for stream output
			
			QA_thresholds = [ QA_thresholds ; threshold_values ];		
						
		end 
		
		
		% desc metric thresholds to QA_threshold stream

		outfile = fullfile(aas_getstudypath(aap), 'QA_thresholds.mat'); 
		save(outfile,'QA_thresholds');
		aap = aas_desc_outputs(aap, 'QA_thresholds', outfile);
 
		% done!

    case 'checkrequirements'
        
    otherwise
        aas_log(aap,1,sprintf('Unknown task %s',task));
		
end

end