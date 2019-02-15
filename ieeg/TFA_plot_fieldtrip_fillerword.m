%%% Plot output of TFA from Fieldtrip.

%% Setup
clear


% Directories
% OSdir = 'E:/';
OSdir = '/media/brian/3ACCF7A0CCF7549B/';
datadir=[OSdir 'ECoG_data/'];
xlsdir = [OSdir 'VI_FSA_ECS_all_pts/'];

sublist={'M15010' 'M16003' 'M18006'};
nsubs=numel(sublist);

project = 'wordtype';
outdir = '/media/brian/hpc/data/ECoG_Analysis/FillerWord/TFA/'; % where to save figures

% for n=1:nsubs
for n=1
    
    subid=ucsdlist{n};
    
%     % Just do the test set.
%     if ~any(strcmp(subid,subset))
%         continue
%     end
    
    %% Load data
    analysisdir = [datadir subid '/' subid '_' project '/analysis/'];
    fname = [subid '_freq'];
    load([analysisdir subid fname '.mat']);                    
    
    %% Plot percent change in each electrode bank
     trialtypes = {'Normal Onset';'Normal Offset';'Filler Onset';'Filler Offset'};     
     banklist = {'A';'B';'C'};
%     for t = 1:numel(freq)
    for t = 1:numel(psd)
        for c=1:numel(freq(t).elec.label)
            banks(c)=find(strcmp(freq(t).elec.label{c}(1),banklist));
        end
        k=0;
        for b = 1:length(unique(banks))
            % Number of channels in the bank
            chans2plot = find(banks==b);
            nchans=numel(chans2plot);
            if nchans>0
                % Determine a nice rectangle shape for the subplots
                area=nchans;
                while sum(isprime([area ceil(area/2)]))>0
                    area=area+1;
                end
                for w=ceil(sqrt(area)):-1:1
                    if ~mod(area,w)
                        % height
                        h = area/w;
                        break
                    end
                end
                % Plot
                figure; set(gcf,'Position',[200 200 2200 1200],'visible','off')
                for c = 1:nchans
                    k=k+1;
                    subplot(h,w,c)
                    hold on;
                    chanind = chans2plot(c);

                    %****** ACTUAL PLOT IS HERE********
%                     imagesc(freq(t).time,freq(t).freq,squeeze(freq_blc(t).powspctrm_perc(chanind,:,:))); axis xy;
                    imagesc(freq(t).time,freq(t).freq,squeeze(mean(psd{t}(:,chanind,:,:)))); axis xy;
                    colormap jet;
                    colorbar
                    %***** Axis labels************
                    if c>(h-1)*w
                        xlabel('Time (seconds)')
                    else
                        set(gca,'XTick',[],'XTickLabels',[])
                    end
                    if sum(c==1:w:(w*h))
%                         ylabel({'Hz';'% Change'}); % color label('Power V^2/Hz')
                        ylabel({'Hz';'Z-score'}); % color label('Power V^2/Hz')
                        %                 else
                        %                     set(gca,'YTick',[],'YTickLabels',[])
                    end

                    %****** Lines, axes, and titles *********
                    % set(gca,'XTick',ft,'XTickLabel',
                    caxis([0 6.00])
                    xlim([freq(t).time(1) freq(t).time(end)]) % x-axis display range
                    ylim([freq(t).freq(1) freq(t).freq(end)])
                    title(strjoin([freq(t).elec.label{chanind} freq(t).elec.hemi{chanind} freq(t).elec.anat{chanind}]))
                    line([0 0],[1 max(freq(t).freq)],'LineStyle','-','Color','k','LineWidth',1)
                    %                 line([1.8 1.8],[1 numel(freq(t).freq)],'LineStyle','--','Color','k','LineWidth',2)
                end
                suplabel(trialtypes{t},'t');
                saveas(gcf,[outdir subid '_tfa_' banklist{b} 'Bank_' trialtypes{t} '.jpg'],'jpg')
                close
            end
        end
    end
end