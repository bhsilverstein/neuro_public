%% Setup
clear


% Directories
% OSdir = 'E:/';
% OSdir = '/media/brian/3ACCF7A0CCF7549B/';
% OSdir = 'C:/Users/bsilv/Documents/Asano Lab/';
% OSdir = '/media/user1/MyHDataStor11/brian/';
OSdir = '/media/brian/hpc/data/';

datadir=[OSdir 'ECoG_data/'];
xlsdir = [OSdir 'VI_FSA_ECS_all_pts/'];
outdir = [OSdir 'ECoG_Analysis/xcorr/roi/zscore/'];

% Subject list
sublist = dir(xlsdir);
sublist={sublist(~[sublist.isdir]).name};
sublist=cellfun(@(x) x(8:end-5),sublist,'UniformOutput',0);
nsubs=numel(sublist);

% Subset for talk
sublist={'M13050'	'M14002'	'M14005'	'M14008'	'M14009'	'M14010'	'M14020'	'M14022'	'M14025'	'M14033' 'M14042'};
ucsdlist = {'M10013' 'M10044' 'M11017' 'M13037' 'M13045' 'M13047' 'M14005' 'M14009' 'M14020' 'M14025' 'M14047' 'M15013'};
nsubs=numel(ucsdlist);

% Event labels
event_labels={'Stim On';'Stim Off';'Response Start';'Response End'};
alltrialtypes = [401 402 501 502];
ntypes = numel(alltrialtypes);

% ROIs
rois = sort({
    'pSTG_post';
    'pSTG_ant';
    'aSTG_post';
    'aSTG_ant';
    'inferiortemporal_pre';
    'inferiortemporal_post';
    'middletemporal_post';
    'middletemporal_pre';
    'precentral_inf';
    'precentral_mid';
    'precentral_sup';
    'postcentral_inf';
    'postcentral_mid';
    'postcentral_sup';
    'parsopercularis';
    'parstriangularis';
    'superiorfrontal_post';
    'superiorfrontal_ant';
    'caudalmiddlefrontal'
    });

rois_short = sort({
    'pSTG_post';
    'pSTG_ant';
    'aSTG_post';
    'aSTG_ant';
    'infTemp_pre';
    'infTemp_post';
    'midTemp_post';
    'midTemp_pre';
    'preCG_inf';
    'preCG_mid';
    'preCG_sup';
    'postCG_inf';
    'postCG_mid';
    'postCG_sup';
    'parsOp';
    'parsTri';
    'supFront_post';
    'supFront_ant';
    'caudMidFront'
    });
roi_pairs = nchoosek(rois,2);
npairs = length(roi_pairs);
nroi = numel(rois);
roi_pair_names = nchoosek(rois_short,2);

%% Load
nsubs=10;
for n=1:nsubs
    % for n=1
    
    subid=ucsdlist{n};
   
    fprintf(['\nLoading ' subid '\n'])
    for t=1:ntypes
        
        if ~exist([outdir 'ieeg_xcorr_aud_' subid '_t' num2str(alltrialtypes(t)) '.mat'],'file')
            warning('No results for %s, trial type: %d. Skipping',subid,alltrialtypes(t))
            xcorr_roi{t}(:,:,:,n) = NaN;
            shuff_m_all{t}(:,:,n)=NaN;
            shuff_s_all{t}(:,:,n)=NaN;
            continue
        end
        load([outdir 'ieeg_xcorr_aud_' subid '_t' num2str(alltrialtypes(t)) '.mat'])
        load([outdir 'ieeg_xcorr_aud_' subid '_t' num2str(alltrialtypes(t)) '_params.mat'])
                      
        % Sort by ROI
        [roi_sort, sortind] = sortrows(roi_pairs);
        xcorr_sort = xcorrdata(:,:,sortind);
        xcorr_roi{t}(:,:,:,n) = xcorr_sort;
        
        shuff_m_all{t}(:,:,n)=shuff_m;
        shuff_s_all{t}(:,:,n)=shuff_s;
    end
end

%% Trigger-relative time vector marking window midpoints.
triggers = [500 900 500 900]; % Double check these
wOnsets = 1:params.stepsize:params.ntime-params.wsize;
wCenter = wOnsets + params.wsize/2;

for t=1:ntypes
    [~,stim_on] = min(abs(wCenter - triggers(t)));
    times(:,t) = wCenter;
    times(:,t) = times(:,t)-times(stim_on,t);
end

%% Average and SD over trials
for t=1:ntypes
    xcorr_m{t} = squeeze(nanmean(xcorr_roi{t},4));
    xcorr_s{t} = squeeze(nanstd(xcorr_roi{t},0,4));
end

%% Lag limits
lagmax=params.maxlag;
lagmin=params.minlag;
% lagmax=100;
% lagmin=-100;
for t=1:ntypes    
    [~,maxind]=min(abs(lagmax-params.lags));
    [~,minind]=min(abs(lagmin-params.lags));
end
lagrange = params.lags(minind:maxind);


%% Peaks and lags
for t = 1:ntypes
    [peaks{t},lags{t}] = max(abs(xcorr_roi{t}(minind:maxind,:,:,:)),[],'omitnan');    
    peaks{t} = squeeze(peaks{t});
    lags{t} = squeeze(lags{t});
    lags{t} = lagrange(lags{t});
    for k=1:npairs
        for n=1:nsubs
            chk=~isnan(peaks{t}(:,k,n));
            if sum(chk(:))==0
                lags{t}(:,k,n)=NaN;
            end
        end
    end
    peaks_m{t} = squeeze(nanmean(peaks{t},3));
    peaks_s{t} = squeeze(nanstd(peaks{t},0,3));
    lags_m{t} = squeeze(nanmean(lags{t},3));
    lags_s{t} = squeeze(nanstd(lags{t},0,3));
end

%% Find peaks that pass threshold using surrogate data
thresh = norminv([0.05 1-0.05]);
thresh = thresh(2);

for t=1:ntypes
    % z-score
    z = (peaks{t} - shuff_m_all{t})./shuff_s_all{t};
    peak_h{t} = peaks{t} > thresh;
    peak_h{t}(isnan(peak_h{t}))=false;
end

%% Compare lags for significant peaks to zero
for t=1:ntypes
    ncomps = sum(peak_h{t}(:));
%     alpha = 0.05/ncomps;
    alpha = 0.001;
    lag_stats = permute(squeeze(lags{t}),[3 1 2]);                
                                       
    [lag_h{t},lag_p{t},lag_ci{t},temp] = ttest(lag_stats,0,'alpha',alpha);
    lag_t{t}=temp.tstat;
%     lag_h{t}=logical(lag_h{t});

    lag_h{t}(isnan(lag_h{t}))=false;
end

%% Plot raw data averaged over trials. Single pair.
% lagrange = [minind:maxind];

k=118;
figure;
for t=1:ntypes
    subplot(2,2,t)
    % heatmap
%     imagesc(times(:,t),lagrange,xcorr_roi{t}(:,:,k,1)); axis xy
    imagesc(times(:,t),lagrange,xcorr_m{t}(:,:,k)); axis xy
    % imagesc(times,lagshort,xcorrz_m(minind:maxind,:,c)); axis xy
    colorbar
    colormap jet
    xlabel('Time (ms)'); ylabel('Lags (ms)');
    title(event_labels{t})    
    hold on;
    line([0 0],[-333 333],'color','k')
    line([times(1,t)-params.wsize times(end,t)+params.wsize],[0 0],'color','k')
    
    % Overlap peaks
    lagplot = lags_m{t}(:,k);
    lagplot(~lag_h{t}(:,k) | ~peak_h{t}(:,k)) = NaN;
    plot(times(:,t),lagplot)
    plot(times(:,t),lags_m{t}(:,k),'--k')
    % boundedline(times,lags_tr_m(:,c),lags_tr_s(:,c),'alpha')    
end
plottitle = [strrep(roi_pair_names{k,1},'_','\_') '->' strrep(roi_pair_names{k,2},'_','\_')];
suplabel(plottitle,'t');

%% Plot by source ROI
% For each ROI, grab the "source" pairs.
% Grab the "target" pairs, and flip the lags. Not sure how to handle
% within-ROI pairs...

% For now, ignore within ROI pairs
for t=1:ntypes
% for t =1
    figure;
    p = 0;
    k = 0;
    for r1 = 1:nroi    
        for r2 = 1:nroi
            p=p+1;
            if r1>=r2
                continue
            end
            k=k+1;
            % Plot
            subplot(nroi,nroi,p)
            imagesc(times(:,t),lagrange,xcorr_m{t}(:,:,k)); axis xy
            colormap jet
            set(gca,'XTick',[],'YTick',[])
            plottitle = strrep(roi_pairs{k},'_','\_');
            title(plottitle,'fontsize',6)
            % Overlap peaks
            hold on;
            line([0 0],[params.minlag params.maxlag],'color','k')
            line([times(1,t)-params.wsize times(end,t)+params.wsize],[0 0],'color','k')
%             plot(times,lags{t}(1,:,k),'--k')
        end
    end
    suplabel(event_labels{t},'t');
end
















