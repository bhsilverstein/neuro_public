%%% Calculate amplitude cross-correlations


%%% amplitude directed phase lag index (adPLI) pipeline.
%%% 1. Filter to 70 - 100 Hz
%%% 2. Hilbert transform, extract amplitude
%%% 3. Cut into events and z-score amplitude to baseline
%%% 4. Identify functionally involved channels based on gamma amplitude
%%% 5a. Plot averaged power spectrum of gamma amplitude signal
%%% 5b. Decide if the amplitude signal needs to be filtered to better
%%% isolate the desired signal. *Also will depend on moving window size.*
%%% 6. Using channel selections, return to step 2 output (continuous
%%% amplitude signal).
%%% 7. Hilbert transform the amplitude signal to extract its phase.
%%% 8. Re-apply event segmentation
%%% 9. Calculate phase differences as done in Joon-Young Moon's original
%%% code.
%%%
%%% This test pipeline is set up to use Fieldtrip output from
%%% preprocess_ieeg_auditory_m002ft.m


%%% *** NOTE FOR INTERPRETATION: NEGATIVE LAG = SIGNAL ONE LEADS SIGNAL TWO


%% Setup
clear


% Directories
% OSdir = 'E:/';
% OSdir = '/media/brian/3ACCF7A0CCF7549B/';
OSdir = '/media/user1/MyHDataStor11/brian/';
datadir=[OSdir 'ECoG_data/'];
xlsdir = [OSdir 'VI_FSA_ECS_all_pts/'];
outdir = [OSdir 'ECoG_Analysis/xcorr/roi/zscore/'];

% Subject list
sublist = dir(xlsdir);
sublist={sublist(~[sublist.isdir]).name};
sublist=cellfun(@(x) x(8:end-5),sublist,'UniformOutput',0);
nsubs=numel(sublist);

% Subset for talk
% sublist={'M13050'	'M14002'	'M14005'	'M14008'	'M14009'	'M14010'	'M14020'	'M14022'	'M14025'	'M14033' 'M14042'};
ucsdlist = {'M10013' 'M10044' 'M11017' 'M12041' 'M13009' 'M13013' 'M13037' 'M13045' 'M13047' 'M14005' 'M14009' 'M14020' 'M14025' 'M14047' 'M15013'};
nsubs=numel(sublist);

% Event labels
event_labels={'Stim On';'Stim Off';'Response Start';'Response End'};
alltrialtypes = [401 402 501 502];

% Select ROIs
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
nroi = numel(rois);
roi_pairs = nchoosek(rois,2);
npairs = length(roi_pairs);

%% Load
for n=1:nsubs
    % for n=2
    
    subid=sublist{n};
      
    
    %% Load data
    analysisdir = [datadir subid '/' subid '_auditory/analysis/'];
    fname = [subid '_ft_continuous'];
    if ~exist([analysisdir fname '.mat'],'file')
        continue
    else
        fprintf(['\nLoading ' subid '\n'])
        load([analysisdir fname '.mat']);
    end
    
    if ~isfield(subdata.elec,'elecpos')
        warning('% doesn''t have electroce coordinates. Skipping.', subid)
        continue
    end
    
    %% Filter
    cfgFilt = [];
    cfgFilt.bpfilter = 'yes';
    cfgFilt.bpfreq = [70 110]; % May need to increase maximum Hz to account for rolloff...
    cfgFilt.bpfilttype = 'firws';
    cfgFilt.bpfiltwintype = 'hamming';
    
    subdata = ft_preprocessing(cfgFilt,subdata);
    
    %% Hilbert and extract amplitude
    cfgHilb = [];
    cfgHilb.hilbert = 'abs'; % or absreal?
    
    subdata = ft_preprocessing(cfgHilb,subdata);
    
    %% Segement amplitude signal into events
    m00dir=[datadir subid '/' subid '_auditory/' subid '_auditory_m00/'];
    event_file=dir([m00dir '*.evt']); event_file=event_file.name;
    cfgTrial = [];
    cfgTrial.eventfile = event_file;
    cfgTrial.dname = m00dir;
    cfgTrial.trialfun = 'ieeg_auditory_ftevents';
    %     cfgTrial.dataset = [m00dir datafile{1}];
    
    % Define trials
    cfgTrial = ft_definetrial(cfgTrial);
    evtdata = ft_redefinetrial(cfgTrial, subdata);
    
    %% Z-score to baseline and select task-active trials
    
    trialtypes = unique(evtdata.trialinfo);
    sigchans = [];
    
    for t = 1:numel(trialtypes)
        
        % Average across trial types
        cfg = [];
        cfg.trials = find(evtdata.trialinfo==trialtypes(t));
        cfg.removemean = 'no';
        
        erp(t) = ft_timelockanalysis(cfg,evtdata);
        
        %% Grab prestim baseline data. Use this for all trial types.
        if trialtypes(t) == 401
            [~, base_end] = min(abs(erp(t).time - -0.1));
            [~, base_start] = min(abs(erp(t).time - -0.4));
            baseline=erp(t).avg(:,base_start:base_end);
        end
        
        %% Calculate z-score
        tmp = erp(t);
        tmp.ztobase = (erp(t).avg-mean(baseline,2))./std(baseline,0,2);
        
        zdata(t)=tmp;
        
        %% Select task-active channels
        % SetThreshold
        [~, stim_on] = min(abs(erp(t).time-0));
        ncomps=numel(erp(t).avg(:,stim_on:end));
        thresh=norminv([0.05/ncomps 1-(0.05/ncomps)]); % Bonferroni
        % Bonferroni is too crazy. Just use a=.001
        %     thresh=norminv([0.001 1-0.001]); % Bonferroni    bp_t = bp;
        
        % Apply threshold
        zthresh = zdata(t).ztobase;
        zthresh(abs(zthresh)<thresh(2))=0;
        
        %% Flag windows that are all above threshold
        % Look for 30 milliseconds of continuous activity.
        wsize = 30;
        sigwins = [];
        for w=1:length(erp(t).time)-wsize
            chk=zthresh(:,w:w+wsize-1)~=0;
            sigwins(:,w)=sum(chk,2)==wsize;
        end
        sigchans(:,t)=sum(sigwins,2);
    end
    
    
    %% Select channels and ROIs to limit channel combos
    % Labels of channels active sometime during the task in a event-related
    % fashion.
    siglabels = zdata(1).elec.label(logical(sum(sigchans,2)));
    siganat = zdata(1).elec.anat(logical(sum(sigchans,2)));
    sighemi = zdata(1).elec.hemi(logical(sum(sigchans,2)));
    
    %     anatInd = ismember(siganat,rois);
    % Find electrodes in the ROIs
    anatInd = ismember(subdata.elec.anat,rois);
    
    % Restrict to left hemisphere for now.
    %     hemiInd = ismember(sighemi,'Lt');
    hemiInd = ismember(subdata.elec.hemi,'Lt');
    
    
    %% Get indexes of channels to analyze
    %     sigchaninds = find(sum(sigchans,2));
    %     sigchaninds = find(sum(roiInd,2));
    %     sigchaninds = find(roiInd);
    sigchaninds = find(logical(sum(sigchans,2)) & hemiInd); % ***All the significant channels in the left hemisphere****
    
    %% See ieeg_psd_hgamp_aud.m for calculation of the power spectrum and ieeg_psd_plot_hgamp_aud.m for plotting the results.
    % Outcome of the above analysis suggests most of the activity is in the
    % 0-3 Hz band, with some very minor activity in the 5 - 20 Hz range.
    % Suggests it might be worth lowpass filtering the amplitude signal
    % down to 5 Hz to improve the SNR.
    
    %% Filter the amplitude signal, if desired.
    
    %     cfgFilt = [];
    %     cfgFilt.lpfilter = 'yes';
    %     cfgFilt.lpfreq = 5; % May need to increase maximum Hz to account for rolloff...
    %     cfgFilt.lpfilttype = 'firws';
    %     cfgFilt.lpfiltwintype = 'hamming';
    %
    %     subdata = ft_preprocessing(cfgFilt,subdata);
    
    %% Segment into trials
    %     evtdata = ft_redefinetrial(cfgTrial, subdata);
    
    %% Run cross correlation
    % Since the power oscillations occur at ~3 Hz, the period is ~1/3 sec.
    % which makes the maximum lag (1/3)/4 = 1/12 sec. Minimum window size
    % should be the length of the period, or event better, 2 periods: (1/3)*2 = 2/3 sec.
    % Min/max lags need some more thought. Some interaction between
    % potential delays based on conduction velocity+distance and limited by
    % the frequency. Since the peak frequency = 3, maybe set the max/min
    % delays at +/- 333 ms. One full cycle.
    
    %% xcorr parameters
    
    % Set up parameters
    params=[];
    params.sigchaninds = sigchaninds;
    params.chanpairinds = nchoosek(sigchaninds,2);                      % Channel pair indexes to analyze
    params.chanpairs = evtdata.elec.label(params.chanpairinds);         % channel pair labels
    params.anatpairs = evtdata.elec.anat(params.chanpairinds);          % Channel pair anatomy
    params.peakf = 3;                                                   % Observed peak power oscillation frequency
    params.wsize = floor(evtdata.fsample/params.peakf);                 % Window size
    params.overlap = 0.90;                                              % Overlap
    params.stepsize = floor((1-params.overlap) * params.wsize);         % Step size
    params.minlag = -floor(evtdata.fsample/params.peakf/2);               % Minimum lag
    params.maxlag = floor(evtdata.fsample/params.peakf/2);                % Maximum lag
    params.ntime = length(evtdata.time{1});                             % Number of timepoints
    params.wnum = numel(1:params.stepsize:params.ntime-params.wsize);   % Number of windows
    params.nshuff = 1000;                                               % Bootstrap number
    params.lags = params.minlag:params.maxlag;                           % Lags
  
    fprintf('\nStarting Cross-Correlations\n')
    fprintf('\n')
    
    %% Loop over trial types
    for t = 1:numel(trialtypes)
        
        % Grab trials of given type
        trials = find(evtdata.trialinfo==trialtypes(t));
        ntrials = numel(trials);
        
        %% Output vars
        xcorrdata=zeros(numel(params.lags),params.wnum,length(params.chanpairinds));
        
        % Window starting indexes
        wins = 1:params.stepsize:params.ntime-params.wsize;        
       
        %% Average within ROIs
        roidata = zeros(nroi,params.ntime);
        for r=1:nroi
            
            % Find channels in the current ROI
            anatInd = ismember(subdata.elec.anat,rois{r});

            % Grab the intersection of ROI, hemisphere, and significance
            roiInd = anatInd & hemiInd & logical(sum(sigchans,2));

            % Average and store
            if sum(roiInd)==1
                roidata(r,:) = erp(t).avg(roiInd,:);
            elseif sum(roiInd)<1
                roidata(r,:)=NaN;
            else
                roidata(r,:) = mean(erp(t).avg(roiInd,:));
            end
        end   
        
        %% Calculate cross-correlations
        
        tic        
        for w = 1:numel(wins)
            
            % Window data
            wdata = roidata(:,wins(w):wins(w)+params.wsize-1);
            
            % Mean-center (as in methods paper)
            wdata = wdata - mean(wdata,2);
            
            % cross-correlation
            [xtemp,lags]=xcorr(wdata',params.maxlag,'coeff');
            
            % Just keep upper triangle for the sake of file size.
            % Okay because ch1->ch2 = -(ch2->ch1)
            auto_ind = find(triu(ones(nroi),1)');
            xtemp = xtemp(:, auto_ind);            
            xcorrdata(:,w,:)=xtemp;
                        
            %% Surrogate data.
            shuff = zeros(numel(auto_ind),params.nshuff);
            parfor s = 1:params.nshuff
                shuffdata = fft_shuff(wdata');
                temp = xcorr(shuffdata,params.maxlag,'coeff');
                r = max(abs(temp(:,auto_ind)));
                shuff(:,s) = r;
            end 
            
            % Store summary variables
            shuff_m(w,:) = mean(shuff);
            shuff_s(w,:) = std(shuff);
                        
        end
        toc
        
        
        %% Save data
        if ~isdir(outdir)
            mkdir(outdir)
        end
        elec_labels=evtdata.elec;
        save([outdir 'ieeg_xcorr_aud_' subid '_t' num2str(trialtypes(t)) '.mat'],'xcorrdata','shuff_m','shuff_s','-v7.3')
        save([outdir 'ieeg_xcorr_aud_' subid '_t' num2str(trialtypes(t)) '_params.mat'],'params','elec_labels')
        
    end
end
























