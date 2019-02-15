%%% TFA with Fieldtrip

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

%% Electrodes of interest per subject. Grabbed from Zahraa's original analysis. Results from this pipeline may vary.
elecs{1} = {'A16'; 'B2'; 'B3';};
elecs{2} = {'B61'; 'C18'; 'B2'; 'B7'; 'B45';};
elecs{3} = {'A41'; 'B19'; 'B1'; 'B9'; 'B27'; 'B34'; 'B21'; 'B62';};
       
%% Analysis
% for n=1:nsubs
for n=2
    
    subid=sublist{n};
       
    %% Load data
    analysisdir = [datadir subid '/' subid '_' project '/analysis/'];
    fname = [subid '_ft_preprocessed'];
    fprintf(['\nLoading ' fname '...\n'])
    load([analysisdir fname '.mat']);

    
    %% Time frequency analysis
    
    % Select trial type
    trialtype = [101 102 301 302 401 402 501 502 601 602 701 702 201];
    
    for t = 1:numel(trialtype)
        
        % Double check the event exists
        if ~any(subdata.trialinfo==trialtype(t))
            warning('Event type %d not found. Moving on.',trialtype(t))
            continue
        end
        
        % Event times
        time = subdata.time{find(subdata.trialinfo==trialtype(t),1)};
        [~, trigger_ind] = min(abs(time-0));        
        
        % Time range to analyze
        wsize = 0.2; % Seconds        
        tstep = 0.01; % Seconds
        tstart = time(1) + wsize/2;
        tend = time(end) - wsize/2 - tstep;
        
        
        % Set up configuration        
        cfg = [];
        cfg.method      = 'mtmconvol';
%         cfg.toi         = -0.41:0.1:0.8;
        cfg.toi         = tstart:tstep:tend;
        cfg.foi         = 15:5:150;
        cfg.t_ftimwin   = ones(length(cfg.foi),1).*wsize;
        cfg.taper       = 'hanning';
        cfg.output      = 'pow';
        cfg.pad         = 'nextpow2';
        cfg.keeptrials  = 'yes';
        cfg.trials = find(subdata.trialinfo==trialtype(t));
        
        tmp = ft_freqanalysis(cfg, subdata);       
        
        tmp.trialtype=trialtype(t);
        freq(t) = tmp;
        
        tmp = cfg;
        tmp.wsize=wsize;
        tmp.tstep=tstep;
        tmp.tstart=tstart;
        tmp.tend=tend;
        
        freqparams(t)=tmp;
    end   
       
    %% Calculate percent change. Store in a structure array.
      
    % Baseline TF
    baseline = freq(trialtype == 201);
    
    for t = 1:numel(trialtype)
        
        % Don't need to run baseline
        if trialtype(t) == 201
            continue
        end
                        
        
        
        % Build output
        tmp = [];
        tmp.label = freq(t).label;     
        tmp.freq = freq(t).freq;
        tmp.time = freq(t).time;
        tmp.elec = freq(t).elec;
        tmp.cfg = freq(t).cfg;
        tmp.trialtype = freq(t).trialtype;
        
        % Baseline correct each event to the previously occuring baseline.
        trial_inds = find(subdata.trialinfo==trialtype(t));
        allbase_inds = find(subdata.trialinfo==201);
        for k = 1:numel(trial_inds)            
            
            % Find matching baseline
            base_trial = find(subdata.trialinfo(1:trial_inds(k))==201,1,'last');            
            base_ind = find(allbase_inds==base_trial);
            
            % Divide by the mean and standard deviation of baseline over the time dimension
            tmp.ztobase(k,:,:,:) = (freq(t).powspctrm(k,:,:,:)-mean(baseline.powspctrm(base_ind,:,:,:),4))./std(baseline.powspctrm(base_ind,:,:,:),0,4);        
        end
        
        freq_blc(t)=tmp;        
        baselineparams = 'Z-scored to the most recent baseline period';
        % Fieldtrip function that does the same thing
        %     cfg              = [];
        %     cfg.baseline     = [-.4 -.1]; 
        %     cfg.baselinetype = 'relchange'; 
        %     freq_perc = ft_freqbaseline(cfg, freq);
        %     baselineparams=cfg;
    end
                
    %% Grab random subset of normal word trials, spread equally across the word types. Combine filler word trial types and normal word trial types.
    ntrials=arrayfun(@(x) size(x.ztobase,1), freq_blc);
    total_filler = sum(ntrials(trialtype == 601 | trialtype == 701));        
    
    % Trial, electrode, frequency, time. Cell = trialtype (1 = normal word onset; 2 =
    % normal word offset; 3 = filler word onset; 4 = fillerword offset)
    psd{1}=[]; psd{2}=[]; psd{3}=[]; psd{4}=[];
    for t=1:numel(ntrials)
        if ~ismember(trialtype(t),[601 602 701 702])
            % Odd trialtypes are onsets, even trial types are offsets
            if ~mod(trialtype(t),2)
                subset_inds = randperm(ntrials(t),round(total_filler/4));
                psd{1} = [psd{1}; freq_blc(t).ztobase(subset_inds,:,:,:)];
            else
                subset_inds = randperm(ntrials(t),round(total_filler/4));
                psd{2} = [psd{2}; freq_blc(t).ztobase(subset_inds,:,:,:)];
            end
        else            
            if ~mod(trialtype(t),2)
                psd{3} = [psd{3}; freq_blc(t).ztobase];
            else
                psd{4} = [psd{4}; freq_blc(t).ztobase];
            end
        end
    end
    
    %% Average over frequency
    [~,f1] = min(abs(freq_blc(1).freq-70));
    [~,f2] = min(abs(freq_blc(1).freq-110));
    for t = 1:numel(psd)
        gamma{t} = squeeze(mean(psd{t}(:,:,f1:f2,:),3));
    end
    
    %% Run t-tests 
    bonf = 0.05/(71*numel(freq_blc(1).time));
    % Window of interest
    [~,t1] = min(abs(freq_blc(1).time - -0.500));
    [~,t2] = min(abs(freq_blc(1).time - 0.500));
    for e = 1:numel(freq_blc(1).label)
        for t = 1:numel(freq_blc(1).time)
            % Onsets
            [~,p(e,t,1),~,~] = ttest2(gamma{1}(:,e,t),gamma{3}(:,e,t));
            % Offsets
            [~,p(e,t,2),~,~] = ttest2(gamma{2}(:,e,t),gamma{4}(:,e,t));
        end
    end
    % FDR Corrections. Just on the window of interest.
    for e = 1:numel(freq_blc(1).label)
%         [~, p_cor(e,:,1)] = fdr(p(e,t1:t2,1),0.05);
%         [~, p_cor(e,:,2)] = fdr(p(e,t1:t2,2),0.05);
        [pthresh(e,1), p_cor(e,:,1)] = fdr(p(e,t1:t2,1),0.05);
        [pthresh(e,2), p_cor(e,:,2)] = fdr(p(e,t1:t2,2),0.05);
    end
       
    % Sig. channels
    [sigchans,~,~] = ind2sub(size(p_cor),find(p_cor));
    sigchans = unique(sigchans);
    
    %% Average over trials
    for t=1:numel(psd)
        gamma_erp{t} = squeeze(mean(gamma{t}));
    end        
    
    
    %% Plot erp
%     elec_inds = zeros(numel(elecs{n}),1);    
%     for e = 1:numel(elecs{n})
%         elec_inds(e) = find(strcmp(elecs{n}{e},freq_blc(1).label));
%     end
    elec_inds = sigchans;
    
    time = freq_blc(1).time;
    
    k=0;
    figure;    
    for e = 1:2:numel(elecs{n})*2
        k=k+1;
        % Onsets
        subplot(numel(elecs{n}),2,e); hold on;
                
        plot(time,gamma_erp{1}(elec_inds(k),:),'--')
        plot(time,gamma_erp{3}(elec_inds(k),:))
        
        % Plot limits
        pmin = min([gamma_erp{1}(elec_inds(k),:) gamma_erp{3}(elec_inds(k),:)]) - 0.5;
        pmax = max([gamma_erp{1}(elec_inds(k),:) gamma_erp{3}(elec_inds(k),:)]) + 0.5;
        
        % Sig. line
        sig = p_cor(elec_inds(k),:,1) .* (pmax - 0.1);
        sig(sig==0) = NaN;
        plot(time(t1:t2),sig,'r')
        
        % Event line
        line([0 0],[pmin pmax],'color','k','linestyle','--')
        
        
        xlim([time(1) time(end)])
        ylim([pmin pmax])
        ylabel('Z-score')
        zlabel('Time')
        legend('Normal','Filler')
        title(['Onsets | ' elecs{n}{k} ' | ' freq_blc(1).elec.anat{elec_inds(k)}{1}])
        
        % Offsets
        subplot(numel(elecs{n}),2,e+1); hold on;
                
        plot(time,gamma_erp{2}(elec_inds(k),:),'--')
        plot(time,gamma_erp{4}(elec_inds(k),:))
        
        pmin = min([gamma_erp{2}(elec_inds(k),:) gamma_erp{4}(elec_inds(k),:)]) - 0.5;
        pmax = max([gamma_erp{2}(elec_inds(k),:) gamma_erp{4}(elec_inds(k),:)]) + 0.5;
        
        
        % Sig. line
        sig = p_cor(elec_inds(k),:,2) .* (pmax - 0.1);
        sig(sig==0) = NaN;
        plot(time(t1:t2),sig,'r')
        
        line([0 0],[pmin pmax],'color','k','linestyle','--')
        xlim([time(1) time(end)])
        ylim([pmin pmax])
        xlabel('Z-score')
        ylabel('Time')
        legend('Normal','Filler')
        
        title(['Offsets | ' elecs{n}{k} ' | ' freq_blc(1).elec.anat{elec_inds(k)}{1}])
        
    end
    
    %% Save data
%     save([analysisdir subid '_freq.mat'],'freqparams','freq_blc','baselineparams','-v7.3')        

%     ft_freqstatistics
    
end



























