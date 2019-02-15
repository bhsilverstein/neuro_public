%%% Import m00 files to fieldtrip

% Call ft_definetrials with ieeg_wordtype_ftevents.m
% Call ft_preprocessing for each m00 files
% Combine loaded data, update the data structure to include all channels
% Clean up channel names
% Drop flag channels to be excluded from the analysis

%% Setup
clear


% Directories
% OSdir = 'D:/';
OSdir = '/media/brian/3ACCF7A0CCF7549B/';
datadir=[OSdir 'ECoG_data/'];
xlsdir = [OSdir 'VI_FSA_ECS_all_pts/'];
T1dir = [OSdir 'T1_data/'];

% Subject list
sublist={'M15010' 'M16003' 'M18006'};

% Subset for talk
% sublist={'M13050'	'M14002'	'M14005'	'M14008'	'M14009'	'M14010'	'M14020'	'M14022'	'M14025'	'M14033' 'M14042'};
% ucsdlist = {'M10013' 'M10022' 'M10044' 'M11017' 'M11024' 'M12041' 'M13009' 'M13013' 'M13037' 'M13045' 'M13047' 'M14005' 'M14009' 'M14020' 'M14025' 'M14047' 'M15013'};
nsubs=numel(sublist);



%% Process data
for n=3
% for n=find(strcmp('M10044',sublist))
 %% Find files
    
    % PICK THE SUBJECT LIST OF INTEREST HERE******
    subid=sublist{n};                     
    
    % Find m00 files
    m00dir=[datadir subid '/' subid '_wordtype/' subid '_wordtype_m00/'];
    if exist(m00dir,'dir')
        cd(m00dir)
    else
        warning('No m00 folder. Skipping %s', subid);
        continue
    end             
    % First check for edf file. Indicates m00 files were replaced by an
    % .foc file. If no edf file is present, use the m00 files. If nothing
    % is found, skip the subject.
    datafile = dir('*.edf');
    datafile={datafile.name};
    numfiles=numel(datafile);
    filetype='edf';
    if numfiles<1
        datafile=dir('*.m00');    
        datafile={datafile.name};    
        if isempty(datafile)
            warning('No m00 or edf files. Skipping %s', subid);
            continue
        else
            numfiles = numel(datafile);
            filetype='m00';
            % Spit out files found
            for f=1:numfiles
                fprintf('\nFound %s \n',datafile{f});
            end
        end            
    else
         % Spit out files found
        for f=1:numfiles
            fprintf('\nFound %s \n',datafile{f});
        end
    end

    % Establish output directory
    savedir=[datadir subid '/' subid '_wordtype/analysis/'];
    if ~isdir(savedir)
        mkdir(savedir)
    end
    
    % Check if it's already been run
    if exist([savedir subid '_ft_preprocessed.mat'],'file')
        warning('Found an existing preprocessed file. Skipping %s',subid)
        continue
    end
    
        
    %% Define event structure    
    event_file=dir([m00dir '*.evt']); event_file=event_file.name;
    cfg = [];
    cfg.eventfile = event_file;
    cfg.dname = m00dir;
    cfg.trialfun = 'ieeg_wordtype_ftevents';    
    cfg.dataset = [m00dir datafile{1}];

    % Define trials
    cfg = ft_definetrial(cfg);
    
    %% Load m00 files
    % data.label = channel labels; data.time = cell array of time vectors for
    % each trial; data.trial = trial-wise cell array ofactual data in channels
    % x time; data.fsample = sample rate; data.sampleinfo = array of start/end
    % trial sample numbers (2 column array); data.trialinfo = 1D vector of
    % trial types;        

    for f=1:numfiles                      
        if f==1
            subdata = ft_preprocessing(cfg);
        else
            % Change metadata to target next datafile
            cfg.dataset=[m00dir datafile{f}];
            cfg.datafile=cfg.dataset;
            cfg.headerfile=cfg.datafile;
            nextfile = ft_preprocessing(cfg);
            
            % Combine data
            subdata.label = [subdata.label; nextfile.label];   
            subdata.trial = cellfun(@(x,y) cat(1,x,y), subdata.trial, nextfile.trial,'UniformOutput',false);    
            subdata.sampleinfo = [subdata.sampleinfo nextfile.sampleinfo];  
            
            % Update header info
            subdata.hdr.nChans = numel(subdata.label);
            subdata.hdr.label = subdata.label;
            subdata.hdr.chantype = [subdata.hdr.chantype; nextfile.hdr.chantype];
            subdata.hdr.chanunit = [subdata.hdr.chanunit; nextfile.hdr.chanunit];
            subdata.cfg.channel = subdata.label;
            if strcmp(filetype,'m00')
                subdata.hdr.dat = [subdata.hdr.dat; nextfile.hdr.dat];
            end
        end
    end

    %% Adjust channel names
    nchan = numel(subdata.label);
    for c = 1:nchan
        if strcmp(subdata.label{c}(1),'P') % Indicates 'POL ' prefix
            subdata.label{c}=subdata.label{c}(5:end);
        end
        if strcmp(subdata.label{c}(end),'V') % Indicates '-AV' suffix
            subdata.label{c}=subdata.label{c}(1:end-3);
        end
        if strcmp(subdata.label{c}(1),'*')
            subdata.label{c}=subdata.label{c}(2:end);
        end
    end
    
    %% Get channel info. Drop N/A channels (those that were disconnected)
    % Load excel file
    xlsfile=([xlsdir 'VI_FSA_' subid '.xlsx']);
    [chaninfo_num,chaninfo_txt,chaninfo_raw]=xlsread(xlsfile);    
    
    % Crop trailing NaNs
    chaninfo_raw=chaninfo_raw(1:size(chaninfo_txt,1),:);    
    
    % Tag functional iEEG channels in the Excel doc
    fsind = find(strcmp('fsaverage',chaninfo_raw(1,:)));
    okchans=~strcmp(chaninfo_txt(2:end,fsind),'N/A');
%     okchans(1)=false; % ignore header line
    
    %% Find bad channels    
    % Check for either 'bad' or 'Aud_bad/artifact'
    badnames={'bad';'Aud_bad/artifact';'artifact/overlap';'artifact+tumor?';'artifact/bad';'bad/overlap';'Bad/Artifact'};
    badind=1;
    while sum(strcmp(chaninfo_raw(1,badind),badnames))<1
        badind=badind+1;
        if badind>numel(chaninfo_raw(1,:))
            error('No bad channel column found. Check Excel sheet for consistency. Aborting %s',subid)
        end
    end    
    badchans = cell2mat(chaninfo_raw(2:end,badind));        
    
    % Check for either 'spiking' or 'spiking/spread'
    spikeind=1;
    while sum(strcmp(chaninfo_raw(1,spikeind),{'spiking';'spiking/spread';'spread'}))<1
        spikeind=spikeind+1;
        if spikeind>numel(chaninfo_raw(1,:))
            error('No spike channel column found. Check Excel sheet for consistency. Aborting %s',subid)
        end
    end
    spikechans = cell2mat(chaninfo_raw(2:end,spikeind));        
    sozchans = cell2mat(chaninfo_raw(2:end,strcmp('onset',chaninfo_raw(1,:))));
    
    allbad = (badchans | sozchans | spikechans | ~okchans);    
    allbadind=find(allbad);
%     reject(allbadind)=true;
    
    %% Grab corresponding channel labels
    allgoodind = find(~allbad);    
    labels_excel_good = chaninfo_txt(allgoodind+1,2); 
    hemi_excel_good =  chaninfo_txt(allgoodind+1,1);  
    anat_excel_good =  chaninfo_txt(allgoodind+1,4);             
    
    %% Identify channels loaded in the m00 files that match those in the Excel sheet (excluding channels labeled N/A).
    % Check to see if each loaded channel is present in the excel file.        
    [isfound, ~] = ismember(subdata.label, labels_excel_good);   
    
    %% Keep only channels on the good list    
    cfg = [];
    cfg.channel = subdata.label(isfound);
    subdata = ft_preprocessing(cfg,subdata);        

    %% Add hemisphere and anatomical location to EEG structure    
    [isfound, loc] = ismember(subdata.label, labels_excel_good);   
    for c = 1:numel(isfound)
        if isfound(c)
            subdata.elec.hemi{c,1}=hemi_excel_good(loc(c));
            subdata.elec.anat{c,1}=anat_excel_good(loc(c));
        end
    end
    subdata.elec.label=subdata.label;
    

    %% ******Add channel coordinates******
% Get coordinates and add them to the EEG structure
    if exist([T1dir subid '/electrode_voxel_coords_pial.mat'],'file')
        load([T1dir subid '/electrode_voxel_coords_pial.mat']);
        % Find channels in the EEG structure that match the loaded
        % coordinates.
        for c=1:numel(subdata.label)
            chk=strcmp(subdata.label{c},elec_label);
            if sum(chk)>0
                [~,ind]=max(chk);
                % Copy the coordinates of matching channels
                % ******** VERIFY ORIENTATION IS CORRECT ************
                % ******** SO FAR IT LOOKS GOOD, BUT YOU NEVER KNOW *****
                subdata.elec.elecpos(c,1) = elec_vert_coords(ind,1);
                subdata.elec.elecpos(c,2) = elec_vert_coords(ind,2);
                subdata.elec.elecpos(c,3) = elec_vert_coords(ind,3);
                % ***************************************************
            end
        end
        subdata.elec.chanpos=subdata.elec.elecpos;
    else
        warning('No coordinate file for %s.',subid)        
    end

    %% Save data
    save([savedir subid '_ft_preprocessed.mat'],'subdata','-v7.3')
    clear subdata nextfile
end










































