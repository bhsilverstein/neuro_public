%%% Test clustering threshold parameter to see how it changes the number of
%%% streamline bundles identified. Stops when it reaches one cluster.
%%% Plot results. Converts tck files to trk files if needed. Can input a
%%% list of tck files to test multiple tracks. Saves generated files in
%%% ./exemplars/ in the directory the tracks are pulled from.
%%% This code depends on the modified DIPY function centroid_test.py.
%%% Make sure it's in the directory where this function is called from or 
%%% modify it's path below. Also depends on convert_tracks_format.m, 
%%% Trk_load.m, get_header.m, and MRtrix I/O functions.
%%%
%%% Args:
%%% tckfiles = cell array of .tck file names
%%% thresholds = array of thresholds to test in mm
%%% Reference image for converting to .trk if needed. *NO .gz compression.*
%%% doplot = boolean flag to plot results or not.
%%%
%%% Example: ncentroids = DWI_exemplars_test_thresholds({<filepath>/<track.tck>; <filepath>/<track.tck>},[5:10 15:20],<fullpath>/<subjT1 or FSA tempalte>,true)

function ncentroids = DWI_exemplars_test_thresholds(tckfiles,thresholds,refimage,doplot)

%% Check inputs
if ~iscell(tckfiles)
    error('tckfiles arg must be a cell array')
end

nfiles = length(tckfiles);
ncentroids = zeros(nfiles,numel(thresholds));

%% Create output directory
[filepath,~,~] = fileparts(tckfiles{1});
if isempty(filepath)
    filepath = '.';
end
outdir = [filepath '/exemplars/'];
if ~isdir(outdir)
    mkdir(outdir);
end
fprintf(['\nSaving exemplars to: ' outdir '\n'])

%% Compute centroids
for f = 1:nfiles

    % Trackfile
    [tckpath,tck,~] = fileparts(tckfiles{f});
    if isempty(tckpath)
        tckpath = './';
    else
        tckpath = [tckpath '/'];
    end
    trk = [tck '.trk'];
    
    %% Convert to trk if needed
    if ~exist(trk,'file')
        fprintf('\nConverting to .trk...\n')
        convert_tracks_format([tckpath tck '.tck'],'tck2trk',refimage);
    else
        fprintf('\n.trk file already exists. Skipping conversion.\n')
    end

    %% Compute centroid with DIPY            
    
    % Compute centroids over a range of thresholds until converging
    % around one cluster    
    clustconv=0;
    t=0;
    while clustconv==0
        t=t+1;
        thresh=num2str(thresholds(t));

        % Exemplar output file
        exemplar = [tck '_exemplar_t' thresh];

        % Run centroid_test.py
        fprintf(['\nThresh = ' thresh '\n'])
        system(['python3 ~/matlab_code/dti/postprocessing/bundles/centroid_test.py ' tckpath trk ' ' outdir exemplar '.trk ' thresh]);

        % Convert back to tck
%         if ~exist([exemplar '_t' thresh '.tck'],'file')
        convert_tracks_format([outdir exemplar '.trk'],'trk2tck',refimage);
%         end

        % Load computed exemplars.
        load_centroid = read_mrtrix_tracks([outdir exemplar '.tck']);
        ncentroids(f,t) = length(load_centroid.data);

        % If the algorithm has reached convergence (just one cluster)
        % then stop
        if ncentroids(f,t)<2
            clustconv=1;
        elseif t==numel(thresholds)
            warning('End of input thresholds. Did not converge to one exemplar.')
            clustconv=1;
        end
    end
end

%% Plot
if doplot
    figure
    for k=1:nfiles
        subplot(ceil(sqrt(nfiles)),ceil(sqrt(nfiles)),k)
        plot(thresholds,ncentroids(k,:),'-o')
        xlim([thresholds(1)-1 thresholds(end)+1])
    %     ylim([0 3])
%         line([33 33],[1 max(ncentroids(k,:))],'LineStyle','--')
        [~,tck,~] = fileparts(tckfiles{f});
        title(tck)    
        xlabel('Threshold')
        ylabel('# of exemplars')
    end

    %%
%     figure; hold on
%     plot(thresholds,ncentroids,'k')
%     plot(thresholds,ncentroids([2 22 23],:),'r')
%     line([thresholds(1) thresholds(end)],[2 2],'LineStyle','--')
% 
%     line([33 33],[1 7],'LineStyle','--')
end



















