%%% Determine the fiber cluster using exemplars. Thanks Minhee for the core
%%% MDF calculations! Takes centroid (exemplars) track files and a
%%% downsampled track file, computes the distance of each streamline to
%%% each centroid. Assigns each streamline to the nearest centroid. Can
%%% handle multiple centroid file/track file pairs. Returns a cell array
%%% of distances to nearest centroid for each streamline and the centroid
%%% to which it's assigned. NOTE: Hardcoded to use 12-point streamlines
%%% (the DIPY default for QuickBundles)
%%%
%%% Args:
%%% exemplarfiles: cell array of computed centroid tck files. Can input
%%% multiple for testing a range of centroid configurations.
%%% tckfile: .tck file containing sreamlines to compare to the loaded
%%% centroids.
%%% threshold: integer or array of integers. Distance (in mm) at which to
%%% threshold streamline to centroid assignment. If a streamline is not
%%% withing <threshold> of a centroid, it is assigned a centroid index of
%%% 0. Threshold of 0 = no threshold, each streamline gets assigned to the
%%% nearest centroid. Default = 0.

function [dist, centroid_ind] = DWI_assign_to_centroids(exemplarfiles, tckfile, threshold)

%% Check inputs

if nargin < 3
    threshold = 0;
end
if ~ischar(tckfile)
    error('tckfile arg must be a string. Input just one .tck file.')
end
if ~iscell(exemplarfiles)
    error('exemplarfiles arg must be a cell array')
end
nfiles = length(exemplarfiles);
nthresh = numel(threshold);

%% Ouputs
dist = cell(nfiles,1);
centroid_ind = dist;


%% Load tracks
[tck_dir, tck, ~] = fileparts(tckfile);
if isempty(tck_dir)
    tck_dir = './';
else
    tck_dir = [tck_dir '/'];
end
downsampled_file    = [tck_dir tck '.tck'];
load_downsampled 	= read_mrtrix_tracks(downsampled_file);
downsample          = load_downsampled.data;
n_points    = size(downsample{1},1);
n_tracts    = size(load_downsampled.data, 2);

%% Check length of streamlines. If it isn't 12, resample to 12.
if n_points ~=12
    fprintf('\n**Streamlines are not 12 points long. Resampling to 12 points.**\n')
    downsampled_file =  [tck_dir tck '_downsampled.tck'];
    system(['tckresample -force -num_points 12 ' tckfile ' ' downsampled_file]);
    load_downsampled 	= read_mrtrix_tracks(downsampled_file);
    downsample          = load_downsampled.data;
    n_points    = size(downsample{1},1);
    n_tracts    = size(load_downsampled.data, 2);
end

%% Compute distances to centroids
for f = 1:nfiles
    
    % Exemplar data files
    [exemplar_dir, exemplar, ~] = fileparts(exemplarfiles{f});
    if isempty(exemplar_dir)
        exemplar_dir = './';
    else
        exemplar_dir = [exemplar_dir '/'];
    end

    % Load centroid
    load_centroid = read_mrtrix_tracks([exemplar_dir exemplar '.tck']);
    class       = size(load_centroid.data, 2);
    index       = zeros(class, n_tracts);

    % Check if the loaded tracks are the same number of points as the
    % centroids
    if n_points ~= size(load_centroid.data{1},1)
        error(['The tracks are ' num2str(n_points) ' points long but the centroids are ' num2str(size(load_centroid.data{1},1)) ' points long.'...
                '\nCheck that the input tracks have been downsampled appropriately.'])
    end
    
    % Calculate distances (MDF)
    for i = 1 : class
        centroids       = load_centroid.data{i};

        % Loop over tracks
        for j = 1 : n_tracts
            Ddirect     = zeros(1, n_points);
            FDdirect    = Ddirect;
            down        = downsample{j};
            flip_v      = flipud(down);

            % Calculate distance from centroids at each point
            for k = 1 : n_points
                Ddirect(k)  = norm(centroids(k, :) - down(k, :));
                FDdirect(k) = norm(centroids(k, :) - flip_v(k, :));
            end

            % Get distance to nearest centroid.
            index(i, j)     = min((1/n_points)*sum(Ddirect), (1/n_points)*sum(FDdirect)) ;
            clear down
        end
%         clear centroids
    end
    dist{f} = index;
    
    %% Assign streamlines to centroids.
    %A track gets associated with the
    % centroid to which it is closest.  If a track is not within
    % threshold for any cluster, it is dropped (cluster ind=0)
    for t = 1:nthresh        
        class_num            = zeros(1,  n_tracts);
        for k = 1 : n_tracts
            if threshold(t)==0
                class_num(k) = find(index(:,k) == min(index(:, k)));
            else
                if sum(index(:,k) < threshold(t))
                    class_num(k) = find(index(:,k) == min(index(:, k)));
                else
                    class_num(k) = 0;
                end
            end
        end
        centroid_ind{f}(:,t) = class_num;
    end        
end

























