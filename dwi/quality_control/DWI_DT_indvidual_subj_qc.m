%%% Dynamic tractography. Pipeline for assessing and cleaning data for
%%% J0044.

% Data in fsaverage space.
subdir = '/media/brian/hpc/data/Jeong-R01-Data/ESM-DTI_study/J0044/PreOP/VN/70deg/fsaverage/';

% Lateral occipital to FEF
tckfile = 'VN_70deg_qclen_fsa_18-2.tck';

%% Calcualte centroids at different distance thresholds to find the optimal number of exemplars.
% centroid_thresh = 4:2:20;
centroid_thresh = 5;
ncentroids = DWI_exemplars_test_thresholds({[subdir tckfile]},centroid_thresh,'/media/brian/hpc/data/Templates/2mm_fsaverage_brain.nii',true);

%% Compute streamline to centroid assignments
downsampled_tcks = [tckfile(1:end-4) '_downsampled.tck'];
exemplar_dir = [subdir 'exemplars/'];
exemplarfiles=cell(numel(centroid_thresh),1);
% Exemplar files calculated above.
for t = 1:numel(centroid_thresh)
    exemplarfiles{t} = [exemplar_dir tckfile(1:end-4) '_exemplar_t' num2str(centroid_thresh(t)) '.tck'];
end
% Centroid assignment thresholds to test
assign_thresh = [5:5:30 0];
[dist, centroid_ind] = DWI_assign_to_centroids(exemplarfiles, [subdir tckfile], assign_thresh);

%% Drop centroids with low numbers of associated streamlines.

% Percent of tracks
perc_cutoff = .1;
% Number of tracks
num_cutoff = 10;

% Centroid threshold
t = 5;

% Distance threshold
% Just check no distance thresholding for now.
dt = find(assign_thresh==0);

% load centroid file
centroids = read_mrtrix_tracks([exemplar_dir tckfile(1:end-4) '_exemplar_t' num2str(t) '.tck']);
centroids_out = centroids;
ncentroids = numel(centroids.data);

% Get number of tracks assigned to each centroid
tckcontrib = zeros(ncentroids,1);
for c = 1:ncentroids
    tckcontrib(c) = sum(centroid_ind{1}(:,dt)==c);
end
tckperc = tckcontrib./ntracks;

% Drop centroids below cutoff
% good_cents = find(tckperc >= perc_cutoff);
% bad_cents = find(tckperc < perc_cutoff);
good_cents = find(tckcontrib >= num_cutoff);
bad_cents = find(tckcontrib < num_cutoff);
centroids_out.data(bad_cents) = [];

% Write out pruned centroid set.
write_mrtrix_tracks(centroids_out,[exemplar_dir tckfile(1:end-4) '_exemplar_pruned_t' num2str(t) '.tck'])
% write_mrtrix_tracks(centroids_out,[exemplar_dir tckfile(1:end-4) '_exemplar_pruned_perc_t' num2str(t) '.tck'])

%% Grab tracks associated with retained centroids
% Load tracks
tcks = read_mrtrix_tracks([subdir tckfile]);
tcks_out = tcks;
% Drop tracks associated with bad centroids
tcks_out.data(ismember(centroid_ind{1}(:,dt),bad_cents)) = [];
% Write pruned dataset
write_mrtrix_tracks(tcks_out,[subdir tckfile(1:end-4) '_pruned.tck']);
% write_mrtrix_tracks(tcks_out,[subdir tckfile(1:end-4) '_pruned_perc.tck']);

%% Plot streamlines
figure; set(gcf,'Position',[200 200 1000 1000]);

% Coronal            
subplot(2,2,1)                        
hold on;
for t = 1:nbad
    plot3(centroids.data{t}(:,1),centroids.data{t}(:,2),centroids.data{t}(:,3),'r');
end
for t = 1:ngood
    plot3(centroids_out.data{t}(:,1),centroids_out.data{t}(:,2),centroids_out.data{t}(:,3),'b');                                
end
view([0 -1 0]);
xlabel('Left-Right'); ylabel('Posterior-Anterior'); zlabel('Inferior-Superior');
axis([-90 90 -90 90 -90 90 ])
title('Coronal')

% Axial
subplot(2,2,2)            
hold on;
for t = 1:nbad
    plot3(centroids.data{t}(:,1),centroids.data{t}(:,2),centroids.data{t}(:,3),'r');
end
for t = 1:ngood
    plot3(centroids_out.data{t}(:,1),centroids_out.data{t}(:,2),centroids_out.data{t}(:,3),'b');                
end
view([0 0 1])
xlabel('Left-Right'); ylabel('Posterior-Anterior'); zlabel('Inferior-Superior');
axis([-70 70 -70 70 -70 70 ])
title('Axial')

% Sagittal
subplot(2,2,3)             
hold on;
% Tracks
for t = 1:nbad
    plot3(centroids.data{t}(:,1),centroids.data{t}(:,2),centroids.data{t}(:,3),'r');
end
for t = 1:ngood
    plot3(centroids_out.data{t}(:,1),centroids_out.data{t}(:,2),centroids_out.data{t}(:,3),'b');
end
view([1 0 0])
xlabel('Left-Right'); ylabel('Posterior-Anterior'); zlabel('Inferior-Superior');
axis([-70 70 -70 70 -70 70 ])
title('Sagittal')

suplabel([subjects{n} ' | ' roi_names_short{r1} '-' roi_names_short{r2} ' | ' thresh],'t');

%% Plot number of streamlines assigned to each centroid

figure;
set(gcf,'Position',[200 0 1100 800])

ntracks = size(centroid_ind{1},1);
for t = 1:numel(centroid_thresh)
    
    % Get number of tracks assigned to each centroid
    tckcontrib = zeros(ncentroids(t),1);
    for c = 1:ncentroids(t)
        % Just check at assign_thresh = 5 for now.
        tckcontrib(c) = sum(centroid_ind{t}(:,dt)==c);
    end
    tckperc = tckcontrib./ntracks;
    [B, I] = sort(tckperc,'descend');
    clust = 1:ncentroids(t);
    
    % Plot
    subplot(numel(centroid_thresh),1,t)
    bar(B)
    %         ylim([0 100])
    xlim([0 numel(clust)])
    title(['Centroid generation threshold = ' num2str(centroid_thresh(t)) 'mm | Cent. assignment threshold = 5mm'])
    if t == numel(centroid_thresh)
        xlabel('Cluster #')
    end
    ylabel('% Tracks')
    set(gca,'XTickLabels',cellstr(num2str(I)))
end


