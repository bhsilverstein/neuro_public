%%% Aggregate streamlines across subjects. Provide a working directory, cell array subject list, tck file
%%% name, output directory for the group merged files, flags to do
%%% normalization to freesurfer atlas or merging across subjects.
%%% Example: DWI_normalize_and_merge([OSdir datadir], '/PreOP/VN/', {'J0001' 'J0002' 'J0003'}, 'VN_70deg_1-2.tck', '[OSdir '/DTI_Analysis/visual_network/group-level/merged/' prefix '/'], 1,)

function DWI_normalize_and_merge(expdir, projectdir, subjects, tckin, normfile, mergedir, normalize)

% Subject count
nsubs=length(subjects);

%% Normailze
if normalize
    for n=1:nsubs
        
        % Input dir
        indir = [expdir subjects{n} projectdir];
        
        % Check for data
        if ~exist([indir tckin],'file')
            warning([subjects{n} ' had no tck file. Skipped.'])
            continue
        end
        
        % Warp Image. If one doesn't exist, make it.  ****THIS SECTION HASN'T REALLY BEEN TESTED.      
        if ~exist([expdir subjects{n} '/PreOP/FOD_invmrtrixwarp_corrected.mif'],'file')
            fprintf(['\nNo warp image found for ' subjects{n} '. Attempting to calculate one\n'])
            % Create identity warp image
            if ~exist('/usr/local/freesurfer/subjects/fsaverage/mri.2mm/brain.mgz','file')
                error('Can''t find the FSA template. Modify this section locally so it directs to the correct template.')
            end
            moving = fullfile(pwd, 'T12B0anatomical.nii.gz');
            fixed  = '/usr/local/freesurfer/subjects/fsaverage/mri.2mm/brain.mgz';
            eval(sprintf('!warpinit -force  %s FOD_ID[].nii;', fixed));
            eval(sprintf('!WarpImageMultiTransform 3 FOD_ID0.nii FOD_invmrtrixwarp0.nii -R %s -i T1B02T1TEMPLATE0GenericAffine.mat T1B02T1TEMPLATE1InverseWarp.nii.gz;', moving));
            eval(sprintf('!WarpImageMultiTransform 3 FOD_ID1.nii FOD_invmrtrixwarp1.nii -R %s -i T1B02T1TEMPLATE0GenericAffine.mat T1B02T1TEMPLATE1InverseWarp.nii.gz;', moving));
            eval(sprintf('!WarpImageMultiTransform 3 FOD_ID2.nii FOD_invmrtrixwarp2.nii -R %s -i T1B02T1TEMPLATE0GenericAffine.mat T1B02T1TEMPLATE1InverseWarp.nii.gz;', moving));
            !warpcorrect -quiet -force FOD_invmrtrixwarp[].nii FOD_invmrtrixwarp_corrected.mif;
        end                    
        
        % Normalized output directory
        normdir = [indir 'fsaverage/'];
        if ~isdir(normdir)
            mkdir(normdir)
        end
                                
        % Normalize. ******* TO DO: UPDATE THIS SECTION SO THAT IT
        % DYNAMICALLY PULLS THE CORRECT WARP IMAGE.
        system(['tcknormalise ' indir tckin ' ' expdir subjects{n} '/PreOP/FOD_invmrtrixwarp_corrected.mif ' normdir normfile ' -quiet -force']);
    end
end

%% Build file list for merge
% Track list
tcklist=[];
k=1;
for n=1:nsubs        
    
    % Normalized dir
    normdir = [expdir subjects{n} projectdir 'fsaverage/'];
    
    % Add to list for group merge if exists
    if exist([normdir normfile],'file')
        tcklist=[tcklist normdir normfile ' '];    
        k=k+1;
    end           
end

%% Merge tracks across subject
% normdir = [OSdir '/DTI_Analysis/visual_network/group-level/merged/' prefix '/'];
if ~isdir(mergedir)
    mkdir(mergedir);
end
  
% Merge
system(['tckedit -force ' tcklist ' ' mergedir normfile]);
