%%% Apply post-processing quality control metrics to the visual network DWI
%%% data.



%% Setup
% Local file structure
OSdir='/media/user1/MyHDataStor11/brian/';
subdir = 'Jeong-R01-Data/ESM-DTI_study/';
expdir = [OSdir subdir];

% Load subject numbers
% All DWI Subjects
fid=fopen('/media/user1/MyHDataStor11/brian/Jeong-R01-Data/DTI_SubjectList.txt');
subjects=textscan(fid,'%s'); subjects=subjects{1};
nsubs=length(subjects);

% Subjects preprocessed by Minhee with 5 million seed whole brain
% tractography
% [~, txt, ~] = xlsread([OSdir 'DTI_Analysis/visual_network/Visual_Network_Subjects.xlsx']);
% subjects = txt(2:end,1);
% nsubs = numel(subjects);

% Subjects with CCEP and DWI data
% fid=fopen([OSdir 'CCEP_data/CCEP_DWI_Subjects.txt']);
% subjects=textscan(fid,'%s'); subjects=subjects{1};
% nsubs=length(subjects);

% Nodes of interest
pairs = [15 18;
         16 18;
         18 3;
         3  4;
         4  2
         15 2
         16 2
         18 2
         3  2
         ];   
    
prefix = '70deg';

for n=1:nsubs  

    datadir = [expdir subjects{n} '/PreOP/VN/' prefix '/'];
    
    if ~exist([datadir 'VN_' prefix '_3-2.tck'],'file')
        continue
    end
    
    %% Apply length threshold
    for p=1:size(pairs,1)
        
        tckin = ['VN_70deg_' num2str(pairs(p,1)) '-' num2str(pairs(p,2)) '.tck'];
        tckout = ['VN_70deg_qclen_' num2str(pairs(p,1)) '-' num2str(pairs(p,2)) '.tck'];
        
        [t, len_m, len_s] = DWI_qc_length(datadir, tckin, tckout, false);
    end
        
    %% Downsample to 100 points and convert to vtk
    for p=1:size(pairs,1)        
        
        % Downsample
        tckin = ['VN_70deg_qclen_' num2str(pairs(p,1)) '-' num2str(pairs(p,2))];
        system(['tckresample -force -num_points 100 ' datadir tckin '.tck ' datadir tckin '_downsampled.tck']);        
        
        % Convert to VTK
        downsampled = ['VN_70deg_qclen_' num2str(pairs(p,1)) '-' num2str(pairs(p,2)) '_downsampled'];
        system(['tckconvert ' datadir downsampled '.tck ' datadir downsampled '.vtk']);        
    end        
    
    %% Compute torsion and curvature using TRAFIC
    for p=1:size(pairs,1)
        trafic = '/media/user1/MyHDataStor11/Data/CNN-Haotian/final_scripts/deep_FL_CL_ATT/preprocessing_tc/TRAFIC.py';
        system(['python trafic.get_torsion_and_curvature ' datadir downsampled '.vtk']);
    end
    
    
end

%% Normalize and merge for plotting
for p=1:size(pairs,1)
% for p=1
    tckin = ['VN_70deg_qclen_' num2str(pairs(p,1)) '-' num2str(pairs(p,2)) '.tck'];
%     normfile = ['VN_70deg_qclen_fsa_' num2str(pairs(p,1)) '-' num2str(pairs(p,2)) '.tck'];
    normfile = ['VN_70deg_' num2str(pairs(p,1)) '-' num2str(pairs(p,2)) '_fsa.tck'];
    DWI_normalize_and_merge(expdir, '/PreOP/VN/70deg/', subjects, tckin, normfile, [OSdir '/DTI_Analysis/visual_network/group-level/merged/' prefix '/'], 0)
end

























