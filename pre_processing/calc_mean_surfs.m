% calc_mean_surfs - generates mean cortical surface FC from all subjects in
% ALSPAC dataset for specified atlases.
%
% Output directories assumed to already exist.
% NaN values only seem to appear in subcortex / cerebellum. Individual NaNs
% are removed while the rest of the FC data is still used (i.e. there is a
% different averaging constant for each parcel).
%
% Author: Chris Norman, Cardiff University
% Email: NormanC4@cardiff.ac.uk
% July 2022; Last revision: 29/07/22
%
%% Paths and settings
Data_Dir = '/scratch/scw1648/proj_cn/data/mica_processed/micapipe';
Results_Dir = '/scratch/scw1648/proj_cn/data/mica_processed/micapipe/mean_surfs';

FileIndicators = importdata('/scratch/scw1648/proj_cn/data/FileIndicators.csv',',');
% Account for missing result if not yet corrected!!
FileIndicators.data(FileIndicators.data(:,1)==1153,3) = 0;

Spaces = {'fsaverage5'};
Atlases_full_list = {...
    'aparc', ...
    'economo', ...
    'glasser-360', ...
    'schaefer-100', ...
    'schaefer-200', ...
    'schaefer-300', ...
    'schaefer-400', ...
    'schaefer-500', ...
    'schaefer-600', ...
    'schaefer-700', ...
    'schaefer-800', ...
    'schaefer-900', ...
    'schaefer-1000', ...
    'vosdewael-100', ...
    'vosdewael-200', ...
    'vosdewael-300', ...
    'vosdewael-400', ...
    };
Atlases = {'glasser-360'}; % or just Atlases_full_list
Descs = {'FC'};
Sessions = {'rest','movie'};

%% Main Loop

for Spacei = 1:length(Spaces); Space = Spaces{Spacei};
    for Atlasi = 1:length(Atlases); Atlas = Atlases{Atlasi};
        for Desci = 1:length(Descs); Desc = Descs{Desci};
            for Sessioni = 1:length(Sessions); Session = Sessions{Sessioni};
                
                % Which FileIndicator column corresponds to this session
                if strcmp(Session,'rest')
                    Indi_col = 3;
                elseif strcmp(Session,'movie')
                    Indi_col = 5;
                else
                    error(['Session "',Session,'" not recognised. Available options are "rest" and "movie".'])
                end
                
                % Begin data averaging
                disp(['Averaging: ',Space,'_',Atlas,'_',Desc,'_',Session,'...'])
                Initialise = 1;
                for IDi = 1:size(FileIndicators.data,1)
                    
                    % Only add subject if data file exists
                    if FileIndicators.data(IDi,Indi_col) == 0
                        continue
                    end
                    
                    % Load data
                    Sub = ['sub-',num2str(FileIndicators.data(IDi,1))];
                    File_Path = fullfile(Data_Dir,Sub,'func',Session,'surfaces',...
                        [Sub,'_rsfmri_space-',Space,'_atlas-',Atlas,'_desc-',Desc,'.txt']);
                    Sub_data = dlmread(File_Path,' ');
                    % Make FC matrix symmetric
                    if strcmp(Desc,'FC')
                        Sub_data = Sub_data + Sub_data' - diag(diag(Sub_data));
                    end
                    
                    % First case, use to initialise mean variable
                    if Initialise
                        mean_surf = zeros(size(Sub_data));
                        N = zeros(size(Sub_data));
                        Initialise = 0;
                    end
                    
                    % Find NaNs, remove from data, update normalisation
                    % constants (N) for each ROI accordingly
                    NaN_ind = isnan(Sub_data);
                    Sub_data(logical(NaN_ind)) = 0;
                    N = N + ones(size(Sub_data))-NaN_ind;
                    
                    % Add subject's data to cumulative variable
                    mean_surf = mean_surf + Sub_data;
                    
                end
                % Normalise to get mean:
                mean_surf = mean_surf ./ N;
                
                % Save output
                writematrix(mean_surf,fullfile(Results_Dir,'func',Session,[Space,'_',Atlas,'_',Desc,'.txt']));
                
                disp('Done.')
                disp(' ');
            end
        end
    end
end
                    

