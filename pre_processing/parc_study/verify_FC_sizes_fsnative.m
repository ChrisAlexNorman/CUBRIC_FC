%% SETTINGS & PATHS
Results_dir = '/scratch/scw1648/proj_cn/data/mica_processed';
Save_fil = '/scratch/scw1648/proj_cn/micapipe/parcellation_investigation/parcellation_sizes_fsnative.mat';
FileIndicators = importdata([Results_dir,'/../FileIndicators.csv'],',');

IDs = FileIndicators.data(:,1);
sessions = {'rest','movie'};
atlases  = {...
    'aparc-a2009s',...
    'aparc',...
    'economo',...
    'glasser-360',...
    'schaefer-100',...
    'schaefer-200',...
    'schaefer-300',...
    'schaefer-400',...
    'schaefer-500',...
    'schaefer-600',...
    'schaefer-700',...
    'schaefer-800',...
    'schaefer-900',...
    'schaefer-1000',...
    'vosdewael-100',...
    'vosdewael-200',...
    'vosdewael-300',...
    'vosdewael-400',...
    };
parcellations = [...
    150,...
    70,...
    88,...
    360,...
    100:100:1000,...
    100:100:400,...
    ];

func_out_sizes = struct;
func_out_sizes.err_1={};
func_out_sizes.err_2={};
func_out_sizes.err_3={};
func_out_sizes.warnings={};
for IDi = 1:length(IDs); ID = IDs(IDi);
    subject = ['sub-',num2str(ID)];
    for sess_i = 1:length(sessions); session = sessions{sess_i};
        if strcmp(session,'rest') && ~FileIndicators.data(IDi,3)
            % Rest doesn't exist
            continue
        elseif strcmp(session,'movie') && ~FileIndicators.data(IDi,5)
            % Movie doesn't exist
            continue
        elseif strcmp(session,'rest') && ID == 1153
            % This one is broken.
            continue
        end
        
        % Subcortex report
        sctx_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/volumetric/',subject,'_space-rsfmri_desc-singleecho_timeseries_subcortical.txt'];
        sctx = importdata(sctx_file_path,' ');
        sctx_size = size(sctx,2);
        func_out_sizes.(['sub_',num2str(ID)]).(session).subcortex=sctx_size;
        if sctx_size ~= 14
            func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,['Subcortex has ',num2str(sctx_size),' cols, expect 14']}];
        end
        
        % Cerebellum report
        cereb_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/volumetric/',subject,'_space-rsfmri_desc-singleecho_timeseries_cerebellum.txt'];
        cereb = importdata(cereb_file_path,' ');
        cereb_size = size(cereb,2);
        func_out_sizes.(['sub_',num2str(ID)]).(session).cerebellum=cereb_size;
        if ~any(cereb_size == 32:34)
            func_out_sizes.err_3 = [func_out_sizes.err_3;{ID,session,atlas,['Cerebellum has ',num2str(cereb_size),' cols, expect 32 - 34, unknown cause']}];
        end
        
        % Cortex parcellation reports
        func_out_sizes.(['sub_',num2str(ID)]).(session).cortex = {'Atlas','Expected # of parcellations','TS_Size','Err_1: TS_Size~=expected number of parcellations + 50','zeros column exists','Err_2: non-cortex cols dont add to 50?'};
        for atlas_i = 1:length(atlases); atlas = atlases{atlas_i};
            func_out_sizes.(['sub_',num2str(ID)]).(session).cortex{atlas_i+1,1} = atlas;
            TS_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/surfaces/',subject,'_rsfmri_space-fsnative_atlas-',atlas,'_desc-timeseries.txt'];
            TS_mat = dlmread(TS_file_path,' ');
            
            % Print expected number of parcellations:
            func_out_sizes.(['sub_',num2str(ID)]).(session).cortex{atlas_i+1,2} = parcellations(atlas_i);
            
            % Print full size of TS
            func_out_sizes.(['sub_',num2str(ID)]).(session).cortex{atlas_i+1,3} = size(TS_mat,2);
            
            % Toggle error 1 if full size of TS ~= expected # of parcellations + 50
            err_1 = (size(TS_mat,2) ~= (parcellations(atlas_i) + 50));
            func_out_sizes.(['sub_',num2str(ID)]).(session).cortex{atlas_i+1,4} = err_1;
            if err_1
                % Expected exemptions:
                if strcmp(atlas,'aparc-a2009s') && size(TS_mat,2) == (parcellations(atlas_i) + 48)
                    func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' has #cols = #parcellations + 48 (not 50)']}];
                elseif strcmp(atlas,'economo') && size(TS_mat,2) == (parcellations(atlas_i) + 48)
                    func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' has #cols = #parcellations + 48 (not 50)']}];
                else
                % Otherwise report error
                    func_out_sizes.err_1 = [func_out_sizes.err_1;{ID,session,atlas}];
                end
            end
            
            % Check if all elements of processed cereb are close to integers
            cereb_int_tog = sum(sum(abs(round(TS_mat(:,(sctx_size+1):(sctx_size+cereb_size)))-TS_mat(:,(sctx_size+1):(sctx_size+cereb_size))))) < numel(TS_mat(:,(sctx_size+1):(sctx_size+cereb_size)))*5e-2;
            % Check if the column afterwards is also close to integers
            extra_int_col_tog = sum(abs(round(TS_mat(:,sctx_size+cereb_size+1))-TS_mat(:,sctx_size+cereb_size+1))) < numel(TS_mat(:,sctx_size+cereb_size+1))*5e-2;
            % In a few cases there are two extra columns added
            dbl_xtra_col_tog = sum(abs(round(TS_mat(:,sctx_size+cereb_size+2))-TS_mat(:,sctx_size+cereb_size+2))) < numel(TS_mat(:,sctx_size+cereb_size+2))*5e-2;
            % Infer consequences
            added_col_tog = (cereb_int_tog && extra_int_col_tog) + dbl_xtra_col_tog;
            func_out_sizes.(['sub_',num2str(ID)]).(session).cortex{atlas_i+1,5} = added_col_tog;
            if added_col_tog == 1
                 func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' has ONE added cerebellum column']}];
            elseif added_col_tog == 2
                 func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' has TWO added cerebellum columns']}];
            elseif cereb_size == 33 || cereb_size == 32
                func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' does not appear to have added cerebellum columns (but probably should)']}];
            end
                
            % Toggle error 2 if non-cortical columns do not add to 50
            err_2 = ((sctx_size + cereb_size + added_col_tog + 2) ~= 50);
            if err_2
                % Expected exemptions:
                if strcmp(atlas,'aparc-a2009s') && size(TS_mat,2) == (parcellations(atlas_i) + 48) && ((sctx_size + cereb_size + added_col_tog + 2) == 48)
                    err_2 = 0;
                    func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' has 48 non-cortical cols']}];
                elseif strcmp(atlas,'economo') && size(TS_mat,2) == (parcellations(atlas_i) + 48) && ((sctx_size + cereb_size + added_col_tog + 2) == 48)
                    func_out_sizes.warnings = [func_out_sizes.warnings;{subject,session,[atlas,' has 48 non-cortical cols']}];
                    err_2 = 0;
                elseif cereb_size == 32
                    % Repeat 
                % Otherwise report error
                    func_out_sizes.err_2 = [func_out_sizes.err_2;{ID,session,atlas}];
                end
            end
            func_out_sizes.(['sub_',num2str(ID)]).(session).cortex{atlas_i+1,6} = err_2;
        end
    end
end

func_out_sizes.Notes.err_1 = 'Error flag for when indicated timeseries does not have #columns = #parcellations + 50 (exceptions generate warning)';
func_out_sizes.Notes.err_2 = 'Error flag for when indicated timeseries does not have 48 non-cortical columns (subcortex + cerebellum + potentially an added column or two)';
func_out_sizes.Notes.err_3 = 'Error flag for when cerebellum has number of columns not in range of 32 - 34 as expected';
func_out_sizes.Notes.warnings{1} = 'Expected errors (with description):';
func_out_sizes.Notes.warnings{2} = 'Note: if Cerebellum does not have 34 parcellations warning expects to be followed by 18 warnings (one for each parcellation) identifying ONE or TWO added columns.';

save(Save_fil,'func_out_sizes');
