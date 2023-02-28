%% SETTINGS & PATHS
% clear all; close all;
Results_dir = '/scratch/scw1648/proj_cn/data/mica_processed';
Fig_dir = '/scratch/scw1648/proj_cn/micapipe/examples';
save_tog = 0; % Save generated figures to Fig_dir?
inc_sctx_cereb = 1; % Include subcortex and cerebellum with labels

subject = 'sub-1213';
atlas = 'aparc'; % 'glasser-360';
session = 'movie';

FC_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/surfaces/',subject,'_rsfmri_space-fsnative_atlas-',atlas,'_desc-FC.txt'];
TS_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/surfaces/',subject,'_rsfmri_space-fsnative_atlas-',atlas,'_desc-timeseries.txt'];
sctx_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/volumetric/',subject,'_space-rsfmri_desc-singleecho_timeseries_subcortical.txt'];
cereb_file_path = [Results_dir,'/micapipe/',subject,'/func/',session,'/volumetric/',subject,'_space-rsfmri_desc-singleecho_timeseries_cerebellum.txt'];

%% LOAD DATA
FC_mat = dlmread(FC_file_path,' ');
FC_mat = FC_mat + FC_mat' - diag(diag(FC_mat));
TS_mat = dlmread(TS_file_path,' ');
sctx = importdata(sctx_file_path,' ');
cereb = importdata(cereb_file_path,' ');
sctx_cereb = [sctx,cereb];

% Note the difference in cerebellum columns.
% Somewhere during spike removal & smoothing (Lines 342 - 392 in micapipe/functions/03_FC.py)
% Values are converted to integers, and a column of zeros is added...

%% CHECK SIZES
sctx_size = size(sctx,2); disp(['Subcortex: ',num2str(sctx_size),' parcellations'])
cereb_size = size(cereb,2); disp(['Cerebellum: ',num2str(cereb_size),' parcellations'])
if all(TS_mat(:,sctx_size+cereb_size+1)==zeros(size(TS_mat,1),1))
    disp(['Column ',num2str(sctx_size+cereb_size+1),' is zeros'])
    disp(['Remaining ',num2str(size(TS_mat,2)-sctx_size-cereb_size-1),' columns belong to the cortex']);
else
    disp(['!!! Column ',num2str(sctx_size+cereb_size+1),' is NOT zeros !!!'])
    disp(['Do the remaining ',num2str(size(TS_mat,2)-sctx_size-cereb_size),' columns all belong to the cortex?']);
end

%% Make figures
% TS
TS_mat_plot = TS_mat';
figure('Position',[1800 -800 800 600]);
colormap hot
imagesc(TS_mat_plot);
hold on; set(gca,'FontSize',14); title('Timeseries'); xlabel('Scan #'); ylabel('Parcellation');
if inc_sctx_cereb
    plot([0.5,size(TS_mat_plot,2)+0.5],[sctx_size+0.5,sctx_size+0.5],'k-','LineWidth',2);
    plot([0.5,size(TS_mat_plot,2)+0.5],[sctx_size+cereb_size+1.5,sctx_size+cereb_size+1.5],'k-','LineWidth',2);
    text(5,(sctx_size+0.5)/2,'Subcortex','FontSize',14,'FontWeight','bold','Rotation',0)
    text(5,(2*sctx_size+cereb_size+1.5)/2,'Cerebellum','FontSize',14,'FontWeight','bold','Rotation',0)
    text(5,sctx_size+cereb_size+3.5,'Cortex','FontSize',14,'FontWeight','bold','Rotation',0,'VerticalAlignment','top')
end
if save_tog
    saveas(gcf,[Fig_dir,'/',subject,'_space-fsnative_atlas-',atlas,'_',session,'-timeseries.png'])
end

% FC
figure('Position',[1800 -800 800 600]);
colormap hot
imagesc(FC_mat);
hold on; set(gca,'FontSize',14); title('FC')
if inc_sctx_cereb
    plot([sctx_size+0.5,sctx_size+0.5],[0.5,size(TS_mat,2)+0.5],'k-','LineWidth',2);
    plot([sctx_size+cereb_size+1.5,sctx_size+cereb_size+1.5],[0.5,size(TS_mat,2)+0.5],'k-','LineWidth',2);
    plot([0.5,size(TS_mat,2)+0.5],[sctx_size+0.5,sctx_size+0.5],'k-','LineWidth',2);
    plot([0.5,size(TS_mat,2)+0.5],[sctx_size+cereb_size+1.5,sctx_size+cereb_size+1.5],'k-','LineWidth',2);
    text((sctx_size+0.5)/2,0.99*size(TS_mat,2),'Subcortex','FontSize',14,'FontWeight','bold','Rotation',90)
    text((2*sctx_size+cereb_size+1.5)/2,0.99*size(TS_mat,2),'Cerebellum','FontSize',14,'FontWeight','bold','Rotation',90)
    text(sctx_size+cereb_size+3.5,0.98*size(TS_mat,2),'Cortex','FontSize',14,'FontWeight','bold','Rotation',0,'VerticalAlignment','middle')
end
if save_tog
    saveas(gcf,[Fig_dir,'/',subject,'_space-fsnative_atlas-',atlas,'_',session,'-FC.png'])
end
%% Quick check
% sctx_cereb is timeseries data BEFORE spike removal & smoothing
% TS_mat is timeseries (with cortex appended) AFTER spike removal & smoothing
% This just checks if the first element of each is equal.
Check_eq = zeros(1,size(sctx_cereb,2));
for j = 1:size(sctx_cereb,2)
    if sctx_cereb(1,j) == TS_mat(1,j)
        Check_eq(j) = 1;
    end
end

%% Replication of excerpt from micapipe/functions/03_FC.py
% Not currently in use

[Vertices_lh, Labels_lh, Colortable_lh] = read_annotation(...
    [Results_dir,'/freesurfer/',subject,'/label/lh.',atlas,'_mics.annot']);
[Vertices_rh, Labels_rh, Colortable_rh] = read_annotation(...
    [Results_dir,'/freesurfer/',subject,'/label/rh.',atlas,'_mics.annot']);

native_parc = zeros(length(Labels_lh)+length(Labels_rh),1);
for x = 1:length(Labels_lh)
    native_parc(x) = find(Colortable_lh.table(:,5)==Labels_lh(x));
end
for x = 1:length(Labels_rh)
    native_parc(x+length(Labels_lh)) = find(Colortable_rh.table(:,5)==Labels_rh(x)) + size(Colortable_lh.table,1);
end
uparcel = unique(native_parc);

% This is the key section which generates a parcellated timeseries from a
% (cleaned) vertex-wise timeseries.
% It will not work because the mapping from vertices to parcellation labels
% (.annot files above) are registered to fsnative space while the cleaned
% vertex-wise timeseries file is registered to conte-69 space...
% 
% data_corr = importdata('/scratch/scw1648/proj_cn/EXAMPLE_mica_outputs/micapipe/sub-1011/func/rest/surfaces/sub-1011_rsfmri_space-conte69-32k_desc-timeseries_clean.txt',' ');
% ts_native_ctx = zeros(size(data_corr,1),length(uparcel));
% for lab = 1:length(uparcel)
%     tmpData = data_corr(:,native_parc == uparcel(lab));
%     ts_native_ctx(:,lab) = mean(tmpData,1);
% end        
