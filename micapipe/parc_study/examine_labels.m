lh_annot_path = './lh.glasser-360_mics.annot'; %'./lh.HCP-MMP1.annot'; % './lh.my_L.annot'; %'./lh.HCP-MMP1_Pipeline_Codes.annot'; % './lh.my_L.annot';
rh_annot_path = './rh.glasser-360_mics.annot'; % './rh.my_R.annot'; %'./rh.HCP-MMP1_Pipeline_Codes.annot'; % './rh.my_R.annot';

micapipe_labeling_path = '/scratch/scw1648/proj_cn/micapipe/parcellation_investigation/HCP_MMP1/micapipe_parcellations/glasser-360_conte69.csv';

[Vertices_lh, Labels_lh, Colortable_lh] = read_annotation(...
    lh_annot_path);
[Vertices_rh, Labels_rh, Colortable_rh] = read_annotation(...
    rh_annot_path);
my_labeling = [Labels_lh;Labels_rh];
my_labeling_edit = my_labeling;
my_labeling_edit(my_labeling_edit==16777215)=0;
my_lab_0s = find(my_labeling_edit==0);

mica_labeling = csvread(micapipe_labeling_path);
mica_lab_0s = find(mica_labeling==0);

equalities = (my_lab_0s == mica_lab_0s);
all_eq = all(equalities);
