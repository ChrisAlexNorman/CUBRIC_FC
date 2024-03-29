# Pre-processing Pipeline 

**Author:** Chris Norman

**Contact:** chris.alex.norman@outlook.com

**micapipe version:** 0.1.2

Structural and functional connectome pre-processing on the ALSPAC dataset was performed using the micapipe pipeline. For simplicity, strucutral and functional pipelines were applied separately and sequentially. The default pipelines, [see documentation](https://micapipe.readthedocs.io/en/latest/pages/01.whatyouneed/index.html), were used, with some notable deviations recorded below.

Check the main scripts *pipeline_struct/micapipe_struct_sub.sh* and *pipeline_func/micapipie_func_sub.sh* to see how the pipeline was adapted to be run on the SCW server on the ALSPAC dataset.

## File Outlines

*pipeline_struct* and *pipeline_func* contain the scripts related to structural and functional processing respectively. Completion of both pipelines was (partially) verified using *pipeline_func/verify_FC_file_sizes.sh*. This returns the size of functional directories output by micapipe which are expected to be populated. The output in *pipeline_func/logs/func_dir_sizes.csv* demonstrates that subject 1153's resting state functional processing did not complete. *pipeline_func* also contains two functions for visualising functional connectivity (FC) and timeseries (TS) output matricies.

*examples* contains raw NIfTI data, processed micapipe outputs, and illustrative .png files, from subjects 1011 and 1028.

*clean_micapipe.sh* calls the [micapipe_cleanup](https://micapipe.readthedocs.io/en/latest/pages/05.micapipe_cleanup/index.html?highlight=micapipe_cleanup#micapipe-cleanup) function with appropriate file mounts and settings for the SCW environment. Note that micapipe_cleanup is sometimes not found on the path using 'singularity exec' (see Known Issues below).

*parc_study* contains the results of an investigation into how micapipe handles parcellation of functional activity in the cortex since this is not clear from the [micapipe documentation](https://micapipe.readthedocs.io/en/latest/pages/01.whatyouneed/index.html) or [source code](https://github.com/MICA-MNI/micapipe).

*get_fsaverage5_FC.py* is an adaptation of 03_FC.py from micapipe which generates cleaned and parcellated timeseries and FC matrices in fsaverage5 space (and the given parcellation).

*get_fsaverage5_FCs.sh* runs *get_fsaverage5_FC.py* for each subject and session in the dataset.

## Deviations From Default Micapipe Pipeline

The default output filesystem is adapted to account for the different recording modalities in the ALSPAC data, resting and movie watching. Subdirectories 'rest' and 'movie' are created in (OUTPUT_DIRECTORY)/micapipe/(SUBJECT)/func and (OUTPUT_DIRECTORY)/micapipe/(SUBJECT)/xfms, and files generated by micapipe's proc_rsfmri pipeline are directed there. **Note:** It is possible to instead specify a particular 'session' through micapipe, however this leads to duplication of structual outputs for the two functional modalities. At ~1GB per subject, and with 220 subjects, it was considered better to avoid this data replication.

Sidecar .json files associated with each scan are found in the same location as the corresponding .nii files. However, micapipe always assumes that .json information is stored at (OUTPUT_DIRECTORY)/task-rest_bold.json, regardless of subject or session. Therefore a dummy file is created, *touch (OUTPUT_DIRECTORY)/task-rest_bold.json*, and that path is re-directed to the correct .json file for each subject and session.

## Known Issues

Subject 1153's resting state functional processing fails for an unkown reason (check logs?). The failure occurs when run as part of the complete batch job, or if run in isolation. This is the same subject where r/w permissions were previously blocked to users except the owner. Interestingly their 'movie' session is processed without error.

**Re-running subjects:** Micapipe's --force option appears to remove **ALL** previously generated subject outputs, not just overwrite those for the selected module. Therefore it is currently neccessary to re-run **ALL** modules from the beginning. The subfunction *micapipe_cleanup* should be used first to avoid errors but it isn't found on the path using 'singularity exec'. Furthermore, manually removing a subject's output directories still results in errors upon re-running (perhaps because of the content of the (OUTPUT_DIRECTORY)/micapipe/micapipe_processed_sub.csv file?). **Working Solution:** Re-run all a subject's modules, output to a new, temporary directory, then copy the results into main output directory. Note that this bypasses addition to micapipe_processed_sub.csv log file, which may then appear incomplete.

*pipeline_func/verify_FC_file_sizes.sh* provides a quick way to check if a set of scans failed processing, but there is not yet any way of verifying if **all** expected file outputs are generated.
