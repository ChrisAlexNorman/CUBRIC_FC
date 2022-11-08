# CN ALSPAC Data Analysis Project

**Author:** Chris Norman

**Contact:** NormanC4@cardiff.ac.uk

This document provides an overview of the methods applied by the author for analysis of the ALSPAC dataset. It is presented, roughly, in chronological order, and may be considered a reference guide to the more detailed, low-level reports and README files.

## Data Processing

Structural scans and functional resting state and movie scans were selected for processing for each subject in *(SCW_PATH)/full_dataset/(SUBJECT)/Connectom/scans*. Missing scans were noted. Incomplete scans (those with very few DICOM files) were omitted. In cases where a scan was repeated then the scan with the largest number of DICOM files was selected, if multiple scans had the same number of DICOM files then the latest recording session was selected. The file names for each subject's selected scans are available in *data/FileLinks_full_dataset.csv*.

The [dcm2niix](https://github.com/rordenlab/dcm2niix) tool was used to converted selected DICOM files to NIfTI format by *utils/get_nifti.sh*. These unprocessed .nii files are stored in *data /NIfTI_Data*. The filesystem and naming conventions are, as far as possible, in line with [BIDS specification](https://bids-specification.readthedocs.io/en/stable/). *data/FileIndicators.csv* uses a binary indicator for whether each subject's scans are available in NIfTI form.

Reverse phase encoding resting state scan were used, where available, to unwarp resting state and movie functional scans for each subject. This step was performed by Hikaru Tsujimura (tsujimuraH@cardiff.ac.uk) using the [topup](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide/) utility. Unwarped data was moved to alongside the original, unwarped NIfTI files.

The [micapipe toolbox](https://micapipe.readthedocs.io/en/latest/pages/01.whatyouneed/index.html) was used to process structural and functional NIfTI data using established pipelines, see *micapipe/README.md*. The outputs are stored in *data/mica_processed*.
