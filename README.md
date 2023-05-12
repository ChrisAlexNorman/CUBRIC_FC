# CUBRIC_FC

**Author:** Chris Norman

**Contact:** chris.alex.norman@outlook.com

This is a sample of code I produced to process and analyse functional connectivity from the raw brain scans of 220 subjects in the Avon Longitudinal Study of Parents and Children (ALSPAC). Subject data is protected and therefore not provided here however the tools are available for open use.

## Pipeline Overview

Structural scans and functional resting state and movie scans were selected for processing for each subject in *(SCW_PATH)/full_dataset/(SUBJECT)/Connectom/scans*. Missing scans were noted. Incomplete scans (those with very few DICOM files) were omitted. In cases where a scan was repeated then the scan with the largest number of DICOM files was selected, if multiple scans had the same number of DICOM files then the latest recording session was selected. The file names for each subject's selected scans are available in *data/FileLinks_full_dataset.csv*.

The [dcm2niix](https://github.com/rordenlab/dcm2niix) tool was used to converted selected DICOM files to NIfTI format by *utils/get_nifti.sh*. These unprocessed .nii files are stored in *data /NIfTI_Data*. The filesystem and naming conventions are, as far as possible, in line with [BIDS specification](https://bids-specification.readthedocs.io/en/stable/). *data/FileIndicators.csv* uses a binary indicator for whether each subject's scans are available in NIfTI form.

Reverse phase encoding resting state scan were used, where available, to unwarp resting state and movie functional scans for each subject using the [topup](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/TopupUsersGuide/) utility. Unwarped data was moved to alongside the original, unwarped NIfTI files.

The [micapipe toolbox](https://micapipe.readthedocs.io/en/latest/pages/01.whatyouneed/index.html) was used to process structural and functional NIfTI data using established pipelines, see *micapipe/README.md*. The outputs are stored in *data/mica_processed*.

Functional gradients were extracted and analysed alongside other data features using tools in *gradients*.

## Environment ##
| *Software* | *Version*     | *Further info* |
|------------|---------------|----------------|
| dcm2niix   | v1.0.20200331 | https://github.com/rordenlab/dcm2niix |
| micapipe   | 0.1.2         | https://micapipe.readthedocs.io/en/latest/ |
| python     | 3.7.15        | https://www.python.org/downloads/ |

See *requirements.txt* for packages.
