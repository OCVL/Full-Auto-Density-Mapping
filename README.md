# Full-Auto-Density-Mapping
The reference implementation of the algorithms described in:

Cooper, R.F., Kalaparambath, S., Aguirre, G.K., Morgan, J.I.W. "The normative human photoreceptor mosaic: Analysis of a publicly available adaptive optics image repository"

This repository implements a fully automated approach for estimating density in a large set montages. It not only assesses full montages, it is also capable of merging confocal, and split-detector modalities into a single dataset based on its "confidence" in its result for a given location. **To recreate the results from the paper, the code must be run in the following sequence:**

1. [Multi_Montage_DFT_Analysis.m](#run-multi_montage_dft_analysism)
2. [DFT_Aggregator.m](#run-dft_aggregatorm)
3. [Aggregate_Analyses.m](#run-aggregate_analysesm)

**Note: This code uses image processing, curve fitting, and parfor functions, which requires the Image Processing, Curve Fitting, and Parallel Toolboxes from MATLAB.**

## Run Multi_Montage_DFT_Analysis.m:

Running this script will ask you to:
  1. Select what the output unit should be. At present, the options are:
     - Microns (using mm^2 for density)
     - Degrees
  2. Select the root folder containing folders of montaged images that you wish to analyze. The script at present expects a folder structure mirroring that used for the paper:
     ```
     -Root Folder/Subject/confocal
     -Root Folder/Subject/split detection
     ```
     - **Note:** Each montage image should be a single image placed within an blank image the size of the full montage.
     - If you only have a PSD of a montage and do not have data in this format, please use [the PSD Layer Exporter](https://github.com/Eurybiadan/PSD_Layer_Export/releases/tag/v1.0). To dump them to disk in this format.     
  3. Select the scaling lookup table. The lookup table allows the software to analyze a folder of images from different subjects/timepoints/conditions. The lookup table itself **must** be a 3 column 'csv' file, where:
     - The **first column** is a common identifier for image/coordinate pairs
     - The **second column** is the axial length (or '24' if the axial length is unknown) of the image/coordinate pairs
     - The **third column** is the pixels per degree of the image/coordinate pairs. Each row must contain a different identifier/axial length/pixels per degree tuple.

     An example common identifier could be a subject number, e.g, when working with the files
     - 1235_dateoftheyear_OD_0004.tif
     - 1235_dateoftheyear_OD_0005.tif
      
     Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD". If all three were placed in a LUT, then the one that matches the most (as determined via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".
      
     If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then _only_ the identifier "1235" would match between all images. However, say the two dates have different scales, then you would want to create two rows in the look up table for each date, with identifiers like: "1235_dateoftheyear" and "1235_differentdateoftheyear".

  4. (**Optional**) Select a foveal image list. This image list simply contains a list (column) of images in a csv that were determined by the user (you) to be of particularly good quality at fovea. These images will be analyzed separately from the other confocal images, and given priority when all data is merged at a later step. This is not required, but we found that it enabled us to achieve better quality at the fovea.

After answering the above prompts, the script will then run, saving the results to disk as a mat file in an "All Analyses" folder in the "Root Folder" specified previously you ran from alongside a mat file that contains the results.


## Run DFT_Aggregator.m

After Multi_Montage_DFT_Analysis.m completes, the next step is to aggregate (or combine) all montages together. This step is comparably simple.

Running this script will ask you to:
  1. Select all montage mat files you wish to combine. As these were all automatically placed into the "All Analyses" folder by the prior step, navigate to that folder, and Ctrl (or Command) -click on all montages you wish to aggregate.
  2. The software will run, and output the aggregated data to a "Results" folder within the All Analyses folder selected.

## Run Aggregate_Analyses.m

This script produces the analyses used for the figures in the paper (density/confidence/total cone plots). Like the above, this step is comparably simple.

Running this script will ask you to:
  1. Select the aggregate file. This script assumes that all of the individual files are one folder above where the aggregate is.
  2. The software will run, and output the analysis figures to the "Results" folder.

