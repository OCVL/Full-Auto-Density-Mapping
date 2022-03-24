# Full-Auto-Density-Mapping
A repository that implements a fully automated approach for estimating density, and spacing both in single images and in montages.


### Run Montage_DFT_Analysis.m:

Running this function will ask you to select the set of montage images that you wish to analyze. Montage images will have a single layer within a montage placed within an image the size of the montage canvas.

**If you only have a PSD of a montage and do not have data in this format, please use [my PSD Layer Exporter](https://github.com/Eurybiadan/PSD_Layer_Export/releases/tag/v1.0). To dump them to disk in this format.

**This function uses parfor loops, which requires the Parallel Toolbox from MATLAB.** If you don't have that toolbox, change the "parfor" on line 93 to a "for" loop.

As above, it will then prompt the user to select what the output unit should be. At present, the options are:
* Microns (using millimeters^2 for density)
* Degrees
* Arcminutes

Once the output unit is select, it will give the user the option to pick a lookup table. The lookup table allows the software to analyze a folder of images from different subjects/timepoints/conditions. The lookup table itself **must** be a 3 column 'csv' file, where the **first column** is a common identifier for image/coordinate pairs, the **second column** is the axial length (or '24' if the axial length is unknown) of the image/coordinate pairs, and the **third column** is the pixels per degree of the image/coordinate pairs. Each row must contain a different identifier/axial length/pixels per degree tuple.

An example common identifier could be a subject number, e.g, when working with the files
- 1235_dateoftheyear_OD_0004.tif
- 1235_dateoftheyear_OD_0005.tif

Common identifiers could be "1235", "1235_dateoftheyear", "1235_dateoftheyear_OD". If all three were placed in a LUT, then the one that matches the most (as determined via levenshtein distance) will be used. In this example, we would use "1235_dateoftheyear_OD".

If we had another date, say: 1235_differentdateoftheyear_OD_0005.tif, then _only_ the identifier "1235" would match between all images. However, say the two dates have different scales, then you would want to create two rows in the look up table for each date, with identifiers like: "1235_dateoftheyear" and "1235_differentdateoftheyear".

**If you do not wish to use a lookup table, then press "cancel", and the software will allow you put in your own scale in UNITS/pixel.**

The software will then run, showing the spacing montage, the density montage, the confidence montage, and the sum map. It will also save them to disk in the same folder you ran from alongside a mat file that contains the results.
