# Scripts made to excract frames from multiple-walker metadynamics according to rbias values in colvar

### Table of Contents
- ***extract_frames directory***: in this directory there is a submission script "**1run_sbatch_extract**" and "**1extract_from_dzvalues_formwf.py**".
For the interval of dz1 and dz2 Minimum and Maximum values indicated by the sbatch submission script, the python script will create a series of directories, each one with frames within a range of dz1 and dz2 values. Frames exctracted in each directory (or dz1/dz2 range) are the ones with largest rbias in that range.
- ***extract_from_dzvalues.ipynb***: this jupyter notebook was made to extract frames in an interactive way.
For the scope it was created, frames are filtered accordingly to chosen dz1 and dz2 ranges (corresponding to different minima);
For each minimum (or each dz1-dz2 pair) N frames can be selected randomly or by being the ones with smallest or largest rbias. 



