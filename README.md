# Scripts made to excract frames from multiple-walker metadynamics according to rbias values in colvar

### Table of Contents

- ***`extract_from_dzvalues.ipynb`***: this jupyter notebook was made to extract frames in an interactive way.
	- For the scope it was created, frames are filtered accordingly to chosen dz1 and dz2 ranges (corresponding to different minima);
	- For each minimum (or each dz1-dz2 pair) N frames can be selected randomly or by being the ones with smallest or largest rbias. 


- ***`extract_frames_as_cmip` directory***. 
	- I this directory there is a submission script _`1run_sbatch_extract`_ and the associated _python_ code  _`1extract_from_dzvalues_ascmip.py`_. For the interval of dz1 and dz2 Minimum and Maximum values indicated by the sbatch submission script, the python script will create a series of directories, each one with frames within a range of dz1 and dz2 values. 
**Note**: Frames extracted in each directory (or dz1/dz2 range) are the ones with largest rbias in that range.
**Note**: the frame with largest rbias in each directory will be named traj0.nc; the frame with second largest rbias will be named traj1.nc and so on ...

	- 2dist.sh script will run _cpptraj_ to calculate distances for each frame in each directory;

	- 3get_dist.py will extract distance data from each directory and grouping them in _data0.dat_, _data1.dat_ corresponding to all traj0.nc, traj1.nc and so on ..

        - 4plot_2dist.py will finally plot selected data files (for example from traj0.nc) in a 2d plot with dz1 and dz2 values as x and y axes, respectively, and colored by distance.

