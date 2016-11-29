# AMM60_paper

Introduction
============
Repository for development of py (and ipynb) files. Since the notebook files have embedded figures in them they do not work well within version control. I use them as dynamic workspace. When stuff works and is sharable and co-editable it goes into a py file, which is loaded into the notebooks.
So, ipynb files are single ownership, multiple readership.
py files are global read and write.

%load file.py # this is a magic command to load a python file into a notebook cell.

%%writefile file.py # this is a magic command to write a python file with the contents of a notebook cell.


**25 Nov**
----------
To many copies of delta_diagnose function and methods. Removed `IT/AMM60_tools.py` and created symbolic link::

    ln -s FASTNEt/amm60_data_tools.py IT/.

Fixed `internaltidemap_AMM60_paper.ipynb` and `pycnocline_mod_obs_virtual_mooring.ipynb` to diagnose delta with function calls. This required treatment of 3D model data or 1D mooring data.
`FASTNEt/amm60_data_tools.py` has therefore been changed (it tests whether the "profile" data is 1D or 3D in space). In particular a **new output field** is added which will break stuff. My stuff is all fixed to work with this.


---
**26 Nov**

* Check `pycnocline_mod_obs_virtual_mooring.ipynb` works
* start edit `diagnose_delta from amm60_data_tools.py`, for running mean
* Implement running mean processing in `diagnose_delta`
* Implement data clipping in `pycnocline_mod_obs_virtual_mooring.ipynb` for pycnocline diagnostics that were computed outside the range that the Doodson filter gave values for no-tide. Note, however, since the running window is computer over +/- 1.5 day window, the first 1.5 days from the first delta_nt will not have from 1.5 increasing to the 3 days of data in the analysis.

---
**27 Nov**

* Relabel `pycnocline_mod_obs_virtual_mooring.ipynb` --> `pycnocline_virtual_mooring.ipynb`
* Create paper ready figure for by pycnocline depth and std diagnostics map onto observational T(z,t).
* Create paper ready figure for intercomparison between pycnocline depth and std for all configuration+obs at all FASTNEt locations.
* Shift core code to python file: `pycnocline_virtual_mooring.py`

---
**28 Nov**

* Cloned repo on SAN. Testing is this is going to break stuff 
* Seems to work. Syncing from SAN now through ssh tunnels.
* Plotting std(delta) vs time. Created py file

---
**29 Nov**

* Testing fresh clone of repo
* Running mean strat. Update py file

---


**To-do**

`internaltidemap_AMM60_paper.ipynb`:
* strat does not vary with time
* grey mask issue not working

`diagnose_delta from amm60_data_tools.py`:
KG: make sure it still works for your FFT code
