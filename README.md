pydass_vasp (or, badass wasp)
=============================
Convenient Python modules and wrapping script executables.

##### Example: plotting band structure

```python
pydass_vasp.plotting.plot_bs(axis_range=[-4,6])
```

![band_structure](http://terencezl.github.io/pydass_vasp/images/band_structure.png)

Returned dictionary:

```python
{'ax': <matplotlib.axes._subplots.AxesSubplot at 0x10a33df50>,
 'data': {'columns': ['k_points', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', 
                        '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32'],
          'data': array(
            [[  0.        , -20.342219  , -16.616756  , ...,   5.849101  ,
                5.855091  ,   6.074841  ],
             [  0.04558028, -20.342181  , -16.616823  , ...,   5.811826  ,
                5.815311  ,   6.060851  ],
             [  0.09116057, -20.34223   , -16.617067  , ...,   5.730248  ,
                5.734556  ,   5.80481   ],
             ..., 
             [  2.49869989, -20.343194  , -16.628521  , ...,   5.172637  ,
                5.204402  ,   5.711173  ],
             [  2.53591604, -20.343228  , -16.6286    , ...,   5.219897  ,
                5.226956  ,   5.730676  ],
             [  2.57313218, -20.34319   , -16.628622  , ...,   5.234177  ,
                5.234205  ,   5.726715  ]])},
 'reciprocal_point_locations': array([ 0.        ,  0.8660254 ,  1.3660254 ,  1.8660254 ,  2.57313218]),
 'reciprocal_points': ['R', 'G', 'X', 'M', 'G']}
```

##### Example: plotting total density of states

```python
pydass_vasp.plotting.plot_tdos()
```

![density_of_states](http://terencezl.github.io/pydass_vasp/images/dos.png)

Returned dictionary:

```python
{'ax': <matplotlib.axes._subplots.AxesSubplot at 0x10acb27d0>,
 'data': {'columns': ['E', 'total', 'integrated'],
          'data': array(
            [[-27.51051771,   0.        ,   0.        ],
             [-27.45851771,   0.        ,   0.        ],
             [-27.40551771,   0.        ,   0.        ],
             ..., 
             [ 14.41548229,   0.        ,  16.        ],
             [ 14.46748229,   0.        ,  16.        ],
             [ 14.52048229,   0.        ,  16.        ]])}}
```

What is it?
===========

This Python package is the result of frustration of searching for an organized, straightforward and flexible approach of plotting, fitting and manipulation of [VASP](https://www.vasp.at/) files (typically `POSCAR`) in a few key strokes as long as you have a terminal, and preferably a (S)FTP client. It has the following features:

* It is a Python package with a straightforward structure. When you `import pydass_vasp`, you have sub-packages `pydass_vasp.plotting`, `pydass_vasp.fitting`, `pydass_vasp.manipulation`, `pydass_vasp.xml_utils`, each containing a few functions to carry out your tasks, with a careful selection of options to choose from. Return values are Python dictionaries, informative enough to be flexible for post-processing. When there is no need for an object, don't bother creating one.

* It also has scripts that utilize the main package. They are simply 'wrappers' to the functions that are used the most frequently, typically `plotting` and `fitting` functions. Instead of `cd`ing into a directory and launching the Python interpreter `python`, making the imports and calling the functions, you just need to stay in your terminal shell and type `plot_tdos.py` to plot the total density of states (DOS), `plot_ldos.py 1` to plot the local projected DOS of the first atom, and `plot_bs.py` to plot the band structure. These scripts accept arguments and options as flexible as their function counterparts.

* The defaults of the functions and wrapping scripts are sane and cater to the most use cases. For example, If you are in a Python interpreter, simply typing in `pydass_vasp.plotting.plot_tdos()` will collect data to plot from `DOSCAR` and parameters necessary under the current VASP job directory, generate and **display** a figure (or more) but **not save** it/them on disk, and at the same time return a dictionary with the extracted data and the [matplotlib](http://matplotlib.org/) axes handler references. It'll automatically obtain critical parameters such as `ISPIN`, `E-fermi` first from VASP **output** files (e.g. `OUTCAR`), then **input** files (e.g. `INCAR`) if the first attempt fails, and decide the number of figures to generate. If you are just in a terminal shell, typing in the script `plot_tdos.py` will do the same, but **rather than display** the figure(s), instead, **save** it/them quietly and use matplotlib's `Agg` backend. This would be particularly helpful for terminal users who don't have X Window forwarding ([Xming](http://www.straightrunning.com/XmingNotes/) for Windows, [XQuartz](http://xquartz.macosforge.org/landing/) for Mac OS X) set up on their own local machine, or the forwarding connection is slow to hold the live generated figures.
 
* The options to the functions and wrapping scripts provide you with room of customization from the beginning. As an example, we consider `pydass_vasp.plotting.plot_tdos()`, shorted as `plot_tdos()`.
	* `plot_tdos(input_file='vasprun.xml')` switches from taking in `DOSCAR` to `vasprun.xml`. It lets you select what file you prefer to use. Any filename containing `'DOSCAR'` is considered to be of `DOSCAR` type, any filename ending with `'.xml'` is considered to be of `vasprun.xml` type.
	* `plot_tdos(axis_range=[-15,5,0,10])` sets the custom axis range.
	* `plot_tdos(ISPIN=2)` lets you manually override the auto-detection of `ISPIN` from files other than `DOSCAR`. The program will skip the corresponding part of work. This is helpful when you only have the major data file `DOSCAR` transferred to you local machine, and do not have the other files necessary to extract the parameters to proceed plotting. To leave no confusion, when `ISPIN` is 2, two figures are generated, one with spin up and down combined, the other with two overlapping curves, denoting spin up and spin down separately.
	* `plot_tdos(on_figs=1)` creates the plot on top of an existing matplotlib figure labeled as `Figure 1`, instead of generating a whole new one.
	* When `ISPIN` is 2, `plot_tdos(on_figs=[1,2])` puts the two plots mentioned before onto `Figure 1` and `Figure 2`. `plot_tdos(on_figs=[1,None])` is also valid, meaning putting the combined curve to `Figure 1`, and the two overlapping curves to a new figure, which you can of course delete on its own.
	* `plot_tdos(display=False, close_figs=True, save_figs=True)` replicates the behavior of the corresponding wrapping script `plot_tdos.py`. There is something more to say about the interactions between the arguments `display` and `close_figs`[^1].
	* `plot_tdos(save_data=True, output_prefix='TDOS')` saves the extracted data to disk, with the prefix `'TDOS'`. The argument `output_prefix` also specifies the filenames for saved figures.
	* The wrapping script `plot_tdos.py` accepts the relevant options in the form of `-i vasprun.xml`, or `--input vasprun.xml`, `-p` or `--display`. For more readily available information, type in `plot_tdos.py -h` to get help direclty from the terminal shell.

* The returned dictionary also leave room for adjustments. Call it

    ```python
	returned_dict = {		
		'data': {'columns': col_names, 'data': data}
		'axes': {'ax': ax}
	}
    ```
		
	`returned_dict['data']` has a 2D numpy array of data, and their column names. This construction is prefered because if you have [pandas](http://pandas.pydata.org/), you can just convert it to a DataFrame by `pd.DataFrame(**returned_dict['data'])`.

	`returned_dict['axes']` contains the matplotlib axes handlers. When `ISPIN` is 2, it has two elements: `'ax_spin_up'` and `'ax_spin_down'`.

* It has a uniform plotting support for the Crystal Orbital Hamilton Populations (COHP) analysis tool [LOBSTER](http://cohp.de/), function `plot_cohp()` and script `plot_cohp.py`.

* If you use matplotlib>=1.4, and you plot with the wrapping scripts, you can automatically enjoy the aesthetics of [ggplot](http://tonysyu.github.io/mpltools/auto_examples/style/plot_ggplot.html) by its newly added sub-package [style](http://matplotlib.org/1.4.2/users/whats_new.html#style-package-added).
	
[^1]: Note, if you are in a python interpreter and matplotlib **interactive** mode (`matplotlib.is_interactive() == True`), the plots will show up anyway regardless of the `display` argument. You can of course close them later. But if you really do not need the plots, and just want to have the convenient access to the extracted data in the returned dictionary, type in `plot_tdos(close_figs=True)` instead, which closes the figures and returns no axes handlers either. The argument `close_figs` specifies the closing behavior when the figures don't show up (`display=False`), or **after** they show up (`display=True`) in **interactive** mode. In **non-interactive** mode, this argument has no effect because the figure blocks/freezes the process, and when you click to close the figure, it is closed anyways. The logical arrangement should be reasonable and easy to comprehend.

Dependencies
============
* Python 2.7 (Python 3 support is currently not considered)
* NumPy
* SciPy
* matplotlib
* IPython (optional, but better to have)

I highly recommend every scientist/researcher who is new to Python to install the scientific superpack [Anaconda](https://store.continuum.io/cshop/anaconda/), if you are using Windows or Mac OS X. Even if you are on Linux, it is still highly recommended if you don't have superuser control over the machine to install packages freely. It is often the case when you have ssh access to a supercomputer. In all these cases, just download the package and do a simple local installation, and you already have everything to start with, not only for `pydass_vasp`, but also for the whole adventure of scientific computing.

Installation
============
~~This package has already been registered on PyPI.~~
So if you have [pip](https://pip.readthedocs.org/en/latest/), which is a must have, and should already have been included in Anaconda,

~~pip install pydass_vasp~~
	
Or if you wish to follow the more updated releases, which should serve you better because small projects can fully enjoy the freedom of updates on GitHub before committing to PyPI,

	pip install git+https://github.com/terencezl/pydass_vasp
	
Alternatively, if you don't have pip, download the `.tar.gz` file from the PyPI page and decompress it in a local directory of your choice, or `git clone https://github.com/terencezl/pydass_vasp`, get into the outer `pydass_vasp` directory and

	python setup.py install
	
	# then you can get out of that directory and just delete it
	cd ..
	rm -r pydass_vasp
	
However, `pip` installation is always recommended, because of the ease of uninstallation,

	pip uninstall pydass_vasp

Getting Help
============
An organized tutorial/documentation has not been ready yet, but the docstrings of functions are fairly complete. If you use [IPython](http://ipython.org/), which is also a must have, and should already have been included in Anaconda,

	# in terminal shell (e.g. Bash)
	$ ipython
	
	# in ipython
	>>> import pydass_vasp
	>>> pydass_vasp.plotting.plot_tdos?
	
You will get help. Experiment on a few options and you'll be quickly on your way.

In addition, the help texts of scripts (option `-h` when you call them in the terminal shell) are available as well.
