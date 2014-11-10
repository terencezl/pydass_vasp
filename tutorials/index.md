---
layout: page
title: Tutorials
---

Here let us get familiar with the sub-packages/modules in `pydass_vasp`. Tutorials are written in the form of 
IPython notebooks, which you can download when you are in the following sections.

These sections correspond to the respective sub-packages.
 
1. [plotting](plotting/) ([download]({{ site.baseurl }}/notebooks/plotting.ipynb))
2. [fitting](fitting/) ([download]({{ site.baseurl }}/notebooks/fitting.ipynb))
3. [manipulation](manipulation/) ([download]({{ site.baseurl }}/notebooks/manipulation.ipynb))
4. [xml_utils](xml_utils/) ([download]({{ site.baseurl }}/notebooks/xml_utils.ipynb))

Before you start, there are some general configurations for the IPython interpreter, including some imports and 
matplotlib customization. Have a look at my config file `ipython.py`. Just copy it all, or make your own changes, and
place it under `~/.ipython/profile_default/startup/` for Mac OS X and Linux, and next time you go into the interpreter,
the lines will be executed. For Windows, in command line prompt type in `ipython profile locate`.
This would be the equivalent part of `~/.ipython/profile_default` for the location in Mac OS X and Linux.

Here is my configuration file.

```python
import os
import sys
try:
    # import NumPy
    import numpy as np
    print("NumPy is imported.")
    np.set_printoptions(suppress=True)
except ImportError:
    print("NumPy is not imported!")

try:
    # import matplotlib
    import matplotlib as mpl
    print("matplotlib is imported.")
    # MacOSX backend has some problems, so change it to TkAgg
    if mpl.get_backend() == 'MacOSX':
        mpl.use('tkagg')
    import matplotlib.pyplot as plt
    print("Using TkAgg backend.")
    # turn on interactive mode
    plt.ion()
    print("Using matplotlib interactive mode.")
    # try using ggplot style for prettier looks
    try:
        plt.style.use('ggplot')
        print("Using ggplot style from matplotlib 1.4.")
    except ValueError:
        print("If matplotlib >= 1.4 is installed, ggplot style will be used for better looks.")

    # figure size and resolution adjustment
    mpl.rcParams['figure.figsize'] = [5.5,4]
    mpl.rcParams['figure.dpi'] = 80
    mpl.rcParams['savefig.dpi'] = 100
except ImportError:
    print("matplotlib is not imported!")
```
