# Installation

### Environment variables
``PASTIS`` requires two environment variables to be set:

- PASTISPATH

- PASTISLIB

### External models
External models are used for the stellar atmospheres, evolution tracks, limb darkening parameters, and interstellar extinction. Besides, a library of photmetric filters are necessary to run PASTIS. All these files are provided as part of the repository in the external_models folder.

**Update** The zipped external models can be found in this [Google Drive](https://drive.google.com/file/d/1oWDlJ45MMGxpTGAqsNxvdBGtnZnjOHiq/view?usp=sharing).

To unpack the library, use ``tar``
```
tar -xzvf pastislib.tgz
```

They have to be placed in the location specified by $PASTISLIB.

### JKTEBOP
``PASTIS`` relies on some subroutines from the ``JKTEBOP`` package by Southworth [ADD LINK].
To make them available to ``PASTIS``, an extension module must be built using ``f2py``:
```
f2py -c task2_v28_components.f -m task2_components
```

If you are running python3, this should be
```
f2py3 -c task2_v28_components.f -m task2_components3
```

In MacOS, you will probably need the command line tools, if they are not installed yet. To do this, run from a terminal:
```
xcode-select --install
```
and follow the instructions.
You also need the ```Python.h''' header file.

The resulting file, ``task2_components.so`` must be located in the directory ``$PASTISLIB/fortran``, which is included in the ``sys.path`` variable at startup.
The original fortran file is found in the ``external_models/fortran`` directory.

### C++
For performance issues the code to fit CCFs using a Gaussian function is written in C++. This is provided in the ``$PASTISLIB/cpp`` directory. The shared library used by ``PASTIS`` is compiled by running

```
icc -O2 -fPIC -shared fitgauss_c.cpp -o fitgauss_c.so
```

If you do not have the (much faster propietary) ``icc`` compiler, you can use ``g++``.

# Reference

If you use this code in a publication, please cite our papers:

```
@ARTICLE{2014MNRAS.441..983D,
       author = {{D{\'\i}az}, R.~F. and {Almenara}, J.~M. and {Santerne}, A. and {Moutou}, C. and {Lethuillier}, A. and {Deleuil}, M.},
        title = "{PASTIS: Bayesian extrasolar planet validation - I. General framework, models, and performance}",
      journal = {\mnras},
     keywords = {methods: statistical, techniques: photometric, techniques: radial velocities, planetary systems, Astrophysics - Earth and Planetary Astrophysics},
         year = 2014,
        month = jun,
       volume = {441},
       number = {2},
        pages = {983-1004},
          doi = {10.1093/mnras/stu601},
archivePrefix = {arXiv},
       eprint = {1403.6725},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2014MNRAS.441..983D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

@ARTICLE{2015MNRAS.451.2337S,
       author = {{Santerne}, A. and {D{\'\i}az}, R.~F. and {Almenara}, J. -M. and {Bouchy}, F. and {Deleuil}, M. and {Figueira}, P. and {H{\'e}brard}, G. and {Moutou}, C. and {Rodionov}, S. and {Santos}, N.~C.},
        title = "{PASTIS: Bayesian extrasolar planet validation - II. Constraining exoplanet blend scenarios using spectroscopic diagnoses}",
      journal = {\mnras},
     keywords = {methods: data analysis, techniques: radial velocities, techniques: spectroscopic, binaries: spectroscopic, planetary systems, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
         year = 2015,
        month = aug,
       volume = {451},
       number = {3},
        pages = {2337-2351},
          doi = {10.1093/mnras/stv1080},
archivePrefix = {arXiv},
       eprint = {1505.02663},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.2337S},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
