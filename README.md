# Installation

### Environment variables
PASTIS requires two environment variables to be set:

- PASTISPATH

- PASTISLIB

### External models
External models are used for the stellar atmospheres, evolution tracks, limb darkening parameters, and interstellar extinction. Besides, a library of photmetric filters are necessary to run PASTIS. All these files are provided as part of the repository in the external_models folder.

### JKTEBOP
PASTIS relies on some subroutines from the JKTEBOP package by Southworth [ADD LINK].
To make them available to PASTIS, an extension module must be built using f2py:
```
f2py -c task2_v28_components.f -m task2_components
```

If you are running python3, this should be
```
f2py3 -c task2_v28_components.f -m task2_components3
```
  
The resulting file, task2_components.so must be located in the directory $PASTISLIB/fortran, which is included in the sys.path variable at startup.
The original fortran file is found in the external_models/fortran directory.
