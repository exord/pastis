# Installation

### Environment variables
PASTIS requires two environment variables to be set:

- PASTISPATH

- PASTISLIB

### Libraries

### JKTEBOP
PASTIS relies on some subroutines from the JKTEBOP package by Southworth [ADD LINK].
To make them available to PASTIS, an extension module must be built using f2py:
```
f2py -c task2_v28_components.f -m task2_components
```
The resulting file, task2_components.so must be located in the directory $PASTISLIB/fortran, which is included in the sys.path variable at startup.

