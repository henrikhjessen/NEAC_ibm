# Evolutionary responses of Atlantic cod (_Gadus morhua_) to concurrent fisheries and climate stressors

This is a repository created to host the source code for the Individual-Based-Model built as part of my PhD Thesis work.
The primary model is the file titled **main.f90**, other FORTRAN files are libraries used in the main body.

Included are also R-scripts used to generate plots.

Recommended compiler is **gfortran**, which was used throughout this work as part of the **gcc**, which was employed on both WSL (Ubuntu) and the MacOS Terminal. In addition to the files provided, anyone intending to compile and run this code will also need to install the NetCDF libraries created for FORTRAN: https://docs.unidata.ucar.edu/netcdf-fortran/current/

No proprietary libraries were used in the model, in an attempt to keep the code open-souce and freely available.
