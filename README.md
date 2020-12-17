# MPhys-RRL

Oh hello! Are you here for the tour? Come on in! Please excuse me while I slowly lose my miiiiind

## final_scripts

Contains (mostly) final versions of the code being used for the project.

### ***BIGSCRIPT***:
Runs only on individual, known RRL, from initial photometry through to plotting PL relations and calculating distances. Currently works for both PAL5 and PAL13 - cluster, channel etc. can be selected in the **Initial conditions** section. Some notes:

- ePSF models need to be generated from either the PAL5 or PAL13 version of ***PHOTOMETRY_SCHMOTOMETRY*** for now; these are saved to file using *pickle* (as are some of the photometry tables generated) and imported back in ***BIGSCRIPT***.
- There are two photometry modes implemented, **aperture** and **PSF**. **aperture** photometry alone was mostly used during testing and writing the scripts, and may not currently work because I don't have the energy to go and tweak it right now. **PSF** photometry is the final approach taken with all the data, and is what's being used to generate final results.
- Still needs to be commented/annotated properly, sorryyyy

### ***PHOTOMETRY_SCHMOTOMETRY***:
Similar to ***BIGSCRIPT***, but this runs on the whole frame, non-variables and variables alike, so it's good for looking at measures of variability and other fun stuffs. There are separate versions for PAL5 and PAL13 currently. These are still a little bloated, and also need commenting properly.

### ***Data***:
- RRL databases gathered mostly from Christine Clement's Catalogue of Variable Stars (and some Gaia data) are in *final_scripts/data/reference/*, and are required when running any of the main scripts.
- Also, though they're not uploaded for storage's sake, *.fits* data files for each cluster should be placed in the *final_scripts/data/[cluster_name]/* folder in the same format and folder structure as provided by Vicky.


## test_photometry

Basically the working directory for most of the project. Includes big photometry scripts, tests for various features, and a few wee ideas I occasionally refer back to, but for the most part can be ignored.
