# BRITE_decor

Dear all,

This is a bundled version of different Python routines, which may prove to be useful for the detrending and studying data taken by the BRITE nano-satellites.

# Requirements
'Standard' Python packages are expected to be installed. These are:
- numpy
- scipy
- pylab / matplotlib
- lmfit (if you want to do a boundend least-squares minimisation; routines marked with lmfit in the name)

In addition, all the plotting routines - currently - have a line stating  "import CheckMatplotlib" in the preambule.  Off course, you should remove this line.

# Examples
So far, I have included three general routines, showing the general use of some of the routines.  These routines are really basic examples, with minor commenting (ToDo), for guidance usage, since each star and each nano-satellite behave slightly different, it is difficult to construct 'general' user-friendly routines.  Therefore, treat them as examples to construct your own dedicated scripts for your current target.
