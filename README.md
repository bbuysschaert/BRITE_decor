# BRITE_decor
We have included a collection of Python routines, which we use to detrend, reduce, analyse, and study the BRITE photometry provided by the BRITE constellation (see http://www.univie.ac.at/brite-constellation/).

At present, every routine has been placed into the 'devel' repository, meaning it is under development.  Feel free to test these routines on your own data, so we understand which routines are mature enough to push them to a 'stable' repository.

# Requirements
Most of the routines were tested under Python 2.7 (Enthought) on a Linux machine.  However, Python in itself is sufficiently versatile, so it should run on Python 3.x.  However, we are dependent on a some common (scientific) Python packages.  We mention these below, together with their respective version, used for testing.  These versions are currently up-to-date.
- numpy (1.9.2)
- scipy (0.16.1)
- pylab / matplotlib (1.5.1)
- lmfit (0.8.0): if you want to do a boundend least-squares minimisation.  These routines are currently not used in the examples, and have very little documentation.  We note that lmfit is no longer under active development.
- statsmodels (0.6.1): easily installed with Canopy / Anaconda Python. Otherwise, you can follow one of both links for information on how to install it: http://statsmodels.sourceforge.net/install.html and http://statsmodels.sourceforge.net/devel/install.html .
 
# Installation
At present, we recommend you to download the routines from the 'devel' repository, using svn. After setting up the proper Python path, everything should run smoothly (unless we made a mistake). Below, we show you how to download the package through svn.

> cd /your/favorite/path/location (-- we recommend /python)

> svn checkout https://github.com/bbuysschaert/BRITE_decor.git/trunk/devel/BRITE_decor

Next, you set up your Python path (in your .bashrc or .bash_profile).

$ export PYTHONPATH=/your/favorite/path/location/:$PYTHONPATH

NOTES:
- You do not include the "~/BRITE_decor" in the Python path.
- Do not forget to resource your .bashrc file ('close the current terminal and / or source .bashprofile').
- Try not to rename the "BRITE_decor", since the routines are explicitly looking for it.

UPDATING:
You can keep your package up-to-date by using 'svn update' while being in your BRITE_decor repository.

# Examples
At present, we have included very minimal example_routines (with some data).  These show some examples of the outlier clipping (step-by-step, some interactive), the PSFdetrending (showing what it can do, most difficult part will be to set up your proper temperaturebins), 'standard' detrending.

You can download this example directory (after installing BRITE_decor) as follows:

svn checkout https://github.com/bbuysschaert/BRITE_decor.git/trunk/examples/

# FAQs

Q: I want to use these routines in my own coding.  How do I properly import them?
A: Currently you can use 

> import BRITE_decor.directory.package as xxxxx (example: import BRITE_decor.inout.load as loadBRITE) 

We are trying the make it as such so you can do an 'import BRITE_decor.directory', and have all the respective routines from the packages.  However, we are still editing the proper directory structure.
