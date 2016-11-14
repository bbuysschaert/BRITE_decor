These tutorials were used during the data reduction session given at the Second BRITE conference at Innsbruck, by H. Pablo.

The documents treat the UBr data (Orion I field) for HD 37043 (iota Ori), highlighting some of the provided routines in the BRITE_decor package.  As such, they should not be interpreted as the best or optimal flow of data preparation and correction, but only as mere examples.

Since full python scripts (.py) cannot be executed on a line-by-line basis, the routines were written as an ipython notebook (.ipynb).  Such files can easily be executed with programs such as Jupyter, which opens a new window where each file can be accessed, edited, and executed.

Three distinct ipython notebooks are provided, one for each major step in the data processing.

-initial_setup.ipynb: The timing correction is applied and outlier rejection is performed
-psf_decorrelation.ipynb: The flux is corrected for the variable PSF of the BRITE satellites
-detrending.ipynb: Classical decorrelation between flux and meta-data.

