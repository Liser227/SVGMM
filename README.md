# SVGMM

This is a unsurpervised PolSAR classification code for training SVGMM, and you can use it in your model to test your own PolSAR datasets.

You can download various PolSAR images from EASA or AFS website: 1) https://earth.esa.int/web/polsarpro/airborne-data-sources. https://uavsar.jpl.nasa.gov/cgi-bin/data.pl. 

The data should be preprocessed via PolSARpro software: https://earth.esa.int/web/polsarpro/home, and saved as C3 form with the dimention of [Sz(1)*Sz(2), D^2], where Sz is the size of image and D is the element number of scattering vector.
