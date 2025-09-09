# IPBA_SU
A new slope unit delineation method based on geomorphon landform classification and invasion percolation-based algorithm 

Readme

This respository contains code and test data for the paper "A new slope unit delineation method based on geomorphon landform classification and invasion percolation-based algorithm" by luo et al. submitted to Computers & Geosciences. WL 9/8/2025.

It contains the following files:

preprocessing.py: this code preprocesses input DEM to produce ridge/valley lines as raster as input for IBPA algorithm; see comments inside the code.

landslides.zip: this contains the full IPBA code after unzip; input data is datasets subfolder; output will be in results subfolder.
example run: "python3 main.py bt bt_valley_rg bt_ridge_rg 1.0 1.0 1.0 32 all1 T"
                where bt = the input DEM in tif, 
				bt_valley_rg = input valley line raster from preprocessing, 
				bt_ridge_rg = input ridge line raster,
				1.0, 1.0, 1.0 = the three parameters p1, p2, p3. p3 should always be set to 1.
				all1 = an identifier in the output file name to indicate the input parameters;
		       	in this case all 3 parameters are 1. you may use 070710 if p1=0.7,p2=0.7,p3=1

postprocessing.py: this code postprocesses the initial slope units from IBPA; see comments inside the code.

bt.tif: the test dem file in .tif format.
