# m-dwarf-flares

This is a model that can generate realistic m - dwarf flare instances and save them to an LCLIB file. Please refer to the instructions below to run the model.

## Setting Up the environment :

Clone the repository and run :

`pip install -r requirements.txt`

We recommend using the a virtual environment to install the requirements and run the project.

## Tables and Data :

**Note:** Please note that all tables are not required for the functioning of the model - some were just used for testing.  

To download and save the light curves and dust maps, run 

`python download_initial_data.py`

Alternatively, `generator.py ` will download the light curves and store them as required. However the dust maps will need to be fetched manually. You can use the code from `download_initial_data.py` to do this.

Table 1 and Table 2 are already included in the repository. The Kepler Input Catalogue will be downloaded by this script too, or can be downloaded manually and placed in the data_files directory for all the functions to work correctly.

____

### Table 1 : AN ALL-SKY CATALOG OF BRIGHT M DWARFS

Article Link:
https://iopscience.iop.org/article/10.1088/0004-6256/142/4/138

Table Link:
https://cfn-live-content-bucket-iop-org.s3.amazonaws.com/journals/1538-3881/142/4/138/1/aj403664t1_mrt.txt?AWSAccessKeyId=AKIAYDKQL6LTV7YY2HIK&Expires=1624826709&Signature=FXqYlHsA%2BakeGNFB8bEzgQCuKMc%3D

___

### Table 2 : The Flaring Activity of M Dwarfs in the Kepler Field

Article Link:
https://iopscience.iop.org/article/10.3847/1538-4357/aa8ea2

Table Link:
https://cfn-live-content-bucket-iop-org.s3.amazonaws.com/journals/0004-637X/849/1/36/1/apjaa8ea2t3_mrt.txt?AWSAccessKeyId=AKIAYDKQL6LTV7YY2HIK&Expires=1624827103&Signature=xknVZ6vjGavgmLNK1FMt4xKQpCw%3D

____

### Table 3 : Kepler Input Catalogue

Download the kepler_kic_v10.csv.gz file from

https://archive.stsci.edu/pub/kepler/catalogs/kepler_kic_v10.csv.gz

And save it to the data_files directory.

___

### Filter Data :

The LSST filter data is already included in the repo. They can be found here:

http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=LSST&asttype=
## Installation :

Inside the repo, run

`python -m pip install -e .`
## Downloading relevant data :

To download all the light curve and dustmap data, run:

`mdwarf-download-data`

## Running :

Run 

`python -m m_dwarf_flare -h`

To see the arguments for the LCLIB file generator. 

___

If you want to create 10 simulations and save them to a file, you can run:

`python -m m_dwarf_flare --flare_count 10 --spectrum_class bb_balmer_jump --dir_name sample --file_name sample.TEXT`

This will save the modeled data to a file named `sample.TEXT` in a directory named sample.
___

If you want to tweak the thresholds for the flare filtering, you can do so in `generator.py`.

## Bugs

We are always working to improve the model! If you run into any issues, please feel free to report them in the `Issues` tab.
