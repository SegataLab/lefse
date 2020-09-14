## script changed from original by for bioconda-lefse by jfg 2017.10.31
## adapted from the original lefse run.sh available at:
## https://bitbucket.org/nsegata/lefse/src/54694b4b0d9e335ff1ecafff8db4f1e0cf7004da/example/run.sh?at=default&fileviewer=file-view-default
##
##
## note: ensure script names below match those in the conda install directory, otherwise will error: 'not found':
## scripts usually visible in ~/miniconda3/envs/lefse/bin/ if was installed to default directory
## 


# Download a 3-classes example (with subclasses and subjects) from huttenhower.sph.harvard.edu
# It is a small subset of the HMP 16S dataset for finding biomarkers characterizing
# different level of oxygen availability in different bodysites
wget https://github.com/biobakery/biobakery/raw/master/demos/biobakery_demos/data/lefse/input/hmp_small_aerobiosis.txt -O hmp_aerobiosis_small.txt

# As using LEfSe through bioconda, need to activate the LEfSe installation:
source activate lefse

# Running the LEfSe commands with -h gives the list of available options

# lefse-format_input.py convert the input data matrix to the format for LEfSe.
#
# In this example we use the class information in the first line (-c 1)
# the subclass in the second line (-s 2) and the subject in the third (-u 3).
# If the subclass or the subject are not present in the data you need to set
# the value -1 for them.
# -o 1000000 scales the feature such that the sum (of the same taxonomic leve)
# is 1M: this is done only for obtaining more meaningful values for the LDA score
lefse-format_input.py hmp_aerobiosis_small.txt hmp_aerobiosis_small.in -c 1 -s 2 -u 3 -o 1000000

# run_lefse.py performs the actual statistica analysis
#
# Apply LEfSe on the formatted data producing the results (to be further processed
# for visualization with the other modules). The option available
# can be listed using the -h option 
run_lefse.py hmp_aerobiosis_small.in hmp_aerobiosis_small.res

# lefse-plot_res.py visualizes the output
#
# Plot the list of biomarkers with their effect size
# Severak graphical options are available for personalizing the output
lefse-plot_res.py hmp_aerobiosis_small.res hmp_aerobiosis_small.png

# plot_cladogram.py visualizes the output on a hierarchical tree
#
# Plot the representation of the biomarkers on the hierarchical tree
# specified in the input data (using | in the name of the features)
# In this case we will obtain the RDP taxonomy.
# This is an early implementation of the module. I'm working on an improved version
# that will be released independently from LEfSe
lefse-plot_cladogram.py hmp_aerobiosis_small.res hmp_aerobiosis_small.cladogram.png --format png

# Create a directory for storing the raw-data representation of the discovered biomarkers
mkdir biomarkers_raw_images

# lefse-plot_features.py visualizes the raw-data features
#
# The module for exporting the raw-data representation of the features.
# With the default options we will obtain the images for all the features that are
# detected as biomarkers
lefse-plot_features.py hmp_aerobiosis_small.in hmp_aerobiosis_small.res biomarkers_raw_images/

## Turn lefse back off
source deactivate lefse

## bonus: seasonal greetings
# echo '~ Oíche Shamhna féile dhuit!'
