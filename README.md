**LEfSe**
==============

LEfSe (Linear discriminant analysis Effect Size) determines the features
(organisms, clades, operational taxonomic units, genes, or functions)
most likely to explain differences between classes by coupling standard
tests for statistical significance with additional tests encoding
biological consistency and effect relevance.

LEfSe is available as a [Galaxy module](http://huttenhower.org/galaxy/),
a Conda formula, a Docker image, and included in bioBakery (VM and
cloud). For additional information, please refer to the [LEfSe
paper](http://www.ncbi.nlm.nih.gov/pubmed/21702898).

## Installation

LEfSe can be installed with Conda or run from a Docker image. Please
note, if you are using bioBakery (Vagrant VM or cloud) you do not need
to install LEfSe because the tool and its dependencies are already
installed.

Install with Conda: `$ conda install -c bioconda lefse`

Install with Docker: `$ docker run -it biobakery/lefse bash`

LEfSe requires R v. 3.6 or higher and the R libraries survival, mvtnorm, modeltools, coin, MASS. 

We provide support for LEfSe users. Please join our [bioBakery Support Forum](https://forum.biobakery.org/c/Downstream-analysis-and-statistics/LEfSe) designated specifically for LEfSe users. 
