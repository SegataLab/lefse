"""
Author: Timothy Tickle
Description: Project constants.
"""

#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#####################################################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = "MIT"
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

##
#Used to test the FileIO class
class ConstantsBreadCrumbs():
    """
    Class to hold project constants.
    """

    #Character Constants
    c_strComma = ','
    c_strColon = ':'
    c_strConfigFileHeaderChar = '['
    c_strConfigFileCommentChar = '#'
    c_strEndline = '\n'
    c_strExtDelim = '.'
    c_cFastaIDLineStart = '>'
    c_strPathDelim = '/'
    c_cPipe = '|'
    c_cQuote = '\"'
    c_cTab = '\t'
    c_strWhiteSpace = ' '
    c_matrixFileDelim = '\t'

    c_strBreadCrumbsSVMSpace = c_strWhiteSpace

    #Default values for missing data in the Abundance Table
    c_strEmptyAbundanceData = "0"
    c_strEmptyDataMetadata = "NA"
    c_strSVMNoSample = "-"

    lNAs = list(set(["NA","na","Na","nA",c_strEmptyDataMetadata]))

    #TODO remove
    #Reference to circlader
    c_strCircladerScript = "circlader/circlader.py"

    #AbundanceTable
    #Suffix given to a file that is check with the checkRawDataFile method
    OUTPUT_SUFFIX = "-checked.pcl"

    #BIOM related
    #PCL File metadata defaults (many of these come from biom file requirements
    #ID
    c_strIDKey = "id"
    c_strDefaultPCLID = None

    #File date
    c_strDateKey = "date"

    #File format type
    c_strFormatKey = "format"
    c_strDefaultPCLFileFormateType = "PCL"

    #File generation source
    c_strSourceKey = "source"
    c_strDefaultPCLGenerationSource = None

    #File type
    c_strTypekey = "type"
    c_strDefaultPCLFileTpe = None

    #Allowable file types for biom files
    c_strOTUType = "OTU"
    c_strOTUBIOMType = "OTU table"
    c_strPathwayType = "Pathway"
    c_strPathwayBIOMType = "Pathway table"
    c_strFunctionType = "Function"
    c_strFunctionBIOMType = "Function table"
    c_strOrthologType = "Ortholog"
    c_strOrthologBIOMType = "Ortholog table"
    c_strGeneType = "Gene"
    c_strGeneBIOMType = "Gene table"
    c_strMetaboliteType = "Metabolite"
    c_strMetaboliteBIOMType = "Metabolite table"
    c_strTaxonType = "Taxon"
    c_strTaxonBIOMType = "Taxon table"
    c_dictFileType = {c_strOTUType:c_strOTUBIOMType, c_strPathwayType:c_strPathwayBIOMType, c_strFunctionType:c_strFunctionBIOMType, c_strOrthologType:c_strOrthologBIOMType, c_strGeneType:c_strGeneBIOMType, c_strMetaboliteType:c_strMetaboliteBIOMType, c_strTaxonType:c_strTaxonType}

    #File URL
    c_strURLKey = "url"
    c_strDefaultPCLURL = None
    c_strFormatUrl =  "format_url"

    #File sparse matrix
    c_strSparsityKey = "sparsity"
    c_fDefaultPCLSparsity = False

    # BIOM related Data
    # Data shape
    c_strDataShapeKey = "shape"

	######################################################################
	# Constants related to biom import and export files                  #
	######################################################################
    # Biom file extension
    c_strBiomFile = "biom"
    c_BiomTaxData = "BiomTaxData"
    c_MetadataID = "column_metadata_id"
    c_Metadata = "Metadata"
    c_metadata_lowercase = "metadata"
    c_sLastMetadata = "sLastMetadata"
    c_columns = "columns"	
    c_rows = "rows"
    c_ascii = "ascii"	
    c_ignore = "ignore"	
    c_Dtype = "Dtype"	
    c_ID = "ID"	
    c_id_lowercase = "id"	
    c_f4 = "f8"		
    c_biom_file_generated_by = "BreadCrumbs"
    c_strPCLFile = "pcl"
    c_taxonomy = "taxonomy"
    c_dRowsMetadata = "dRowsMetadata"
    c_BiomFileInfo = "BiomFileInfo"
    c_MatrixTtype = "matrix_type"
    c_GeneratedBy = "generated_by"
    c_MetadataEntriesTotal = "MetadataEntriesTotal"
    c_MaximumLength = "MaximumLength"


    def __init__(self):
      pass
