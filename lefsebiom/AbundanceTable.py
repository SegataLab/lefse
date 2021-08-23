"""
Author: Timothy Tickle
Description: Class to abstract an abundance table and methods to run on such a table.
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

import csv
import sys
from .CClade import CClade
from .ConstantsBreadCrumbs import ConstantsBreadCrumbs
import copy
from datetime import date
import numpy as np
import os
import re
import scipy.stats
import string
from .ValidateData import ValidateData


#*************************************************************
#*  import biom                                              *
#* If not found - abnormally exit                            *
#*************************************************************
try:
    from biom.parse import *
    from biom.table import *
except ImportError:
    sys.stderr.write("************************************************************************************************************ \n")
    sys.stderr.write("* Error:   biom environment required to process biom files - Not found - run abnormally terminated         * \n")
    sys.stderr.write("* See http://http://biom-format.org/                                                                       * \n")
    sys.stderr.write("************************************************************************************************************ \n")
    exit(1)

from biom.parse import *
from biom.table import *

c_dTarget    = 1.0
c_fRound    = False
c_iSumAllCladeLevels = -1
c_fOutputLeavesOnly = False

class RowMetadata:
    """
    Holds the row (feature) metadata and associated functions.
    """

    def __init__(self, dictRowMetadata, iLongestMetadataEntry=None, lsRowMetadataIDs=None):
        """ Constructor requires a dictionary or row metadata.
        :param dictRowMetadata:    The row metadata values with the ids as the keys, must be stable (keep order)
        :type:            {string feature id: {'metadata': {'taxonomy': [list of metadata values]}}}
        """

        self.dictRowMetadata = dictRowMetadata
        self.iLongestMetadataEntry = iLongestMetadataEntry
        self.lsRowMetadataIDs = lsRowMetadataIDs

        self.dictMetadataIDs = {}
        # Get the ids for the metadata
        if self.dictRowMetadata:
            for dictMetadata in self.dictRowMetadata.values():
                dictMetadata = dictMetadata.get(ConstantsBreadCrumbs.c_metadata_lowercase, None)

                if dictMetadata:
                    for key,value in dictMetadata.items():
                        if self.dictMetadataIDs.get(key, None):
                            self.dictMetadataIDs[key] = max(self.dictMetadataIDs[key],len(dictMetadata[key]))
                        else:
                            self.dictMetadataIDs[key] = len(dictMetadata[key])

    def funcMakeIDs(self):
        """ There should be a one to one mapping of row metadata ids and the values associated here with a feature ID.
            If not make ids from the key by appending numbers.
        """

        # If there exists a ids list already return (this allows ID order to be given and preserved)
        # Else make a list of IDs
        if self.lsRowMetadataIDs:
            return self.lsRowMetadataIDs

        lsIDs = []
        lsKeys = []

        for key, value in self.dictMetadataIDs.items():
            lsKeys.append( key )
            if value > 1:
                lsIDs.extend( [ "_".join( [ key, str( iIndex ) ] ) for iIndex in range( value ) ] )
            else:
                lsIDs.append( key )
        return [ lsIDs, lsKeys ]

    def funGetFeatureMetadata(self, sFeature, sMetadata):
        """
        Returns a list of values in the order of row metadta ids for a microbial feature given an id.

        :param sFeature:    Feature id to get metadata
        :type:            string
        :param sMetadata:    Metadata to get
        :type:            string
        :return:        list of metadata associated with the metadata
        """
        lsMetadata = []
        if self.dictRowMetadata:
            dictFeature = self.dictRowMetadata.get( sFeature, None )
            if dictFeature:
                dictFeatureMetadata = dictFeature.get(ConstantsBreadCrumbs.c_metadata_lowercase, None)
                if dictFeatureMetadata:
                    lsMetadata = dictFeatureMetadata.get(sMetadata, None)
        return lsMetadata


class AbundanceTable:
    """
    Represents an abundance table and contains common function to perform on the object.

    This class is made from an abundance data file. What is expected is a text file delimited by
    a character (which is given to the object). The first column is expected to be the id column
    for each of the rows. Metadata is expected before measurement data. Columns are samples and
    rows are features (bugs).

    This object is currently not hashable.
    """

    def __init__(self, npaAbundance, dictMetadata, strName, strLastMetadata, rwmtRowMetadata = None, dictFileMetadata = None, lOccurenceFilter = None, cFileDelimiter = ConstantsBreadCrumbs.c_cTab, cFeatureNameDelimiter = "|"):
        """
        Constructor for an abundance table.

        :param    npaAbundance:    Structured Array of abundance data (Row=Features, Columns=Samples)
        :type:    Numpy Structured Array abundance data (Row=Features, Columns=Samples)
        :param    dictMetadata:    Dictionary of metadata {"String ID":["strValue","strValue","strValue","strValue","strValue"]}
        :type:    Dictionary    Dictionary
        :param    npaRowMetdata    Structured Array of row (feature) metadata (optional)
        :type:    Numpy Structured Array abundance data (Row=Features, Columns=Feature metadata)
         :param    strName:    The name of the metadata that serves as the ID for the columns (For example a sample ID)
        :type:    string
        :param    strLastMetadata: The string last metadata name
              :type:    string
        :param    lOccurenceFilter: List of integers used in an occurence filter. [Min abundance, Min sample]
        :type:    List of integers
        :param    cFileDelimiter:    Character used as the delimiter of the file that is read in to create the abundance table.
                                Will also be used to write the abudance table file to a file to keep file consistency.
        :type:    Character delimiter for reading the data in (default = TAB)
        :param    cFeatureNameDelimiter:    Character used as the delimiter of the feature names (column 1). This is useful if the name are complex, for instance consensus lineages in metagenomics.
        :type:    Character delimiter for feature names (default = |)
        """

        ### File Metadata

        #Date
        self.dateCreationDate = dictFileMetadata.get(ConstantsBreadCrumbs.c_strDateKey,None) if dictFileMetadata else None

        #Indicates if the table has been filtered and how
        self._strCurrentFilterState = ""

        #The delimiter from the source file
        self._cDelimiter = cFileDelimiter

        #The feature name delimiter
        self._cFeatureDelimiter = cFeatureNameDelimiter

        #File type
        self.strFileFormatType = dictFileMetadata.get(ConstantsBreadCrumbs.c_strFormatKey,None) if dictFileMetadata else None

        #File generation source
        self.strFileGenerationSource = dictFileMetadata.get(ConstantsBreadCrumbs.c_strSourceKey,None) if dictFileMetadata else None

        #File type
        self.strFileType = dictFileMetadata.get(ConstantsBreadCrumbs.c_strTypekey,None) if dictFileMetadata else None

        #File url
        self.strFileURL = dictFileMetadata.get(ConstantsBreadCrumbs.c_strURLKey,None) if dictFileMetadata else None

        #The id of the file
        self.strId = dictFileMetadata.get(ConstantsBreadCrumbs.c_strIDKey,None) if dictFileMetadata else None

        #The lastmetadata name (which should be preserved when writing the file)
        # Can be a None if biom file is read in.
        self._strLastMetadataName = strLastMetadata

        #The original number of features in the table
        self._iOriginalFeatureCount = -1

        #The name of the object relating to the file it was read from or would have been read from if it exists
        #Keeps tract of changes to the file through the name
        #Will be used to write out the object to a file as needed
        self._strOriginalName = strName

        #The original number of samples in the table
        self._iOriginalSampleCount = -1

        #Data sparsity type
        self.fSparseMatrix = dictFileMetadata.get(ConstantsBreadCrumbs.c_strSparsityKey,False) if dictFileMetadata else False

        ### Data metadata
        #The column (sample) metdata
        self._dictTableMetadata = dictMetadata

        #The row (feature) metadata (Row Metadata object)
        self.rwmtRowMetadata = rwmtRowMetadata

        ### Data

        #The abundance data
        self._npaFeatureAbundance = npaAbundance


        ### Logistical

        #Clade prefixes for biological samples
        self._lsCladePrefixes = ["k__","p__","c__","o__","f__","g__","s__"]

        #This is not a hashable object
        self.__hash__ = None


        ### Prep the object

        self._fIsNormalized = self._fIsSummed = None
        #If contents is not a false then set contents to appropriate objects
        # Checking to see if the data is normalized, summed and if we need to run a filter on it.
        if len(self._npaFeatureAbundance) and self._dictTableMetadata:
            self._iOriginalFeatureCount = self._npaFeatureAbundance.shape[0]
            self._iOriginalSampleCount = len(self.funcGetSampleNames())
        
            self._fIsNormalized = ( max( [max( list(a)[1:] or [0] ) for a in self._npaFeatureAbundance] or [0] ) <= 1 )

            lsLeaves = AbundanceTable.funcGetTerminalNodesFromList( [a[0] for a in self._npaFeatureAbundance], self._cFeatureDelimiter )
            self._fIsSummed = ( len( lsLeaves ) != len( self._npaFeatureAbundance ) )

            #Occurence filtering
            #Removes features that do not have a given level iLowestAbundance in a given amount of samples iLowestSampleOccurence
            if ( not self._fIsNormalized ) and lOccurenceFilter:
                iLowestAbundance, iLowestSampleOccurrence = lOccurenceFilter
                self.funcFilterAbundanceBySequenceOccurence( iLowestAbundance, iLowestSampleOccurrence )
#      else:
#        sys.stderr.write( "Abundance or metadata was None, should be atleast an empty object\n" )

    @staticmethod
    def funcMakeFromFile(xInputFile, cDelimiter = ConstantsBreadCrumbs.c_cTab, sMetadataID = None, sLastMetadataRow = None, sLastMetadata = None,
       lOccurenceFilter = None, cFeatureNameDelimiter="|", xOutputFile = None, strFormat = None):
        """
        Creates an abundance table from a table file.

        :param    xInputFile:    Path to input file.
        :type:    String        String path.
        :param    cDelimiter:    Delimiter for parsing the input file.
        :type:    Character    Character
        :param    sMetadataID:    String ID that is a metadata row ID (found on the first column) and used as an ID for samples
        :type:    String        String ID
        :param sLastRowMetadata: The id of the last (most right column) row metadata
        :type: String    String ID
        :param    sLastMetadata:    The ID of the metadata that is the last metadata before measurement or feature rows.
        :type:    String        String ID
        :param    lOccurenceFilter: List of integers used in an occurence filter. [Min abundance, Min sample]
        :type:    List of integers
        :param    cFeatureNameDelimiter:    Used to parse feature (bug) names if they are complex.
                        For example if they are consensus lineages and contain parent clade information.
        :type:    Character    Delimiting letter
        :param    xOutputFile:    File to output the abundance table which was read in.
        :type:    FileStream or String file path
        :return    AbundanceTable:    Will return an AbundanceTable object on no error. Returns False on error.
        """
        
        #Get output file and remove if existing
        outputFile = open( xOutputFile, "w" ) if isinstance(xOutputFile, str) else xOutputFile
        
        #################################################################################
        #    Check if file is a biom file - if so invoke the biom routine               #
        #################################################################################
        strFileName = xInputFile if isinstance(xInputFile, str) else xInputFile.name

                # Determine the file read function by file extension
        if strFileName.endswith(ConstantsBreadCrumbs.c_strBiomFile) or (strFormat == ConstantsBreadCrumbs.c_strBiomFile):
            BiomCommonArea = AbundanceTable._funcBiomToStructuredArray(xInputFile)
            if  BiomCommonArea:
                lContents = [BiomCommonArea[ConstantsBreadCrumbs.c_BiomTaxData],
                    BiomCommonArea[ConstantsBreadCrumbs.c_Metadata],
                    BiomCommonArea[ ConstantsBreadCrumbs.c_dRowsMetadata],
                    BiomCommonArea[ConstantsBreadCrumbs.c_BiomFileInfo]
                    ]

                # Update last metadata and id if given
                if not sLastMetadata: 
                    strLastMetadata = BiomCommonArea[ConstantsBreadCrumbs.c_sLastMetadata]
            else:
                # return false on failure
                lContents = False
        elif( strFileName.endswith(ConstantsBreadCrumbs.c_strPCLFile) or (strFormat == ConstantsBreadCrumbs.c_strPCLFile) ):    
            #Read in from text file to create the abundance and metadata structures
            lContents = AbundanceTable._funcTextToStructuredArray(xInputFile=xInputFile, cDelimiter=cDelimiter,
                sMetadataID = sMetadataID, sLastMetadataRow = sLastMetadataRow, sLastMetadata = sLastMetadata, ostmOutputFile = outputFile)
        else:    
            print("I do not understand the format to read and write the data as, please use the correct file extension or indicate a type.")
            return( false )

        #If contents is not a false then set contents to appropriate objects
        return AbundanceTable(npaAbundance=lContents[0], dictMetadata=lContents[1], strName=str(xInputFile), strLastMetadata=sLastMetadata, rwmtRowMetadata = lContents[2],
        dictFileMetadata = lContents[3], lOccurenceFilter = lOccurenceFilter, cFileDelimiter=cDelimiter, cFeatureNameDelimiter=cFeatureNameDelimiter) if lContents else False

    #Testing Status: Light happy path testing
    @staticmethod
    def funcCheckRawDataFile(strReadDataFileName, iFirstDataIndex = -1, sLastMetadataName = None, lOccurenceFilter = None, strOutputFileName = "", cDelimiter = ConstantsBreadCrumbs.c_cTab):
        """
        Check the input abundance table.
        Currently reduces the features that have no occurence.
        Also inserts a NA for blank metadata and a 0 for blank abundance data.
        Gives the option to filter features through an occurence filter (a feature must have a level of abundance in a minimal number of samples to be included).
        Either iFirstDataIndex or sLastMetadataName must be given

        :param    strReadDataFileName:    File path of file to read and check.
        :type:    String    File path.
        :param    iFirstDataIndex:    First (row) index of data not metadata in the abundance file.
        :type:    Integer    Index starting at 0.
        :param    sLastMetadataName:    The ID of the last metadata in the file. Rows of measurements should follow this metadata.
        :type:    String
        :param    lOccurenceFilter:    The lowest number of occurences in the lowest number of samples needed for a feature to be kept
        :type:    List[2]    List length 2 [lowest abundance (not normalized), lowest number of samples to occur in] (eg. [2.0,2.0])
        :param    strOutputFileName:    File path of out put file.
        :type:    String    File path.
        :param    cDelimiter:    Character delimiter for reading and writing files.
        :type:    Character    Delimiter.
        :return    Output Path:    Output path for written checked file.
        """

        #Validate parameters
        if (iFirstDataIndex == -1) and (sLastMetadataName == None):
            print("AbundanceTable:checkRawDataFile::Error, either iFirstDataIndex or sLastMetadataNamemust be given.")
            return False

        #Get output file and remove if existing
        outputFile = strOutputFileName
        if not strOutputFileName:
            outputFile = os.path.splitext(strReadDataFileName)[0]+ConstantsBreadCrumbs.OUTPUT_SUFFIX

        #Read input file lines
        #Drop blank lines
        readData = ""
        with open(strReadDataFileName,'rU') as f:
            readData = f.read()
        readData = list(filter(None,readData.split(ConstantsBreadCrumbs.c_strEndline)))

        #Read the length of each line and make sure there is no jagged data
        #Also hold row count for the metadata
        iLongestLength = len(readData[0].split(cDelimiter))
        iMetadataRow = -1
        if not sLastMetadataName:
            sLastMetadataName = "None"
        for iIndex, strLine in enumerate(readData):
            sLineElements = strLine.split(cDelimiter)
        if sLineElements[0] == sLastMetadataName:
            iMetadataRow = iIndex
        iLongestLength = max(iLongestLength, len(sLineElements))

        #If not already set, set iFirstDataIndex
        if iFirstDataIndex < 0:
            iFirstDataIndex = iMetadataRow + 1

        #Used to substitute . to -
        reSubPeriod = re.compile('\.')

        #File writer
        with open(outputFile,'w') as f:

            #Write metadata
            #Empty data is changed to a default
            #Jagged ends are filled with a default
            for strDataLine in readData[:iFirstDataIndex]:
                lsLineElements = strDataLine.split(cDelimiter)
                for iindex, sElement in enumerate(lsLineElements):
                    if not sElement.strip():
                        lsLineElements[iindex] = ConstantsBreadCrumbs.c_strEmptyDataMetadata
                if len(lsLineElements) < iLongestLength:
                    lsLineElements = lsLineElements + ([ConstantsBreadCrumbs.c_strEmptyDataMetadata]*(iLongestLength-len(lsLineElements)))
                f.write(cDelimiter.join(lsLineElements)+ConstantsBreadCrumbs.c_strEndline)

            #For each data line in the table
            for line in readData[iFirstDataIndex:]:
                writeToFile = False
                cleanLine = list()
                #Break line into delimited elements
                lineElements = line.split(cDelimiter)

                #Clean feature name
                sCleanFeatureName = reSubPeriod.sub("-",lineElements[0])

                #For each element but the first (taxa name)
                #Element check to see if not == zero
                #If so add to output
                for element in lineElements[1:]:
                    if(element.strip() in string.whitespace):
                        cleanLine.append(ConstantsBreadCrumbs.c_strEmptyAbundanceData)
                    #Set abundance of 0 but do not indicate the line should be saved
                    elif(element == "0"):
                        cleanLine.append(element)
                    #If an abundance is found set the line to be saved.
                    else:
                        cleanLine.append(element)
                        writeToFile = True

                #Occurence filtering
                #Removes features that do not have a given level iLowestAbundance in a given amount of samples iLowestSampleOccurence
                if lOccurenceFilter:
                    iLowestAbundance, iLowestSampleOccurence = lOccurenceFilter
                    if iLowestSampleOccurence > sum([1 if float(sEntry) >= iLowestAbundance else 0 for sEntry in cleanLine]):
                        writeToFile = False

                #Write to file
                if writeToFile:    
                    f.write(sCleanFeatureName+cDelimiter+cDelimiter.join(cleanLine)+ConstantsBreadCrumbs.c_strEndline)
        return outputFile

    def __repr__(self):
        """
        Represent or print object.
        """
        return "AbundanceTable"

    def __str__(self):
      """
      Create a string representation of the Abundance Table.
      """

      return "".join(["Sample count:", str(len(self._npaFeatureAbundance.dtype.names[1:])),
      os.linesep+"Feature count:", str(len(self._npaFeatureAbundance[self._npaFeatureAbundance.dtype.names[0]])),
      os.linesep+"Id Metadata:", self._npaFeatureAbundance.dtype.names[0],
      os.linesep+"Metadata ids:", str(list(self._dictTableMetadata.keys())),
      os.linesep+"Metadata count:", str(len(list(self._dictTableMetadata.keys()))),
      os.linesep+"Originating source:",self._strOriginalName,
      os.linesep+"Original feature count:", str(self._iOriginalFeatureCount),
      os.linesep+"Original sample count:", str(self._iOriginalSampleCount),
      os.linesep+"Is normalized:", str(self._fIsNormalized),
      os.linesep+"Is summed:", str(self._fIsSummed),
      os.linesep+"Current filtering state:", str(self._strCurrentFilterState),
      os.linesep+"Feature delimiter:", self._cFeatureDelimiter,
      os.linesep+"File delimiter:",self._cDelimiter])

    def __eq__(self, objOther):
        """
        Check if an object is equivalent in data to this object
        Check to make sure that objOther is not None
        Check to make sure objOther is the correct class type
        Check to make sure self and other internal data are the same (exclusing file name)
        Check data and make sure the npa arrays are the same
        Check the metdata to make sure the dicts are the same 
        (will need to sort the keys of the dicts before comparing, they do not guarentee any order.
        """
                # Check for none
        if objOther is None:
            return False

                #Check for object type
        if isinstance(objOther,AbundanceTable) != True:
            return False
        
        #Check feature delimiter
        if self.funcGetFeatureDelimiter() != objOther.funcGetFeatureDelimiter():
            return False

        #Check file delimiter
        if self.funcGetFileDelimiter() != objOther.funcGetFileDelimiter():
            return False

            
            
            
        #************************************************** 
        #* Commented out                                  *  
        #**************************************************             
        #Check name  - Commented out by GW on 2013/09/14 because
        #If we import pcl file into biom file and then export to pcl, the file names might be different but the tables might be the same
                #Check name
        #if self.funcGetName() != objOther.funcGetName():
            #return  False
        

            #Check sample metadata
        #Go through the metadata
        result1 = self.funcGetMetadataCopy()
        result2 = objOther.funcGetMetadataCopy()
        if sorted(result1.keys()) != sorted(result2.keys()):
            return False
        for strKey in result1.keys():
            if strKey not in result2:
                return False
            if result1[strKey] != result2[strKey]:
                return False

        #TODO check the row (feature) metadata

        #TODO check the file metadata
        #Check the ID
        if self.funcGetFileDelimiter() != objOther.funcGetFileDelimiter():
            return False

        #Check the date
        if self.dateCreationDate != objOther.dateCreationDate:
            return False

        #Check the format
        if self.strFileFormatType != objOther.strFileFormatType:
            return False

            

        #************************************************** 
        #* Commented out                                  *  
        #**************************************************             
        #Check source  - Commented out by GW on 2013/09/14 because
        #If we import pcl file into biom file and then export to pcl, the file names might be different but the tables might be the same            
        #Check the source
        #if self.strFileGenerationSource != objOther.strFileGenerationSource:
            #return False

        #Check the type
        if self.strFileType != objOther.strFileType:
            return False

        #Check the URL
        if self.strFileURL != objOther.strFileURL:
            return False

        #Check data
        #TODO go through the data
        #TODO also check the data type
        result1 = self.funcGetAbundanceCopy()
        result2 = objOther.funcGetAbundanceCopy()    
        if len(result1) != len(result2):
            return False

        sorted_result1 = sorted(result1, key=lambda tup: tup[0])
        sorted_result2 = sorted(result2, key=lambda tup: tup[0])        
            
        if sorted_result1 != sorted_result2 :
            return  False
                

        #************************************************** 
        #* Commented out                                  *  
        #**************************************************             
        #Check AbundanceTable.__str__(self)  - Commented out by GW on 2013/09/14 because
        #If we import pcl file into biom file and then export to pcl, the file names might be different but the tables might be the same                        
                
        #Check string representation
        #if AbundanceTable.__str__(self) !=  AbundanceTable.__str__(objOther):
                #return False
                    
        #Check if sample ids are the same and in the same order
        if self.funcGetSampleNames() != objOther.funcGetSampleNames():
            return  False

        return  True
    

    def __ne__(self, objOther):
        return not self == objOther

      
    #Testing Status: Light happy path testing
        #TODO: Tim change static to class methods
    @staticmethod
    def _funcTextToStructuredArray(xInputFile = None, cDelimiter = ConstantsBreadCrumbs.c_cTab, sMetadataID = None, sLastMetadataRow = None, sLastMetadata = None, ostmOutputFile = None):

        """
        Private method
        Used to read in a file that is samples (column) and taxa (rows) into a structured array.

        :param    xInputFile:    File stream or path to input file.
        :type:    String        File stream or string path.
        :param    cDelimiter:    Delimiter for parsing the input file.
        :type:    Character    Character.
        :param    sMetadataID:    String ID that is a metadata row ID (found on the first column) and used as an ID for samples. 
                    If not given it is assumed to be position 0
        :type: String        String ID
        :param sLastMetadataRow: String ID that is the last row metadat id (id of the most right column with row/feature metadata)
        :type:    String        String ID
        :param    sLastMetadata:    The ID of the metadata that is the last metadata before measurement or feature rows.
        :type:    String        String ID
        :param    ostmOutputFile:    Output File to write to if needed. None does not write the file.
        :type:    FileStream or String
        :return    [taxData,metadata,rowmetadata]: Numpy Structured Array of abundance data and dictionary of metadata.
                                         Metadata is a dictionary as such {"ID", [value,value,values...]}
                                         Values are in the order thety are read in (and the order of the sample names).
                                         ID is the first column in each metadata row.
-                                        rowmetadata is a optional Numpy strucured array (can be None if not made)
+                                        rowmetadata is a optional RowMetadata object (can be None if not made)
                                         The rowmetadata and taxData row Ids should match
-                                        [Numpy structured Array, Dictionary, Numpy structured array]
+                                        The last dict is a collection of BIOM fielparameters when converting from a BIOM file
+                                        [Numpy structured Array, Dictionary, Numpy structured array, dict]
        """

        # Open file from a stream or file path
        istmInput = open( xInputFile, 'rU' ) if isinstance(xInputFile, str) else xInputFile
        # Flag that when incremented will switch from metadata parsing to data parsing
        iFirstDataRow = -1
        # Sample id row
        namesRow = None
        # Row metadata names
        lsRowMetadataIDs = None
        # Index of the last row metadata
        iIndexLastMetadataRow = None
        # Holds metadata {ID:[list of values]}
        metadata = dict()
        # Holds the data measurements [(tuple fo values)]
        dataMatrix = []
        # Holds row metadata { sID : [ list of values ] }
        dictRowMetadata = {}
        # Positional index
        iIndex = -1
        # File handle
        csvw = None

        # Read in files
        if ostmOutputFile:
            csvw = csv.writer( open(ostmOutputFile,'w') if isinstance(ostmOutputFile, str) else ostmOutputFile, csv.excel_tab, delimiter = cDelimiter )
        # For each line in the file, and assume the tax id is the first element and the data follows
        for lsLineElements in csv.reader( istmInput, dialect = csv.excel_tab, delimiter = cDelimiter ):
            iIndex += 1
            taxId, sampleReads = lsLineElements[0], lsLineElements[1:]

            # Read through data measurements
            # Process them as a list of tuples (needed for structured array)
            if iFirstDataRow > 0:
                try:
                    # Parse the sample reads, removing row metadata and storing row metadata if it exists
                    if lsRowMetadataIDs:
                        # Build expected dict for row metadata dictionary {string feature id: {'metadata': {metadatakey: [list of metadata values]}}}
                        dictFeature = dict([ [sID, [sKey]] for sID, sKey in zip( lsRowMetadataIDs, sampleReads[ 0 : iIndexLastMetadataRow ]) ])
                        if len( dictFeature ):
                            dictRowMetadata[ taxId ] = { ConstantsBreadCrumbs.c_metadata_lowercase: dictFeature }
                        dataMatrix.append(tuple([taxId] + [( float(s) if s.strip( ) else 0 ) for s in sampleReads[ iIndexLastMetadataRow: ]]))
                    else:
                        dataMatrix.append(tuple([taxId] + [( float(s) if s.strip( ) else 0 ) for s in sampleReads]))
                except ValueError:
                    sys.stderr.write( "AbundanceTable:textToStructuredArray::Error, non-numerical value on data row. File:" + str(xInputFile) +
                        " Row:" + str(lsLineElements) + "\n" )
                    return False
            # Go through study measurements
            else:
                # Read in metadata values, if the entry is blank then give it the default empty metadata value.
                for i, s in enumerate( sampleReads ):
                    if not s.strip( ):
                        sampleReads[i] = ConstantsBreadCrumbs.c_strEmptyDataMetadata

                # If no id metadata (sample ids) is given then the first row is assumed to be the id row, otherwise look for the id for the metadata.
                # Add the metadata to the containing dict
                if ( ( not sMetadataID ) and ( iIndex == 0 ) ) or ( taxId == sMetadataID ):
                    namesRow = lsLineElements
                    # Remove the row metadata ids, these names are for the column ID and the samples ids
                    if sLastMetadataRow:
                        iIndexLastMetadataRow = lsLineElements.index(sLastMetadataRow)
                        lsRowMetadataIDs = namesRow[ 1 : iIndexLastMetadataRow + 1 ]
                        namesRow = [ namesRow[ 0 ] ] + namesRow[ iIndexLastMetadataRow + 1: ]

                        # If the sample metadata dictionary already has entries then remove the row metadata info from it.
                        if len( metadata ) and len( lsRowMetadataIDs ):
                            for sKey, lsValues in metadata.items():
                                metadata[ sKey ] = lsValues[ iIndexLastMetadataRow: ]

                # Set the metadata without row metadata entries
                metadata[taxId] = sampleReads[ iIndexLastMetadataRow: ] if (lsRowMetadataIDs and len( lsRowMetadataIDs )) else sampleReads

                # If the last metadata was just processed switch to data processing
                # If the last metadata name is not given it is assumed that there is only one metadata
                if ( not sLastMetadata ) or ( taxId == sLastMetadata ):
                    iFirstDataRow = iIndex + 1

            # If writing out the data write back out the line read in.
            # This happens at the end so that the above cleaning is captured and written.
            if csvw:
                csvw.writerow( [taxId] + sampleReads )

        if sLastMetadata and ( not dataMatrix ):
            sys.stderr.write( "AbundanceTable:textToStructuredArray::Error, did not find the row for the last metadata ID. File:" + str(xInputFile) +
                " Identifier:" + sLastMetadata + "\n" )
            return False

        # Make sure the names are found
        if namesRow == None:
            sys.stderr.write( "AbundanceTable:textToStructuredArray::Error, did not find the row for the unique sample/column. File:" + str(xInputFile) +
                " Identifier:" + sMetadataID + "\n" )
            return False

        # Now we know the longest taxId we can define the first column holding the tax id
        # Gross requirement of Numpy structured arrays, a = ASCII followed by max # of characters (as a string)
        longestTaxId = max( len(a[0]) for a in dataMatrix )
        dataTypeVector = [(namesRow[0],'a' + str(longestTaxId*2))] + [(s, "f4") for s in namesRow[1:]]
        # Create structured array
        taxData = np.array(dataMatrix,dtype=np.dtype(dataTypeVector))

        # Returns a none currently because the PCL file specification this originally worked on did not have feature metadata
         # Can be updated in the future.
        # [Data (structured array), column metadata (dict), row metadata (structured array), file metadata (dict)]
        return [taxData, metadata, RowMetadata(dictRowMetadata = dictRowMetadata, lsRowMetadataIDs = lsRowMetadataIDs), {
                    ConstantsBreadCrumbs.c_strIDKey:ConstantsBreadCrumbs.c_strDefaultPCLID,
                    ConstantsBreadCrumbs.c_strDateKey:str(date.today()),
                    ConstantsBreadCrumbs.c_strFormatKey:ConstantsBreadCrumbs.c_strDefaultPCLFileFormateType,
                    ConstantsBreadCrumbs.c_strSourceKey:ConstantsBreadCrumbs.c_strDefaultPCLGenerationSource,
                    ConstantsBreadCrumbs.c_strTypekey:ConstantsBreadCrumbs.c_strDefaultPCLFileTpe,
                    ConstantsBreadCrumbs.c_strURLKey:ConstantsBreadCrumbs.c_strDefaultPCLURL,
                    ConstantsBreadCrumbs.c_strSparsityKey:ConstantsBreadCrumbs. c_fDefaultPCLSparsity}]

#    def funcAdd(self,abndTwo,strFileName=None):
#        """
#        Allows one to add an abundance table to an abundance table. They both must be the same state of normalization or summation
#        or they will be summed or normalized if one of the two are.
#
#        :param    abndTwo:    AbundanceTable object 2
#        :type:    AbundanceTable
#        :return    AbudanceTable:
#        """
#
#        #Check summation and normalization
#        if(self.funcIsSummed() or abndTwo.funcIsSummed()):
#            self.funcSum()
#            abndTwo.funcSum()
#        if(self.funcIsNormalized() or abndTwo.funcIsNormalized()):
#            self.funcNormalize()
#            abndTwo.funcNormalize()
#
#        #Normalize Feature names
#            #Get if the abundance tables have clades
#            fAbndInputHasClades = self.funcHasFeatureHierarchy()
#            fAbndCompareHasClades = abndTwo.funcHasFeatureHierarchy()
#
#            if(fAbndInputHasClades or fAbndCompareHasClades):
#            #If feature delimiters do not match, switch
#            if not self.funcGetFeatureDelimiter() == abndTwo.funcGetFeatureDelimiter():
#                abndTwo.funcSetFeatureDelimiter(self.funcGetFeatureDelimiter())
#
#            #Add prefixes if needed.
#                    self.funcAddCladePrefixToFeatures()
#                abndTwo.funcAddCladePrefixToFeatures()
#
#        #Get feature Names
#        lsFeatures1 = self.funcGetFeatureNames()
#        lsFeatures2 = abndTwo.funcGetFeatureNames()
#
#        #Make one feature name list
#        lsFeaturesCombined = list(set(lsFeatures1+lsFeature2))
#
#        #Add samples by features (Use 0.0 for empty data features, use NA for empty metadata features)
#        
#
#        #Combine metadata
#        dictMetadata1 = self.funcGetMetadataCopy()
#        dictMetadata2 = abndTwo.funcGetMetadataCopy()
#
#        #Get first table metadata and add NA for metadata it is missing for the length of the current metadata
#        lsMetadataOnlyInTwo = list(set(dictMetadata2.keys())-set(dictMetadata1.keys()))
#        dictCombinedMetadata = dictMetadata1
#        lsEmptyMetadata = ["NA" for i in range(self.funcGetSampleCount())]
#        for sKey in lsMetadataOnlyInTwo:
#            dictCombinedMetadata[sKey]=lsEmptyMetadata
#        #Add in the other metadata dictionary
#        lsCombinedKeys = dictCombinedMetadata.keys()
#        lsEmptyMetadata = ["NA" for i in range(abndTwo.funcGetSampleCount())]
#        for sKey in lsCombinedKeys():
#            if sKey in dictMetadata2:
#                dictCombinedMetadata[sKey] = dictCombinedMetadata[sKey]+dictMetadata2[sKey]
#            else:
#                dictCombinedMetadata[sKey] = dictCombinedMetadata[sKey]+lsEmptyMetadata
#
#        #Make Abundance table
#        return AbundanceTable(npaAbundance=npaAbundance,
#                dictMetadata = dictCombinedMetadata,
#                strName = strFileName if strFileName else os.path.splitext(self)[0]+"_combined_"+os.path.splitext(abndTwo)[0],
#                strLastMetadata = self.funcGetLastMetadataName(),
#                cFileDelimiter = self.funcGetFileDelimiter(), cFeatureNameDelimiter=self.funcGetFeatureDelimiter())

    #TODO This does not adjust for sample ordering, needs to
    def funcAddDataFeature(self, lsNames, npdData):
        """
        Adds a data or group of data to the underlying table.
        Names should be in the order of the data
        Each row is considered a feature (not sample).

        :param lsNames:    Names of the features being added to the data of the table
        :type: List    List of string names
        :param npdData: Rows of features to add to the table
        :type:    Numpy array accessed by row.
        """
        if ( self._npaFeatureAbundance == None ):
            return False

        # Check number of input data rows
        iDataRows = npdData.shape[0]
        if (len(lsNames) != iDataRows):
            print("Error:The names and the rows of data features to add must be of equal length")

        # Grow the array by the neccessary amount and add the new rows
        iTableRowCount = self.funcGetFeatureCount()
        iRowElementCount = self.funcGetSampleCount()
        self._npaFeatureAbundance.resize(iTableRowCount+iDataRows)
        for iIndexData in range(iDataRows):
            self._npaFeatureAbundance[iTableRowCount+iIndexData] = tuple([lsNames[iIndexData]]+list(npdData[iIndexData]))

        return True

    #TODO This does not adjust for sample ordering, needs to
    def funcAddMetadataFeature(self,lsNames,llsMetadata):
        """
        Adds metadata feature to the underlying table.
        Names should be in the order of the lists of metadata
        Each internal list is considered a metadata and paired to a name
        """
        if ( self._dictTableMetadata == None ):
            return False

        # Check number of input data rows
        iMetadataCount = len(llsMetadata)
        if (len(lsNames) != iMetadataCount):
            print("Error:The names and the rows of metadata features to add must be of equal length")

        # Add the metadata
        for tpleMetadata in zip(lsNames,llsMetadata):
            self._dictTableMetadata[tpleMetadata[0]]=tpleMetadata[1]
        return True

    #2 test Cases
    def funcSetFeatureDelimiter(self, cDelimiter):
        """
        Changes the feature delimiter to the one provided.
        Updates the feature names.

        :param    cDelimiter:    Character feature delimiter
        :type:    Character
        :return    Boolean:    Indicator of success or not (false)
        """
        if ( self._npaFeatureAbundance == None ):
            return False
        cDelimiterCurrent = self.funcGetFeatureDelimiter()
        if ( not cDelimiter or not cDelimiterCurrent):
            return False

        #Make new feature names
        lsNewFeatureNames = [sFeatureName.replace(cDelimiterCurrent,cDelimiter) for sFeatureName in self.funcGetFeatureNames()]
        
        #Update new feature names to abundance table
        if (not self.funcGetIDMetadataName() == None):
            self._npaFeatureAbundance[self.funcGetIDMetadataName()] = np.array(lsNewFeatureNames)

        #Update delimiter
        self._cFeatureDelimiter = cDelimiter
        return True

    #Happy path tested
    def funcGetSampleNames(self):
        """
        Returns the sample names (IDs) contained in the abundance table.

        :return    Sample Name:    A List of sample names indicated by the metadata associated with the sMetadataId given in table creation.
                                A list of string names or empty list on error as well as no underlying table.
        """

        return self._npaFeatureAbundance.dtype.names[1:] if ( len(self._npaFeatureAbundance) > 1 ) else []

    #Happy Path Tested
    def funcGetIDMetadataName(self):
        """
        Returns the metadata id.

        :return    ID:    The metadata id (the sample Id).
                      Returns none on error.
        """

        return self._npaFeatureAbundance.dtype.names[0] if ( len(self._npaFeatureAbundance) > 1  ) else None

    #Happy path tested
    def funcGetAbundanceCopy(self):
        """
        Returns a deep copy of the abundance table.

        :return    Numpy Structured Array:    The measurement data in the Abundance table. Can use sample names to access each column of measurements.
                                       Returns none on error.
        """

        return self._npaFeatureAbundance.copy() if ( self._npaFeatureAbundance != None ) else None

    #Happy path tested
    def funcGetAverageAbundancePerSample(self, lsTargetedFeatures):
        """
        Averages feature abundance within a sample.
    
        :param    lsTargetedFeatures:    String names of features to average
        :type:    List of string names of features which are measured
        :return    List: List of lists or boolean (False on error). One internal list per sample indicating the sample and the feature's average abudance
            [[sample,average abundance of selected taxa]] or False on error
        """

        #Sample rank averages [[sample,average abundance of selected taxa]]
        sampleAbundanceAverages = []
        
        sampleNames = self.funcGetSampleNames()
        allTaxaNames = self.funcGetFeatureNames()
        #Get an abundance table compressed to features of interest
        abndReducedTable = self.funcGetFeatureAbundanceTable(lsTargetedFeatures)
        if abndReducedTable == None:
            return False  

        #If the taxa to be selected are not in the list, Return nothing and log
        lsMissing = []
        for sFeature in lsTargetedFeatures:
            if not sFeature in allTaxaNames:
                lsMissing.append(sFeature)
            else:
                #Check to make sure the taxa of interest is not average abundance of 0
                if not abndReducedTable.funcGetFeatureSumAcrossSamples(sFeature):
                    lsMissing.append(sFeature)
        if len(lsMissing) > 0:
            sys.stderr.write( "Could not find features for averaging: " + str(lsMissing) )
            return False

        #For each sample name get average abundance
        for sName in sampleNames:
            npaFeaturesSample = abndReducedTable.funcGetSample(sName)
            sampleAbundanceAverages.append([sName,sum(npaFeaturesSample)/float(len(npaFeaturesSample))])

        #Sort based on average
        return sorted(sampleAbundanceAverages, key = lambda sampleData: sampleData[1], reverse = True)

    #Happy path tested 1
    def funcGetAverageSample(self):
        """
        Returns the average sample of the abundance table.
        This average sample is made of the average of each feature.
        :return list: A list of averages in the order of the feature names.
        """

        ldAverageSample = []
        #If there are no samples then return empty list.
        if len(self.funcGetSampleNames()) < 1:
            return ldAverageSample

        #If there are samples return the average of each feature in the order of the feature names.
        for sFeature in self._npaFeatureAbundance:
            npFeaturesAbundance = list(sFeature)[1:]
            ldAverageSample.append(sum(npFeaturesAbundance)/float(len(npFeaturesAbundance)))

        return ldAverageSample

    #Tested 2 cases
    def funcHasFeatureHierarchy(self):
        """
        Returns an indicator of having a hierarchy in the features indicated by the existance of the
        feature delimiter.

        :return    Boolean:    True (Has a hierarchy) or False (Does not have a hierarchy)
        """

        if ( self._npaFeatureAbundance == None ):
            return None
        cDelimiter = self.funcGetFeatureDelimiter()
        if ( not cDelimiter ):
            return False

        #For each feature name, check to see if the delimiter is in the name
        for sFeature in self.funcGetFeatureNames():
            if cDelimiter in sFeature:
                return True
        return False

    def funcGetCladePrefixes(self):
        """
        Returns the list of prefixes to use on biological sample hierarchy

        :return    List:    List of strings
        """
        return self._lsCladePrefixes

    #3 test cases
    def funcAddCladePrefixToFeatures(self):
        """
        As a standardized clade prefix to indicate biological clade given hierarchy.
        Will not add a prefix to already prefixes feature names.
        Will add prefix to feature names that do not have them or clades in a feature name that
        do not have them while leaving ones that do as is.

        :return    Boolean:    True (Has a hierarchy) or False (Does not have a hierarchy)
        """

        if ( self._npaFeatureAbundance == None ):
            return None
        cDelimiter = self.funcGetFeatureDelimiter()
        lsPrefixes = self.funcGetCladePrefixes()
        iPrefixLength = len(lsPrefixes)
        if ( not cDelimiter ):
            return False

        #Append prefixes to feature names
        lsUpdatedFeatureNames = []
        lsFeatureNames = self.funcGetFeatureNames()

        for sFeatureName in lsFeatureNames:
            lsClades = sFeatureName.split(cDelimiter)
            #If there are not enough then error
            if(len(lsClades) > iPrefixLength):
                print("Error:: Too many clades given to be biologically meaningful")
                return False
            lsUpdatedFeatureNames.append(cDelimiter.join([lsPrefixes[iClade]+lsClades[iClade] if not(lsClades[iClade][0:len(lsPrefixes[iClade])]==lsPrefixes[iClade]) else lsClades[iClade] for iClade in range(len(lsClades))]))

        #Update new feature names to abundance table
        if not self.funcGetIDMetadataName() == None:
            self._npaFeatureAbundance[self.funcGetIDMetadataName()] = np.array(lsUpdatedFeatureNames)

        return True

    #Happy Path Tested
    def funcGetFeatureAbundanceTable(self, lsFeatures):
        """
        Returns a copy of the current abundance table with the abundance of just the given features.

        :param    lsFeatures:    String Feature IDs that are kept in the compressed abundance table.
        :type:    List of strings    Feature IDs (found as the first entry of a filter in the input file.
        :return    AbundanceTable:    A compressed version of the abundance table.
                  On an error None is returned.
        """
        
        if ( self._npaFeatureAbundance == None ) or ( lsFeatures == None ):
            return None

        #Get a list of boolean indicators that the row is from the features list
        lfFeatureData = [sRowID in lsFeatures for sRowID in self.funcGetFeatureNames()]
        #compressed version as an Abundance table
        lsNamePieces = os.path.splitext(self._strOriginalName)
        abndFeature = AbundanceTable(npaAbundance=np.compress(lfFeatureData, self._npaFeatureAbundance, axis = 0),
                    dictMetadata = self.funcGetMetadataCopy(),
                    strName = lsNamePieces[0] + "-" + str(len(lsFeatures)) +"-Features"+lsNamePieces[1],
                    strLastMetadata=self.funcGetLastMetadataName(),
                    cFileDelimiter = self.funcGetFileDelimiter(), cFeatureNameDelimiter= self.funcGetFeatureDelimiter())
        #Table is no longer normalized
        abndFeature._fIsNormalized = False
        return abndFeature

    #Happy path tested
    def funcGetFeatureDelimiter(self):
        """
        The delimiter of the feature names (For example to use on concensus lineages).

        :return    Character:    Delimiter for the feature name pieces if it is complex.
        """

        return self._cFeatureDelimiter

    #Happy path tested
    def funcGetFeatureCount(self):
        """
        Returns the current feature count.

        :return    Count:    Returns the int count of features in the abundance table.
                        Returns None on error.
        """

        return self._npaFeatureAbundance.shape[0] if not self._npaFeatureAbundance is None else 0

    #Happy path tested
    def funcGetFeatureSumAcrossSamples(self,sFeatureName):
        """
        Returns float sum of feature values across the samples.

        :param    sFeatureName: The feature ID to get the sum from.
        :type:    String.
        :return    Double:    Sum of one feature across samples.
        """
        return sum(self.funcGetFeature(sFeatureName))

    def funcGetFeature(self,sFeatureName):
        """
        Returns feature values across the samples.

        :param    sFeatureName: The feature ID to get the sum from.
        :type:    String.
        :return    Double:    Feature across samples.
        """

        for sFeature in self._npaFeatureAbundance:
            if sFeature[0] == sFeatureName:
                return list(sFeature)[1:]
        return None

    #Happy path tested
    def funcGetFeatureNames(self):
        """
        Return the feature names as a list.

        :return    Feature Names:    List of feature names (or IDs) as strings.
                                As an error returns empty list.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance[self.funcGetIDMetadataName()]
        return []

    #Happy path tested
    def funcGetFileDelimiter(self):
        """
        The delimiter of the file the data was read from and which is also the delimiter which would be used to write the data to a file.

        :return    Character:    Delimiter for the parsing and writing the file.
        """

        return self._cDelimiter

    def funcGetLastMetadataName(self):
        """
        Get the last metadata name that seperates abundance and metadata measurements.

        :return string:    Metadata name
        """
        return self._strLastMetadataName

    #Happy path tested
    def funcGetSample(self,sSampleName):
        """
        Return a copy of the feature measurements of a sample.

        :param    sSampleName:    Name of sample to return.    
        :type:    String    
        :return    Sample: Measurements    Feature measurements of a sample.
                Empty numpy array returned on error.
        """

        if (not self._npaFeatureAbundance == None):
            return self._npaFeatureAbundance[sSampleName].copy()
        return np.array([])

    #Happy path tested
    def funcGetMetadata(self, strMetadataName):
        """
        Returns a list of metadata that is associated with the given metadata name (id).

        :param    strMetadataName:    String metadata ID to be returned
        :type:    String    ID
        :return    Metadata:    List of metadata
        """
        
        return copy.deepcopy( self._dictTableMetadata.get(strMetadataName) ) \
            if self._dictTableMetadata else None

    #Happy path tested
    def funcGetMetadataCopy(self):
        """
        Returns a deep copy of the metadata.

        :return    Metadata copy:    {"ID":[value,value...]}
        """

        return copy.deepcopy(self._dictTableMetadata)
        
    #Happy path tested
    def funcGetName(self):
        """
        Returns the name of the object which is the file name that generated it.
        If the object was generated from an Abundance Table (for instance through stratification)
        the name is still in the form of a file that could be written to which is informative
        of the changes that have occurred on the data set.
        :return string: Name
        """
        return self._strOriginalName

    #Happy path tested. could do more
    def funcGetTerminalNodes(self):
        """
        Returns the terminal nodes given the current feature names in the abundance table. The 
        features must contain a consensus lineage or all will be returned.
        :return List:    List of strings of the terminal nodes given the abundance table.
        """
        return AbundanceTable.funcGetTerminalNodesFromList(lsNames=self.funcGetFeatureNames(),cNameDelimiter=self.funcGetFeatureDelimiter())

    #Tested 2 test cases
    @staticmethod
    def funcGetTerminalNodesFromList(lsNames,cNameDelimiter):
        """
        Returns the terminal nodes given the current feature names in the abundance table. The 
        features must contain a consensus lineage or all will be returned.

        :param    lsNames:    The list of string names to parse and filter.
        :type:    List of strings
        :param    cNameDelimiter:    The delimiter for the name of the features.
        :type:    Character    Delimiter
        :return list:    A list of terminal elements in the list (given only the list).
        """

        #Build hash
        dictCounts = dict()
        for strTaxaName in lsNames:
            #Split into the elements of the clades
            lsClades = list(filter(None,strTaxaName.split(cNameDelimiter)))
            #Count clade levels
            iCladeLength = len(lsClades)

            #Evaluate first element
            sClade = lsClades[0]
            dictCounts[sClade] = sClade not in dictCounts

            #Evaluate the rest of the elements
            if iCladeLength < 2:
                continue
            for iIndex in range(1,iCladeLength):
                prevClade = sClade
                sClade = cNameDelimiter.join([sClade,lsClades[iIndex]])
                if sClade in dictCounts:
                    dictCounts[sClade] = dictCounts[prevClade] = False
                else:
                    dictCounts[sClade] = True
                    dictCounts[prevClade] = False

        #Return only the elements that were of count 1
        return list(filter( lambda s: dictCounts[s] == True, dictCounts ))

    #Happy path tested
    def funcIsNormalized(self):
        """
        Returns if the data has been normalized.

        :return    Boolean:    Indicates if the data is normalized.
                           True indicates it the data is normalized.
        """

        return self._fIsNormalized

    #Happy path tested
    def funcIsPrimaryIdMetadata(self,sMetadataName):
        """
        Checks the metadata data associated with the sMetadatName and returns if the metadata is unique.
        This is important to some of the functions in the Abundance Table specifically when translating from one metadata to another.
        
        :param    sMetadataName:    ID of metadata to check for uniqueness.
        :type:    String    Metadata ID.
        :return    Boolean:    Returns indicator of uniqueness.
                            True indicates unique.
        """

        lMetadata = self.funcGetMetadata(sMetadataName)
        if not lMetadata:
            return False
        return (len(lMetadata) == len(set(lMetadata)))

    #Happy path tested
    def funcIsSummed(self):
        """
        Return is the data is summed.

        :return    Boolean:    Indicator of being summed. True indicates summed.
        """

        return self._fIsSummed

    #Happy path tested
    def funcFilterAbundanceByPercentile(self, dPercentileCutOff = 95.0, dPercentageAbovePercentile=1.0):
        """
        Filter on features.
        A feature is removed if it's abundance is not found in the top X percentile a certain percentage of the samples.

        :param    dPercentileCutOff:    The percentile used for filtering.
        :type:    double    A double between 0.0 and 100.0
        :param    dPercentageAbovePercentile:    The percentage above the given percentile (dPercentileCutOff) that must exist to keep the feature.
        :type:    double    Between 0.0 and 100.0
        :return    Boolean:    Indicator of filtering occuring without error. True indicates filtering occuring.
        """

        #No need to do anything
        if(dPercentileCutOff==0.0) or (dPercentageAbovePercentile==0.0):
            return True

        #Sample names
        lsSampleNames = self.funcGetSampleNames()

        #Scale percentage out of 100
        dPercentageAbovePercentile = dPercentageAbovePercentile/100.0

        #Sample count
        iSampleCount = len(lsSampleNames)

        #Get a threshold score of the value at the specified percentile for each sample
        #In the order of the sample names
        ldScoreAtPercentile = [scipy.stats.scoreatpercentile(self._npaFeatureAbundance[lsSampleNames[iIndex]],dPercentileCutOff) for iIndex in range(iSampleCount)]

        #Record how many entries for each feature have a value equal to or greater than the dPercentileCutOff
        #If the percentile of entries passing the criteria are above the dPercentageAbovePercentile put index in list to keep
        liKeepIndices = []
        iSampleCount = float(iSampleCount)
        for iRowIndex, npaRow in enumerate(self._npaFeatureAbundance):
            iCountPass = sum([1 if dValue >= ldScoreAtPercentile[iValueIndex] else 0 for iValueIndex, dValue in enumerate(list(npaRow)[1:])])
            if (iCountPass / iSampleCount) >= dPercentageAbovePercentile:
                liKeepIndices.append(iRowIndex)

        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepIndices,:]

        #Update filter state
        self._strCurrentFilterState += ":dPercentileCutOff=" + str(dPercentileCutOff) + ",dPercentageAbovePercentile=" + str(dPercentageAbovePercentile)

        #Table is no longer normalized
        self._fIsNormalized = False

        return True

    def funcFilterAbundanceByMinValue(self, dMinAbundance = 0.0001, iMinSamples = 3):
        """
        Filter abundance by requiring features to have a minimum relative abundance in a minimum number of samples.
        Will evaluate greater than or equal to the dMinAbundance and iMinSamples.

        :param    dMinAbundance:    Minimum relative abundance.
        :type:    Real    Number Less than 1.
        :param    iMinSamples:    Minimum samples to have the relative abundnace or greater in.
        :type:    Integer    Number greater than 1.
        :return    Boolean:    Indicator of the filter running without error. False indicates error.
        """

        #No need to do anything
        if(dMinAbundance==0) or (iMinSamples==0):
            return True

        #This normalization requires the data to be relative abundance
        if not self._fIsNormalized:
            #sys.stderr.write( "Could not filter by sequence occurence because the data is already normalized.\n" )
            return False

        #Holds which indexes are kept
        liKeepFeatures = []
        for iRowIndex, dataRow in enumerate( self._npaFeatureAbundance ):
            #See which rows meet the criteria and keep the index if needed.
            if len( list( filter( lambda d: d >= dMinAbundance, list(dataRow)[1:] ) ) ) >= iMinSamples:
                liKeepFeatures.append(iRowIndex)

        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]
        #Update filter state
        self._strCurrentFilterState += ":dMinAbundance=" + str(dMinAbundance) + ",iMinSamples=" + str(iMinSamples)

        return True

    #Happy path tested
    def funcFilterAbundanceBySequenceOccurence(self, iMinSequence = 2, iMinSamples = 2):
        """
        Filter occurence by requiring features to have a minimum sequence occurence in a minimum number of samples.
        Will evaluate greater than or equal to the iMinSequence and iMinSamples.

        :param    iMinSequence:    Minimum sequence to occur.
        :type:    Integer    Number Greater than 1.
        :param    iMinSamples:    Minimum samples to occur in.
        :type:    Integer    Number greater than 1.
        :return    Boolean:    Indicator of the filter running without error. False indicates error.
        """

        #No need to do anything
        if(iMinSequence==0) or (iMinSamples==0):
            return True

        #This normalization requires the data to be reads
        if self._fIsNormalized:
            #sys.stderr.write( "Could not filter by sequence occurence because the data is already normalized.\n" )
            return False

        #Holds which indexes are kept
        liKeepFeatures = []
        for iRowIndex, dataRow in enumerate( self._npaFeatureAbundance ):
            #See which rows meet the criteria and keep the index if needed.
            if len( list( filter( lambda d: d >= iMinSequence, list(dataRow)[1:] ) ) ) >= iMinSamples:
                liKeepFeatures.append(iRowIndex)

        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]
        #Update filter state
        self._strCurrentFilterState += ":iMinSequence=" + str(iMinSequence) + ",iMinSamples=" + str(iMinSamples)

        return True
   
    #1 Happy path test
    def funcFilterFeatureBySD(self, dMinSDCuttOff = 0.0):
        """
        A feature is removed if it's abundance is not found to have standard deviation more than the given dMinSDCutoff.

        :param    dMinSDCuttOff:    Standard deviation threshold.
        :type:    Double    A double greater than 0.0.
        :return    Boolean:    Indicator of success. False indicates error.
        """

        #No need to do anything
        if(dMinSDCuttOff==0.0):
            return True

        #Holds which indexes are kept
        liKeepFeatures = []

        #Evaluate each sample
        for iRowIndex, dataRow in enumerate(self._npaFeatureAbundance):
            if(np.std(list(dataRow)[1:])>=dMinSDCuttOff):
                liKeepFeatures.append(iRowIndex)
        
        #Compress array
        self._npaFeatureAbundance = self._npaFeatureAbundance[liKeepFeatures,:]

        #Update filter state
        self._strCurrentFilterState += ":dMinSDCuttOff=" + str(dMinSDCuttOff)

        #Table is no longer normalized
        self._fIsNormalized = False

        return True

        #Happy path tested 2 tests
    def funcGetWithoutOTUs(self):
        """
        Remove features that are terminal otus. Terminal otus are identified as being an integer.
        """

        #Get the feature names
        lsFeatures = self.funcGetFeatureNames()

        #Reduce, filter the feature names
        lsFeatures = [sFeature for sFeature in lsFeatures if not (ValidateData.funcIsValidStringInt(sFeature.split(self.funcGetFeatureDelimiter())[-1]))]

        return self.funcGetFeatureAbundanceTable(lsFeatures)

    #Happy path tested
    def funcNormalize(self):
        """
        Convenience method which will call which ever normalization is approriate on the data.
        :return Boolean: Indicator of success (true).
        """

        if self._fIsSummed:
            return self.funcNormalizeColumnsWithSummedClades()
        else:
            return self.funcNormalizeColumnsBySum()

    #Testing Status: Light happy path testing
    def funcNormalizeColumnsBySum(self):
        """
        Normalize the data in a manner that is approrpiate for NOT summed data.
        Normalize the columns (samples) of the abundance table.
        Normalizes as a fraction of the total (number/(sum of all numbers in the column)).
        Will not act on summed tables.

        :return    Boolean:    Indicator of success. False indicates error.
        """

        if self._fIsNormalized:
#            sys.stderr.write( "This table is already normalized, did not perform new normalization request.\n" )
            return False

        if self._fIsSummed:
            sys.stderr.write( "This table has clades summed, this normalization is not appropriate. Did not perform.\n" )
            return False

        #Normalize
        for columnName in self.funcGetSampleNames():
            column = self._npaFeatureAbundance[columnName]
            columnTotal = sum(column)
            if(columnTotal > 0.0):
                column = column/columnTotal
            self._npaFeatureAbundance[columnName] = column

        #Indicate normalization has occured
        self._fIsNormalized = True

        return True

    #Happy path tested
    def funcNormalizeColumnsWithSummedClades(self):
        """
        Normalizes a summed Abundance Table.
        If this is called on a dataset which is not summed and not normalized.
        The data will be summed first and then normalized.
        If already normalized, the current normalization is kept.

        :return    Boolean:    Indicator of success. False indicates error.
        """

        if self._fIsNormalized:
#            sys.stderr.write( "This table is already normalized, did not perform new summed normalization request.\n" )
            return False

        if not self._fIsSummed:
            sys.stderr.write( "This table does not have clades summed, this normalization is not appropriate until the clades are summed. The clades are being summed now before normalization.\n" )
            self.funcSumClades()

        #Load a hash table with root data {sKey: npaAbundances}
        hashRoots = {}
        for npaRow in self._npaFeatureAbundance:

            curldAbundance = np.array(list(npaRow)[1:])
            curFeatureNameLength = len(npaRow[0].split(self._cFeatureDelimiter))
            curlRootData = hashRoots.get(npaRow[0].split(self._cFeatureDelimiter)[0])

            if not curlRootData:
                hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]] = [curFeatureNameLength, curldAbundance]
            elif curlRootData[0] > curFeatureNameLength:
                hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]] = [curFeatureNameLength, curldAbundance]

        #Normalize each feature by thier root feature
        dataMatrix = list()
        for npaRow in self._npaFeatureAbundance:

            curHashRoot = list(hashRoots[npaRow[0].split(self._cFeatureDelimiter)[0]][1])
            dataMatrix.append(tuple([npaRow[0]]+[npaRow[i+1]/curHashRoot[i] if curHashRoot[i] > 0 else 0 for i in range(len(curHashRoot))]))

        self._npaFeatureAbundance = np.array(dataMatrix,self._npaFeatureAbundance.dtype)

        #Indicate normalization has occured
        self._fIsNormalized = True

        return True
    
    def _funcRankAbundanceHelper( self, aaTodo, iRank, lRankAbundance ):
        """
        Helper method for ranking abudance which are tied.

        :params aaTodo: List of tied ranks to change to a rank.
        :type:    List of Enumerates of samples.
        :params iRank: Current Rank
        :type:    Integer
        :params lRankAbundance: Sample of abundance
        :type:    List of integers
        """

        # Subtract one from iRank (each time) to account for next loop iteration
        # Then average it with itself minus (the length of aaTodo + 1)
        dRank = ( iRank + iRank - len( aaTodo ) - 1 ) / 2.0
        for a in aaTodo:
            lRankAbundance[a[0]] = dRank

    #1 Happy path test
    def funcRankAbundance(self):
        """
        Rank abundances of features with in a sample.

        :return    AbundanceTable:    Abundance table data ranked (Features with in samples).
                              None is returned on error.
        """

        if self._npaFeatureAbundance == None:
            return None

        lsSampleNames = self.funcGetSampleNames()
        npRankAbundance = self.funcGetAbundanceCopy()
        liRanks = []
        #For each sample name get the ranks
        for sName in lsSampleNames:
            #Enumerate for order and sort abundances
            lfSample = list(enumerate(npRankAbundance[sName]))
            lfSample = sorted(lfSample, key = lambda a: a[1], reverse = True)

            # Accumulate indices until a new value is encountered to detect + handle ties
            aaTodo = []
            for i, a in enumerate( lfSample ):
                if ( not aaTodo ) or ( a[1] == aaTodo[-1][1] ):
                    aaTodo.append( a )
                else:
            # Make multiple tied ranks = average of first and last
                    self._funcRankAbundanceHelper( aaTodo, i, npRankAbundance[sName] )
                    aaTodo = [a]
            self._funcRankAbundanceHelper( aaTodo, i + 1, npRankAbundance[sName] )

        abndRanked = AbundanceTable(npaAbundance=npRankAbundance, dictMetadata=self.funcGetMetadataCopy(),
            strName= self.funcGetName() + "-Ranked",
            strLastMetadata=self.funcGetLastMetadataName(),
            cFileDelimiter=self.funcGetFileDelimiter(),
            cFeatureNameDelimiter=self.funcGetFeatureDelimiter())

        #Table is no longer normalized
        abndRanked._fIsNormalized = False
        return abndRanked

    def funcGetSampleCount(self):
        """
        Returns the sample count of the abundance table.
        """
        return len(self.funcGetSampleNames())

    #Happy Path Tested
    def funcReduceFeaturesToCladeLevel(self, iCladeLevel):
        """
        Reduce the current table to a certain clade level.

        :param    iCladeLevel:    The level of the clade to trim the features to.
        :type:    Integer    The higher the number the more clades are presevered in the consensus lineage contained in the feature name.
        :return    Boolean:    Indicator of success. False indicates error.
        """

        if iCladeLevel < 1: return False
        if not self._npaFeatureAbundance == None:
            liFeatureKeep = []
            [liFeatureKeep.append(tplFeature[0]) if (len(tplFeature[1][0].split(self.funcGetFeatureDelimiter())) <= iCladeLevel) else 0
             for tplFeature in enumerate(self._npaFeatureAbundance)]
            #Compress array
            self._npaFeatureAbundance = self._npaFeatureAbundance[liFeatureKeep,:]

            #Update filter state
            self._strCurrentFilterState += ":iCladeLevel=" + str(iCladeLevel)
            return True
        else:
            return False

    #Happy path tested
    def funcRemoveSamples(self,lsSampleNames):
        """
        Removes the samples given in the list.

        :param    lsSampleNames:    A list of string names of samples to remove.
        :type:    List of strings    Unique values
        :return Boolean: Indicator of success (True = success, no error)
        """

        #Samples to remove
        setSamples = set(lsSampleNames)

        #Get orignal sample count
        iOriginalCount  = self._iOriginalSampleCount

        #The samples to keep
        lsKeepSamples = [sSample for sSample in self.funcGetSampleNames() if not sSample in setSamples]
        #The sample to keep as boolean flags for compressing the metadata
        lfKeepSamples = [not sSample in setSamples for sSample in self.funcGetSampleNames()]
        
        #Reduce the abundance data and update
        self._npaFeatureAbundance = self._npaFeatureAbundance[[self.funcGetIDMetadataName()]+lsKeepSamples]

        #Reduce the metadata and update
        for sKey in self._dictTableMetadata:
            self._dictTableMetadata[sKey] = [value for iindex, value in enumerate(self._dictTableMetadata[sKey]) if lfKeepSamples[iindex]]

        #Update sample number count
        self._iOriginalSampleCount = len(self.funcGetSampleNames())

        return self._iOriginalSampleCount == (iOriginalCount-len(setSamples))

    #Happy path tested
    def funcRemoveSamplesByMetadata(self, sMetadata, lValuesToRemove):
        """
        Removes samples from the abundance table based on values of a metadata.
        If a metadata has any value given the associated sample is removed.

        :param    sMetadata:    ID of the metdata to check the given values.
        :type:    String    Metadata ID
        :param    lValuesToRemove:    A list of values which if equal to a metadata entry indicate to remove the associated sample.
        :type:    List of values:    List
        :return    Boolean:    Indicator of success (True = success, no error)
        """

        lsSampleNames = self.funcGetSampleNames()
        return self.funcRemoveSamples([lsSampleNames[iindex] for iindex, sValue in enumerate(self.funcGetMetadata(sMetadata)) if sValue in lValuesToRemove])

    #Happy path testing
    def funcSumClades(self):
        """
        Sums abundance data by clades indicated in the feature name (as consensus lineages).

        :return    Boolean:    Indicator of success.
                    False indicates an error.
        """

        if not self.funcIsSummed():

            #Read in the data
            #Find the header column (iCol) assumed to be 1 or 2 depending on the location of "NAME"
            #Create a list (adSeq) that will eventually hold the sum of the columns of data
            astrHeaders = iCol = None
            adSeqs = np.array([0] * len(self.funcGetSampleNames()))
            pTree = CClade( )
            aastrRaw = []

            #For each row in the npaAbundance
            #Get the feature name, feature abundances, and sum up the abudance columns
            #Keep the sum for later normalization
            #Give a tree the feature name and abundance
            for dataRow in self._npaFeatureAbundance:
                
                sFeatureName = dataRow[0]
                ldAbundances = list(dataRow)[1:]

                #Add to the sum of the columns (samples)
                adSeqs = adSeqs + np.array(list(dataRow)[1:])

                #Build tree
                pTree.get( sFeatureName.split(self._cFeatureDelimiter) ).set( ldAbundances )

            #Create tree of data
            #Input missing data
            #Fill hashFeatures with the clade name (key) and a blist of values (value) of the specified level interested.
            pTree.impute( )
            hashFeatures = {}
            pTree.freeze( hashFeatures, c_iSumAllCladeLevels, c_fOutputLeavesOnly )
            setstrFeatures = list(hashFeatures.keys( ))

            #Remove parent clades that are identical to child clades
            for strFeature, adCounts in hashFeatures.items( ):
                    astrFeature = strFeature.strip( ).split( "|" )
                    while len( astrFeature ) > 1:
                        astrFeature = astrFeature[:-1]
                        strParent = "|".join( astrFeature )
                        adParent = hashFeatures.get( strParent )
                        if adParent == adCounts:
                            del hashFeatures[strParent]
                            setstrFeatures.remove( strParent )

            #Sort features to be nice
            astrFeatures = sorted( setstrFeatures )

            #Change the hash table to an array
            dataMatrix = list()
            for sFeature in astrFeatures:
                dataMatrix.append(tuple([sFeature]+list(hashFeatures[sFeature])))
            self._npaFeatureAbundance=np.array(dataMatrix,self._npaFeatureAbundance.dtype)

            #Indicate summation has occured
            self._fIsSummed = True

        return True

    #Happy path tested
    def funcStratifyByMetadata(self, strMetadata, fWriteToFile=False):
        """
        Stratifies the AbundanceTable by the given metadata.
        Will write each stratified abundance table to file
        if fWriteToFile is True the object will used it's internally stored name as a file to write to
        if fWriteToFile is a string then it should be a directory and end with "." This will rebase the file
        and store it in a different directory but with an otherwise unchanged name.
        Note: If the metadata used for stratification has NAs, they will be segregated to thier own table and returned.

        :param    strMetadata:    Metadata ID to stratify data with.
        :type:    String    ID for a metadata.
        :param    fWriteToFile:    Indicator to write to file.
        :type:    Boolean    True indicates to write to file.
        :return    List:    List of AbundanceTables which are deep copies of the original.
                        Empty list on error.
        """

        if self._npaFeatureAbundance is None or self._dictTableMetadata is None:
            return []

        #Get unique metadata values to stratify by
        lsMetadata = self._dictTableMetadata.get(strMetadata,[])
        setValues = set(lsMetadata)
        #If there is only one metadata value then no need to stratify so return the original in the list (and write if needed)
        if len(setValues) == 0:
          return []

        retlAbundanceTables = []
        dictAbundanceBlocks = dict()
        #Given here there are multiple metadata values, continue to stratify
        lsNames = self.funcGetSampleNames()
        #Get index of values to break up
        for value in setValues:
            lfDataIndex = [sData==value for sData in lsMetadata]
            #Get abundance data for the metadata value
            #The true is added to keep the first column which should be the feature id
            npaStratfiedAbundance = self._npaFeatureAbundance[[self.funcGetIDMetadataName()]+list(np.compress(lfDataIndex,lsNames))]

            #Get metadata for the metadata value
            dictStratifiedMetadata = dict()
            for metadataType in self._dictTableMetadata:
                dictValues = self.funcGetMetadata(metadataType)
                dictStratifiedMetadata[metadataType] = np.compress(lfDataIndex,dictValues).tolist()

            #Make abundance table
            #Add abundance table to the list
            lsNamePieces = os.path.splitext(self._strOriginalName)
            objStratifiedAbundanceTable = AbundanceTable(npaAbundance=npaStratfiedAbundance, dictMetadata=dictStratifiedMetadata,
                strName=lsNamePieces[0] + "-StratBy-" + value+lsNamePieces[1],
                strLastMetadata=self.funcGetLastMetadataName(),
                cFeatureNameDelimiter=self._cFeatureDelimiter, cFileDelimiter = self._cDelimiter)
            if fWriteToFile:
                objStratifiedAbundanceTable.funcWriteToFile(lsNamePieces[0] + "-StratBy-" + value+lsNamePieces[1])
            #Append abundance table to returning list
            retlAbundanceTables.append(objStratifiedAbundanceTable)

        return retlAbundanceTables

    #Happy Path Tested
    def funcTranslateIntoMetadata(self, lsValues, sMetadataFrom, sMetadataTo, fFromPrimaryIds=True):
        """
        Takes the given data values in one metadata and translates it to values in another
        metadata of the sample samples holding the values of the first metadata
        FPrimaryIds, if true the sMetadataFrom are checked for unique values,
        If FPrimaryIds is not true, duplicate values can stop the preservation of order
        Or may cause duplication in the "to" group. This is not advised.
        if the sMetadataFrom has any duplicates the function fails and return false.

        :param    lsValues:    Values to translate.
        :type:    List    List of values.
        :param    sMetadataFrom:    The metadata the lsValues come from.
        :type:    String    ID for the metadata.
        :param    sMetadataTo:    The metadata the lsValues will be translated into keeping the samples the same.
        :type:    String    ID for the metadata.
        :param    fFromPrimaryIds:    The metadata that are in the from metadata list must be unique in each sample.
        :type:    Boolean    True indicates the metadata list should be unique in each sample. Otherwise a false will return.
        :return List:    List of new values or False on error.
        """

        #Get metadata
        lFromMetadata = self.funcGetMetadata(sMetadataFrom)
        if not lFromMetadata:
                sys.stderr.write( "Abundancetable::funcTranlateIntoMetadata. Did not receive lFromMetadata.\n" )
                return False

        lToMetadata = self.funcGetMetadata(sMetadataTo)
        if not lToMetadata:
                sys.stderr.write( "Abundancetable::funcTranlateIntoMetadata. Did not receive lToMetadata.\n" )
                return False

        #Check to see if the values are unique if indicated to do so
        if fFromPrimaryIds:
            if not len(lFromMetadata) == len(set(lFromMetadata)):
                sys.stderr.write( "Abundancetable::funcTranlateIntoMetadata. sMetadataFrom did not have unique values.\n" )
                return False

        #Translate over
        if lFromMetadata and lToMetadata:
            return [lToMetadata[iIndex] for iIndex in [lFromMetadata.index(value) for value in lsValues]]

        return False

    #Happy path tested
    def funcToArray(self):
        """
        Returns a numpy array of the current Abundance Table.
        Removes the first ID head column and the numpy array is
        Made of lists, not tuples.

        :return Numpy Array:    np.array([[float,float,...],[float,float,...],[float,float,...]])
                                None is returned on error.
        """

        if not self._npaFeatureAbundance == None:
            return np.array([list(tplRow)[1:] for tplRow in self._npaFeatureAbundance],'float')
        return None

    #Happy Path tested
    def funcWriteToFile(self, xOutputFile, cDelimiter=None, cFileType=ConstantsBreadCrumbs.c_strPCLFile):
        """
        Writes the AbundanceTable to a file strOutputFile.
        Will rewrite over a file as needed.
        Will use the cDelimiter to delimit columns if provided.

        :param    xOutputFile:    File stream or File path to write the file to.
        :type:    String    File Path
        :param    cDelimiter:    Delimiter for the output file.
        :type:    Character    If cDlimiter is not specified, the internally stored file delimiter is used.
        """

        if not xOutputFile:
            return
        # Check delimiter argument
        if not cDelimiter:
            cDelimiter = self._cDelimiter

        #  Check file type: If pcl: Write pcl file; If biom: write biom file;  If None - write pcl file
        if(cFileType == None):        
            cFileType == ConstantsBreadCrumbs.c_strPCLFile 
                
        if(cFileType == ConstantsBreadCrumbs.c_strPCLFile):
            # Write as a pcl file
            self._funcWritePCLFile(xOutputFile, cDelimiter=cDelimiter)
        elif(cFileType == ConstantsBreadCrumbs.c_strBiomFile):
            #Write as a biom  file
            self._funcWriteBiomFile(xOutputFile)
        return

    def _funcWritePCLFile(self, xOutputFile, cDelimiter=None):
        """
        Write an abundance table object as a PCL file.

        :param    xOutputFile:    File stream or File path to write the file to.
        :type:    String    File Path
        :param    cDelimiter:    Delimiter for the output file.
        :type:    Character    If cDlimiter is not specified, the internally stored file delimiter is used.
        """

        f = csv.writer(open( xOutputFile, "w" ) if isinstance(xOutputFile, str) else xOutputFile, csv.excel_tab, delimiter=cDelimiter)
        
        # Get Row metadata id info (IDs for column header, keys that line up with the ids)
        lsRowMetadataIDs, lsRowMetadataIDKeys = self.rwmtRowMetadata.funcMakeIDs() if self.rwmtRowMetadata else [[],[]]

        #Write Ids
        f.writerows([[self.funcGetIDMetadataName()]+lsRowMetadataIDs+list(self.funcGetSampleNames())])

        #Write column metadata
        lsKeys = list(set(self._dictTableMetadata.keys())-set([self.funcGetIDMetadataName(),self.funcGetLastMetadataName()]))
        lMetadataIterations = list(set(lsKeys+[self.funcGetLastMetadataName()] ))

        f.writerows([[sMetaKey]+([ConstantsBreadCrumbs.c_strEmptyDataMetadata]*len(lsRowMetadataIDs))+self.funcGetMetadata(sMetaKey) for sMetaKey in lMetadataIterations if sMetaKey != self.funcGetIDMetadataName() and not sMetaKey is None]) 

        #Write abundance
        lsOutput = list()
        curAbundance = self._npaFeatureAbundance.tolist()

        for curAbundanceRow in curAbundance:
            # Make feature metadata, padding with NA as needed
            lsMetadata = []
            for sMetadataId in lsRowMetadataIDKeys:
                lsMetadata = lsMetadata + self.rwmtRowMetadata.funGetFeatureMetadata( curAbundanceRow[0], sMetadataId )
                lsMetadata = lsMetadata + ( [ ConstantsBreadCrumbs.c_strEmptyDataMetadata ] * 
                    ( self.rwmtRowMetadata.dictMetadataIDs.get( sMetadataId, 0 ) - len( lsMetadata ) ) )
            f.writerows([[curAbundanceRow[0]]+lsMetadata+[str(curAbundanceElement) for curAbundanceElement in curAbundanceRow[1:]]])
        return

    def _funcWriteBiomFile(self, xOutputFile):
        """
        Write an abundance table object as a Biom file.
        :param    xOutputFile:    File stream or File path to write the file to.
        :type:    String    File Path    
        """
        
        #**************************
        # Get Sample Names        *
        #**************************
        lSampNames = list(self.funcGetSampleNames())

        #**************************
        # Metadata Names          *
        #**************************

        dictMetadataCopy = self.funcGetMetadataCopy()
        lMetaData = list()
        iKeysCounter = 0
        for lMetadataCopyEntry in dictMetadataCopy.items():
            iKeysCounter +=1
            sMetadataName = lMetadataCopyEntry[0]
            lMetadataEntries = lMetadataCopyEntry[1]
            iMetadataEntryCounter =  -1
            for sMetadataEntry in lMetadataEntries:
                iMetadataEntryCounter+=1
                dictMetadataNames = dict()
                dictMetadataNames[sMetadataName ] = sMetadataEntry
                if iKeysCounter == 1:
                    lMetaData.append(dictMetadataNames)
                else:
                    lMetaData[iMetadataEntryCounter][sMetadataName ] = sMetadataEntry


        #**************************
        # Observation Ids         *
        # and row metadata        *
        #**************************
        bTaxonomyInRowsFlag = False    
        if  self.rwmtRowMetadata.dictRowMetadata  is not None:
                bTaxonomyInRowsFlag = True    
            
        lObservationMetadataTable = list()

        lObservationIds = list()
        lFeatureNamesResultArray = self.funcGetFeatureNames()
        for sFeatureName  in  lFeatureNamesResultArray:
            lObservationIds.append(sFeatureName)

            if self.rwmtRowMetadata and self.rwmtRowMetadata.dictRowMetadata:
                RowMetadataEntry = self.rwmtRowMetadata.dictRowMetadata[sFeatureName][ConstantsBreadCrumbs.c_metadata_lowercase]    
                lObservationMetadataTable.append( RowMetadataEntry )

        #**************************
        # Data                    *
        #**************************
        
        lData = list()
        lAbundanceCopyResultArray = self.funcGetAbundanceCopy()

        for r in lAbundanceCopyResultArray:
            lr = list(r)
            lr.pop(0)    #Remove metadata
            lAbundanceValues = list()
            for AbundanceEntry in lr:
                flAbundanceEntry = float(AbundanceEntry)
                lAbundanceValues.append(flAbundanceEntry)
            lData.append(lAbundanceValues)
        arrData = array(lData)  #Convert list to array

        

        #**************************
        # Invoke the              *
        # biom table factory      *     
        #**************************
        if  bTaxonomyInRowsFlag == False:
            BiomTable = table_factory(arrData,
                              lSampNames,
                              lObservationIds,
                              lMetaData,
                              constructor=SparseOTUTable)
        else:                #There was metadata in the rows
            BiomTable = table_factory(arrData,
                              lSampNames,
                              lObservationIds,
                              lMetaData,
                              lObservationMetadataTable if len(lObservationMetadataTable) > 0 else None,
                              constructor=SparseOTUTable)    
      
        #**************************
        # Generate biom Output    *   
        #**************************
        f = open( xOutputFile, "w" ) if isinstance(xOutputFile, str) else xOutputFile
        f.write(BiomTable.getBiomFormatJsonString(ConstantsBreadCrumbs.c_biom_file_generated_by))
        f.close()
        return

    #Testing Status: 1 Happy path test
    @staticmethod
    def funcPairTables(strFileOne, strFileTwo, strIdentifier, cDelimiter, strOutFileOne, strOutFileTwo, lsIgnoreValues=None):
        """
        This method will read in two files and abridge both files (saved as new files)
        to just the samples in common between the two files given a common identifier.
        ***If the identifier is not unique in each data set, the first sample with the pairing id is taken so make sure the ID is unique.
        Expects the files to have the sample delimiters.

        :param    strFileOne:    Path to file one to be paired.
        :type:    String    File path.
        :param    strFileTwo:    Path to file two to be paired.
        :type:    String    File path.
        :param    strIdentifier:    Metadata ID that is used for pairing.
        :type:    String    Metadata ID.
        :param    cDelimiter:    Character delimiter to read the files.
        :type:    Character    Delimiter.
        :param    strOutFileOne:    The output file for the paired version of the first file.
        :type:    String    File path.
        :param    strOutFileTwo:    The output file for the paired version of the second file.
        :type:    String    File path.
        :param    lsIgnoreValues:    These values are ignored even if common IDs between the two files.
        :type:    List    List of strings.
        :return    Boolean:    Indicator of no errors.
                              False indicates errors.
        """

        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strFileOne)):
            sys.stderr.write( "AbundanceTable:checkRawDataFile::Error, file not valid. File:" + strFileOne + "\n" )
            return False
        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strFileTwo)):
            sys.stderr.write( "AbundanceTable:checkRawDataFile::Error, file not valid. File:"+ strFileTwo + "\n" )
            return False

        #Make file one
        #Read in file
        istm = csv.reader(open(strFileOne,'rU'), csv.excel_tab, delimiter=cDelimiter)
        lsContentsOne = [lsRow for lsRow in istm]

        #Get the file identifier for file one
        fileOneIdentifier = None
        for sLine in lsContentsOne:
            if sLine[0] == strIdentifier:
                fileOneIdentifier = sLine
                break

        #Make file two
        #Read in file
        istm = csv.reader(open(strFileTwo,'rU'), csv.excel_tab, delimiter=cDelimiter)
        lsContentsTwo = [lsRow for lsRow in istm]

        #Get the file identifier for file two
        fileTwoIdentifier = None
        for sLine in lsContentsTwo:
            if sLine[0] == strIdentifier:
                fileTwoIdentifier = sLine
                break

        #Get what is in common between the identifiers
        #And find which columns to keep in the tables based on the common elements
        setsCommonIdentifiers = set(fileOneIdentifier) & set(fileTwoIdentifier)
        if lsIgnoreValues:
            setsCommonIdentifiers = setsCommonIdentifiers - set(lsIgnoreValues)

        #Get positions of common identifiers in each data set, if the identifier is not unique in a date set just take the first index
        lfFileOneIDIndexes = [fileOneIdentifier.index(sCommonID) for sCommonID in setsCommonIdentifiers]
        lfFileTwoIDIndexes = [fileTwoIdentifier.index(sCommonID) for sCommonID in setsCommonIdentifiers]

        #Convert index list to list of boolean
        lfFileOneElements = [iIndex in lfFileOneIDIndexes for iIndex, sIdentifier in enumerate(fileOneIdentifier)]
        lfFileTwoElements = [iIndex in lfFileTwoIDIndexes for iIndex, sIdentifier in enumerate(fileTwoIdentifier)]

        #Write out file one
        ostm = csv.writer(open(strOutFileOne,'w'), csv.excel_tab, delimiter=cDelimiter)
        (ostm.writerows([np.compress(lfFileOneElements,sLine) for sLine in lsContentsOne]))

        #Write out file two
        ostm = csv.writer(open(strOutFileTwo,'w'), csv.excel_tab, delimiter=cDelimiter)
        (ostm.writerows([np.compress(lfFileTwoElements,sLine) for sLine in lsContentsTwo]))

        return True

    #Testing Status: Light happy path testing
    @staticmethod
    def funcStratifyAbundanceTableByMetadata(strInputFile = None, strDirectory = "", cDelimiter = ConstantsBreadCrumbs.c_cTab, iStratifyByRow = 1, llsGroupings = []):
        """
        Splits an abundance table into multiple abundance tables stratified by the metadata

        :param    strInputFile:    String file path to read in and stratify.
        :type:    String    File path.
        :param    strDirectory:    Output directory to write stratified files.
        :type:    String    Output directory path.
        :param    cDelimiter:    The delimiter used in the adundance file.
        :type:    Character    Delimiter.
        :param    iStratifyByRow:    The row which contains the metadata to use in stratification.
        :type:    Integer    Positive integer index.
        :param    llsGroupings:    A list of string lists where each string list holds values that are equal and should be grouped together.
                                So for example, if you wanted to group metadata "1", "2", and "3" seperately but "4" and "5" together you would
                                Give the following [["4","5"]].
                                If you know what "1" and "3" also together you would give [["1","3"],["4","5"]]
        :type    List    List of list of strings
        :return    Boolean:    Indicator of NO error.
                            False indicates an error.
        """

        #Validate parameters
        if(not ValidateData.funcIsValidFileName(strInputFile)):
            sys.stderr.write( "AbundanceTable:stratifyAbundanceTableByMetadata::Error, file not valid. File:" + strInputFile + "\n" )
            return False
        if(not ValidateData.funcIsValidStringType(cDelimiter)):
            sys.stderr.write( "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Delimiter is not a valid string/char type. Delimiter =" + cDelimiter + "\n" )
            return False
        if(not ValidateData.funcIsValidPositiveInteger(iStratifyByRow, tempZero = True) and (not ValidateData.funcIsValidString(iStratifyByRow))):
            sys.stderr.write( "AbundanceTable:stratifyAbundanceTableByMetadata::Error, Stratify by row is not a positive integer or string keyword. Row =" +
                str(iStratifyByRow) + ".\n" )
            return False

        #Get the base of the file path
        #This is dependent on the given output directory and the prefix of the file name of the input file
        #If no output file is given then the input file directory is used.
        baseFilePath = strDirectory
        lsFilePiecesExt = os.path.splitext(strInputFile)
        if baseFilePath:
            baseFilePath = baseFilePath + os.path.splitext(os.path.split(strInputFile)[1])[0]
        else:
            baseFilePath = lsFilePiecesExt[0]

        #Read in file
        istm = csv.reader(open(strInputFile,'rU'), csv.excel_tab, delimiter=cDelimiter)
        sFileContents = [lsRow for lsRow in istm]

        #Collect metadata
        metadataInformation = dict()

        #If the tempStratifyRow is by key word than find the index
        if ValidateData.funcIsValidString(iStratifyByRow):
            for iLineIndex, strLine in enumerate(sFileContents):
                if strLine[0].strip("\"") == iStratifyByRow:
                    iStratifyByRow = iLineIndex
                    break

        #Stratify by metadata row
        #Split metadata row into metadata entries
        #And put in a dictionary containing {"variable":[1,2,3,4 column index]}
        stratifyByRow = sFileContents[iStratifyByRow]
        for metaDataIndex in range(1,len(stratifyByRow)):
            metadata = stratifyByRow[metaDataIndex]
            #Put all wierd categories, none, whitespace, blank space metadata cases into one bin
            if not metadata or metadata in string.whitespace:
                metadata = "Blank"
            #Remove any extraneous formatting
            metadata = metadata.strip(string.whitespace)
            #Store processed metadata with column occurence in dictionary
            if(not metadata in metadataInformation):
                metadataInformation[metadata] = []
            metadataInformation[metadata].append(metaDataIndex)

        #For each of the groupings
        #Use the first value as the primary value which the rest of the values in the list are placed into
        #Go through the dict holding the indices and extend the list for the primary value with the secondary values
        #Then set the secondary value list to empty so that it will be ignored.
        if llsGroupings:
            for lSKeyGroups in llsGroupings:
                if len(lSKeyGroups) > 1:
                    for sGroup in lSKeyGroups[1:]:
                        if sGroup in metadataInformation:
                            metadataInformation[lSKeyGroups[0]].extend(metadataInformation[sGroup])
                            metadataInformation[sGroup] = []

        #Stratify data
        stratifiedAbundanceTables = dict()
        for tableRow in sFileContents:
            if(len(tableRow)> 1):
                for metadata in metadataInformation:
                    #[0] includes the taxa line
                    columns = metadataInformation[metadata]
                    if columns:
                        columns = [0] + columns
                        lineList = list()
                        for column in columns:
                            lineList.append(tableRow[column])
                        stratifiedAbundanceTables.setdefault(metadata,[]).append(lineList)

        #Write to file
        lsFilesWritten = []
        for metadata in stratifiedAbundanceTables:
            sOutputFile = baseFilePath+"-by-"+metadata.strip("\"")+lsFilePiecesExt[1]
            f = csv.writer(open(sOutputFile,'w'), csv.excel_tab, delimiter = cDelimiter )
            f.writerows(stratifiedAbundanceTables[metadata])
            lsFilesWritten.append(sOutputFile)

        return lsFilesWritten
        
        
                
    #*******************************************
    #* biom interface functions:               *
    #* 1. _funcBiomToStructuredArray           *
    #* 2. _funcDecodeBiomMetadata              *
    #*******************************************    
    @staticmethod
    def _funcBiomToStructuredArray(xInputFile = None):    
        """
        Reads the biom input file and builds a "BiomCommonArea"  that contains:
        1.BiomCommonArea['sLastMetadata'] - This is the name of the last Metadata (String)
        2.BiomCommonArea['BiomTaxData']- dict() - going to be used as  lcontents[0]==TaxData 
         3.BiomCommonArea['Metadata']   - dict() -  going to be used as lcontents[1]==MetaData
        4.BiomCommonArea['BiomFileInfo'] - dict() - going to be used as lcontents[2]==FileInfo (id, format:eg. Biological Observation Matrix 0.9.1) etc.
        5.BiomCommonArea['column_metadata_id'] - This is a string which is the name of the column id
          :param    xInputFile:    File path of biom file to read.
        :type:    String    File path.
        :return:   BiomCommonArea  (See description above)
        :type:    dict()        
        """    
 
        #*******************************************
        #* Build the metadata                      *
        #*******************************************
        try:

            BiomTable = load_table(xInputFile) if isinstance(xInputFile, str) else xInputFile    #Import the biom file
        except:
            print("Failure decoding biom file - please check your input biom file and rerun")
            BiomCommonArea = None
            return BiomCommonArea
 
        BiomCommonArea = dict()        
        dBugNames = list()            #Bug Names Table
        dRowsMetadata = None        #Initialize the np.array of the Rows metadata
        BiomElements  =  json.loads(BiomTable.to_json(''))
        for BiomKey, BiomValue in BiomElements.items():
        #****************************************************
        #*     Checking the different keys:  format,        *
        #*     rows, columns, date, generated_by            *
        #****************************************************
            if (BiomKey == ConstantsBreadCrumbs.c_strFormatKey  
            or BiomKey == ConstantsBreadCrumbs.c_strFormatUrl  
            or BiomKey == ConstantsBreadCrumbs.c_MatrixTtype
            or BiomKey == ConstantsBreadCrumbs.c_strTypekey
            or BiomKey == ConstantsBreadCrumbs.c_strIDKey #Same as below
            or BiomKey == ConstantsBreadCrumbs.c_GeneratedBy  #<---Need to follow up with Biom as always BiomValue = "" even though in the file has a value
            or BiomKey == ConstantsBreadCrumbs.c_strDateKey):  #Same as above
                BiomCommonArea = AbundanceTable._funcInsertKeyToCommonArea(BiomCommonArea, BiomKey, BiomValue)


            if BiomKey == ConstantsBreadCrumbs.c_rows:
                iMaxIdLen = 0 
                for iIndexRowMetaData in range(0, len(BiomValue)):
                    if ConstantsBreadCrumbs.c_id_lowercase in BiomValue[iIndexRowMetaData]:
                        sBugName = BiomValue[iIndexRowMetaData][ConstantsBreadCrumbs.c_id_lowercase]
                        dBugNames.append(sBugName)      #Post to the bug table
                        if len(sBugName) > iMaxIdLen:    #We  are calculating dynamically the length of the ID
                            iMaxIdLen  =  len(sBugName)
        
                if ConstantsBreadCrumbs.c_metadata_lowercase in BiomValue[0] and BiomValue[0][ConstantsBreadCrumbs.c_metadata_lowercase] != None :
                     dRowsMetadata = AbundanceTable._funcBiomBuildRowMetadata(BiomValue,  iMaxIdLen )

 
            if BiomKey == ConstantsBreadCrumbs.c_columns:
                BiomCommonArea = AbundanceTable._funcDecodeBiomMetadata(BiomCommonArea, BiomValue, iMaxIdLen)    #Call the subroutine to Build the metadata
 
     
        #*******************************************
        #* Build the TaxData                       *
        #*******************************************
    
        BiomTaxDataWork = list()            #Initlialize TaxData
        BiomObservations = BiomTable.iter(axis='observation')        #Invoke biom method to fetch data from the biom file
        for BiomObservationData in BiomObservations:
            sBugName = BiomObservationData[1]
            BiomTaxDataEntry = [sBugName]
            BiomObservationsValues = BiomObservationData[0]
            BiomTaxDataEntry.extend(BiomObservationsValues.tolist())
            BiomTaxDataWork.append(tuple(BiomTaxDataEntry))    
    
        BiomCommonArea[ConstantsBreadCrumbs.c_BiomTaxData] = np.array(BiomTaxDataWork,dtype=np.dtype(BiomCommonArea[ConstantsBreadCrumbs.c_Dtype]))
        BiomCommonArea[ConstantsBreadCrumbs.c_dRowsMetadata] = RowMetadata(dRowsMetadata)
        del(BiomCommonArea[ConstantsBreadCrumbs.c_Dtype])            #Not needed anymore
 
        return BiomCommonArea
    

    @staticmethod
    def _funcDecodeBiomMetadata(BiomCommonArea,  BiomValue = None,  iMaxIdLen=0 ):    
        """
        Decode the Biom Metadata and build:  
            1. BiomCommonArea['Metadata'] 
            2. BiomCommonArea['Dtype']
            3. BiomCommonArea['sLastMetadata']
            4. BiomCommonArea['column_metadata_id'] - This is a string which is the name of the column id
            These elements will be formatted and passed down the line to build the AbundanceTable
         :param    BiomValue:    The "columns" Metadata from the biom file (Contains the Metadata information)
        :type:    dict()     
        :param    iMaxIdLen:     The maximum length of a row ID
        :type:    Integer    
        :return:   BiomCommonArea 
        :type:    dict()        
        """

        BiomCommonArea[ConstantsBreadCrumbs.c_sLastMetadata] = None    #Initialize the LastMetadata element 
        BiomCommonArea['dRowsMetadata'] = None                #Initialize for cases that there is no metadata in the rows

        strLastMetadata = None
        strIDMetadata = None

        lenBiomValue = len(BiomValue)
        BiomMetadata = dict()                 
        for cntMetadata in range(0, lenBiomValue):
            BiomMetadataEntry = BiomValue[cntMetadata]
 
            for key, value in BiomMetadataEntry.items():         #Loop on the entries
                 if key == ConstantsBreadCrumbs.c_id_lowercase:        #If id - process it
                    strIDMetadata = ConstantsBreadCrumbs.c_ID
                    if  ConstantsBreadCrumbs.c_ID  not in BiomMetadata:    #If ID  not in the common area - initalize it
                        BiomMetadata[ConstantsBreadCrumbs.c_ID] = [None] * lenBiomValue #Initialize a list
                    BiomMetadata[ConstantsBreadCrumbs.c_ID][cntMetadata] = value
 
                 if  key == ConstantsBreadCrumbs.c_metadata_lowercase:        #If key = metadata
                    if  not value is None:                    #And value is not empty
                        MetadataDict = value                #Initialize a dictionary and post the values
                        for MDkey, MDvalue in MetadataDict.items():
                            MDkey = MDkey 
                            MDvalue = MDvalue 
                            
                            if  len(MDkey) > 0:        #Search for the last metadata
                                    if not strIDMetadata:
                                        strIDMetadata = MDkey
                                    BiomCommonArea[ConstantsBreadCrumbs.c_sLastMetadata] =  MDkey #Set the last Metadata
                            if  MDkey  not in BiomMetadata:
                                BiomMetadata[MDkey] = [None] * lenBiomValue
                            BiomMetadata[MDkey][cntMetadata] = MDvalue 
 

        BiomCommonArea[ConstantsBreadCrumbs.c_Metadata] = BiomMetadata
        BiomCommonArea[ConstantsBreadCrumbs.c_MetadataID] = strIDMetadata
        
        #**********************************************
        #*    Build dtype                             *
        #**********************************************

        BiomDtype = list()
        iMaxIdLen+=10 #Increase it by 10
        FirstValue = ConstantsBreadCrumbs.c_ID
        SecondValue = "U" + str(iMaxIdLen)
        BiomDtypeEntry = tuple([FirstValue, SecondValue])
        BiomDtype.append(BiomDtypeEntry)

        BiomDtype.extend([tuple( [a, ConstantsBreadCrumbs.c_f4] ) for a in BiomMetadata[ConstantsBreadCrumbs.c_ID]])
                
        BiomCommonArea[ConstantsBreadCrumbs.c_Dtype] = BiomDtype
        return BiomCommonArea

    @staticmethod
    def _funcBiomBuildRowMetadata( BiomValue, iMaxIdLen ):    
        """
        Builds the row metadata from a BIOM value

          :param    BiomValue:    BIOM Value from the BIOM JSON parsing
        :type:            Complex dict of string pairs and dicts
        :param    iMaxIdLen:    Maximum length of all the IDs
        :type:            int
        :return:        dictRowsMetadata - np Array containing the rows metadata
        :type:            {string feature id: {'metadata': {'taxonomy': [list of metadata values]}}}    
        """    
        # Build the input dict for RowMetadata from a dict of dicts from a BIOM file 
        dictRowsMetadata = dict()
        for iIndexRowMetaData in range(0, len(BiomValue)):
            dictRowsMetadata[str(BiomValue[iIndexRowMetaData][ConstantsBreadCrumbs.c_id_lowercase])] = dict()
            RowMetadataEntryFromTable = BiomValue[iIndexRowMetaData][ConstantsBreadCrumbs.c_metadata_lowercase]
            dMetadataTempDict = dict()
            for key, value in RowMetadataEntryFromTable.items():
                dMetadataTempDict[key] = value
            dictRowsMetadata[str(BiomValue[iIndexRowMetaData][ConstantsBreadCrumbs.c_id_lowercase])][ConstantsBreadCrumbs.c_metadata_lowercase] = dMetadataTempDict
        return dictRowsMetadata

    @staticmethod
    def _funcInsertKeyToCommonArea(BiomCommonArea, BiomKey, BiomValue):
        """
        Inserts the keys into the BiomCommonArea["BiomFileInfo"]
          :param    BiomCommonArea   - The common area that has been built before
        :type:    dict()
        :param    BiomKey - The current key (eg. format, date, generated by)
        :type:    str
        :param    BiomValue - The current value of the key (eg. for format: "Biological Observation Matrix 0.9.1")
        :type:    str
        :return:   BiomCommonArea  - The updated common area
        :type:    dict()        
        """    
    
        if ConstantsBreadCrumbs.c_BiomFileInfo not in BiomCommonArea:
                BiomCommonArea[ConstantsBreadCrumbs.c_BiomFileInfo] = dict()
            
        strInsertKey = BiomKey            #Set Default - But it is now always the same... (eg. URL is not: format_url -->url and others)
        PostBiomValue = BiomValue        #The default value to be posted 
        if  BiomKey == ConstantsBreadCrumbs.c_strFormatUrl:
            strInsertKey = ConstantsBreadCrumbs.c_strURLKey
            
        if  BiomKey == ConstantsBreadCrumbs.c_MatrixTtype:
            strInsertKey = ConstantsBreadCrumbs.c_strSparsityKey
            
        if  BiomKey == ConstantsBreadCrumbs.c_GeneratedBy:
            PostBiomValue = None

        if  BiomKey == ConstantsBreadCrumbs.c_strDateKey:
            PostBiomValue = None            
            
        BiomCommonArea[ConstantsBreadCrumbs.c_BiomFileInfo][strInsertKey] = PostBiomValue
        return BiomCommonArea
        
