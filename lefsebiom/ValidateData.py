"""
Author: Timothy Tickle
Description: Validate Data containing methods for testing variables.
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

#Import local code
from types import *
import decimal
import os
import re
import string

class ValidateData:

    #Tested 5
    @staticmethod
    def funcIsValidBoolean(parameterValue):
        """
        Validates a parameter as a valid boolean.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid boolean.
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is BooleanType:
            return False
        return True

    #Tested 6
    @staticmethod
    def funcIsTrue(parameterValue):
        """
        Validates a parameter as true.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is True.
        :type	Boolean
        """

        if(ValidateData.funcIsValidBoolean(parameterValue)):
            if(parameterValue == True):
                return True
        return False

    #Tested 6
    @staticmethod
    def funcIsFalse(parameterValue):
        """
        Validates a parameter as false.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is False.
        :type	Boolean
        """

        if(ValidateData.funcIsValidBoolean(parameterValue)):
            if(parameterValue == False):
                return True
        return False

    #Tested 5
    @staticmethod
    def funcIsValidInteger(parameterValue):
        """
        Validates a parameter as an integer.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is an integer.
        :type	Boolean
        """

        #Check to make sure it is not null
        if (parameterValue == None):
            return False

        #Check to make sure it is an integer
        if not type(parameterValue) is IntType:
            return False

        return True

    #Tested 5
    @staticmethod
    def funcIsValidPositiveInteger(parameterValue, tempZero = False):
        """
        Validates a parameter as false.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :param	tempZero:	Allows one to set what the value for zero should return.
        :type	Boolean	The return value for zero.
        :return	Boolean:	True indicates the parameter is a positive integer.
        :type	Boolean
        """

        #Check to make sure it is not null
        if not ValidateData.funcIsValidInteger(parameterValue):
            return False

        #Check to see it is positive
        if (parameterValue < 0):
            return False

        #Check for zero value
        if(parameterValue == 0):
            return tempZero
        return True

    #Tested 14
    @staticmethod
    def funcIsValidNumeric(parameterValue):
        """
        Validates a parameter as an integer.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a numeric.
        :type	Boolean
        """

        #Check to make sure it is not null
        if (parameterValue == None):
            return False
        #Check to make sure it is an integer
        if((type(parameterValue) == IntType)or(type(parameterValue) == LongType)or(type(parameterValue) == FloatType)or(type(parameterValue) == ComplexType)or(str(type(parameterValue)) == "<type 'numpy.float64'>")):
            if(not type(parameterValue) == BooleanType):
                return True
        return False

    #Tested 5
    @staticmethod
    def funcIsValidStringType(parameterValue):
        """
        Validates a parameter as a string. This allows the string to be blank or empty.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a string type.
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is StringType:
            return False

        return True

    #Tested 5
    @staticmethod
    def funcIsValidString(parameterValue):
        """
        Validates a parameter as a string. Does NOT allow string to be blank or empty.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a string.
        :type	Boolean
        """

        #Type check
        if not ValidateData.funcIsValidStringType(parameterValue):
            return False

        #Check to see it is not blank
        if parameterValue.strip() == "":
            return False
        return True

    @staticmethod
    def funcIsValidStringInt(parameterValue):
        """
        Validates a parameter that is a string as a format which is an integer.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        """

        #Type string check
        if not ValidateData.funcIsValidStringType(parameterValue):
            return False

        #Check to see if the string can be converted to an integer
        try:
            int(parameterValue)
        except:
            return False
        return True

    @staticmethod
    def funcIsValidStringFloat(parameterValue):
        """
        Validates a parameter that is a string as a format which is a numeric.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        """

        #Type string check
        if not ValidateData.funcIsValidStringType(parameterValue):
            return False

        #Check to see if the string can be converted to a double
        try:
            float(parameterValue)
        except:
            return False
        return True

    #Tested 6
    @staticmethod
    def funcIsValidFormatString(parameterValue):
        """
        Validates a parameter as a valid format string.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        lettersValid = False
        if ValidateData.funcIsValidString(parameterValue):
            validChars = "BbcdfHhIiLlPpsx0123456789"
            for letter in parameterValue:
                lettersValid = letter in validChars
                if(not lettersValid):
                    break
        return lettersValid

    #Tested 5
    @staticmethod
    def funcIsValidChar(parameterValue):
        """
        Validates a parameter as a valid character.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        return ValidateData.funcIsValidString(parameterValue)

    #Tested 13
    @staticmethod
    def funcIsValidPositiveNumberChar(parameterValue):
        """
        Validates a parameter as a valid character representing a number.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        #Check to make sure is a valid string
        if not ValidateData.funcIsValidString(parameterValue):
            return False

        #Try to convert to decimal
        try:
            decimalConversion = decimal.Decimal(parameterValue)
            if decimalConversion < 0:
                return False
        except:
            return False
        return True

    #Tested 9
    @staticmethod
    def funcIsValidFlagChar(parameterValue):
        """
        Validates a parameter as a valid character representing a boolean.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        if parameterValue == '0' or parameterValue == "0" or parameterValue == '1' or parameterValue == "1":
            return True
        return False

    #Tested 15
    @staticmethod
    def funcIsValidBoundedIntegerChar(parameterValue, iValueOne, iValueTwo):
        """
        Validates a parameter as a valid characater that represents an integer inclusively bounded by two given values.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :param	iValueOne:	One bound for the value.
        :type	Integer
        :param	iValueTwo:	The other bound for the data.
        :type	Integer
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        #Check to make sure is a valid string
        if not ValidateData.funcIsValidString(parameterValue):
            return False

        #Check to make sure is a valid integer
        if not ValidateData.funcIsValidInteger(iValueOne):
            return False

        #Check to make sure is a valid integer
        if not ValidateData.funcIsValidInteger(iValueTwo):
            return False

        #Try to convert to decimal
        try:
            intConversion = int(parameterValue)
            if(iValueOne < iValueTwo):
                if ((intConversion >= iValueOne) and (intConversion <= iValueTwo)):
                    return True
                return False
            if(iValueTwo < iValueOne):
                if ((intConversion >= iValueTwo) and (intConversion <= iValueOne)):
                    return True
                return False
            if(iValueOne == iValueTwo):
                if (intConversion == iValueOne):
                    return True
                return False
        except:
            return False

    #Tested 9
    @staticmethod
    def funcIsValidList(parameterValue):
        """
        Validates a parameter as a list.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a list
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a list
        if not type(parameterValue) is ListType:
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in range(0,listSize):
            if parameterValue[i] == None:
                return False
            if type(parameterValue[i]) is ListType:
                if ValidateData.funcIsValidList(parameterValue[i]) == False:
                    return False
        return True

    #Tested 9
    @staticmethod
    def funcIsValidTuple(parameterValue):
        """
        Validates a parameter as a tuple.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a tuple
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is TupleType:
            return False

        #Check elements
        tupleSize = len(parameterValue)
        for i in range(0,tupleSize):
            if parameterValue[i] == None:
                return False
            if type(parameterValue[i]) is TupleType:
                if ValidateData.funcIsValidTuple(parameterValue[i]) == False:
                    return False
        return True

    #Tested 7
    @staticmethod
    def funcIsValidNumericList(parameterValue):
        """
        Validates a parameter as a list of numeric values.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a list of numeric values.
        :type	Boolean
        """

        #Check is valid list
        if(not ValidateData.funcIsValidList(parameterValue)):
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in range(0,listSize):
            if(not ValidateData.funcIsValidNumeric(parameterValue[i])):
                return False
        return True

    #Tested 7
    @staticmethod
    def funcIsValidStringList(parameterValue):
        """
        Validates a parameter as a list of string values.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a list of string values.
        :type	Boolean
        """

        #Check is valid list
        if(not ValidateData.funcIsValidList(parameterValue)):
            return False

        #Check elements
        listSize = len(parameterValue)
        for i in range(0,listSize):
            if(not ValidateData.funcIsValidString(parameterValue[i])):
                return False
        return True

    #Tested 4
    @staticmethod
    def funcIsValidNPArray(parameterValue):
        """
        Validates a parameter as a numpy array.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a numpy array.
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a structure array
        if not str(type(parameterValue)) == "<type 'numpy.ndarray'>":
            return False

        return True

    #Tested 9
    @staticmethod
    def funcIsValidDictionary(parameterValue):
        """
        Validates a parameter as a dictionary.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a dictionary.
        :type	Boolean
        """

        #Check to make sure it is not null
        if parameterValue == None:
            return False

        #Check to make sure it is a string
        if not type(parameterValue) is DictType:
            return False

        #Check key elements
        keyList = list(parameterValue.keys())
        keyListSize = len(keyList)
        for i in range(0,keyListSize):
            if keyList[i] == None:
                return False
            if type(keyList[i]) is ListType:
                if validateData.funcIsValidList(keyList[i]) == False:
                    return False

        #Check key elements
        itemList = list(parameterValue.values())
        itemListSize = len(itemList)

        for i in range(0,itemListSize):
            if itemList[i] == None:
                return False
            if type(itemList[i]) is ListType:
                if ValidateData.funcIsValidList(itemList[i]) == False:
                    return False
        return True

    #Tested 18
    @staticmethod
    def funcIsValidDNASequence(parameterValue):
        """
        Validates a parameter as a valid DNA sequence.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        if ValidateData.funcIsValidString(parameterValue):
            expression = re.compile(r'[^atcgATCG]')
            if not None == expression.search(parameterValue):
                return False
            return True
        return False

    #Tested 15
    @staticmethod
    def funcIsValidNucleotideBase(parameterValue):
        """
        Validates a parameter as a character which is a valid nucleotide representation.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        if (ValidateData.funcIsValidDNASequence(parameterValue) or (parameterValue == 'u') or (parameterValue == "U")):
            if (len(parameterValue) == 1):
                return True
        return False

    #Testing 4
    @staticmethod
    def funcIsValidFileName(parameterValue):
        """
        Validates a parameter as a valid file name.

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid file path.
        :type	Boolean
        """

        if parameterValue is None:
            return False
        elif(ValidateData.funcIsValidString(parameterValue)):
            return os.path.exists(parameterValue)
        return False

    #Tested 5
    @staticmethod
    def funcIsValidClass(parameterValue, strCorrectName):
        """
        Validates a parameter as a valid class (of a specifc type given by name).

        :param	parameterValue:	Value to be evaluated.
        :type	Unknown
        :param	strCorrectName:	Name of te class the parameter should be.
        :type	Unknown
        :return	Boolean:	True indicates the parameter is a valid value.
        :type	Boolean
        """

        if(parameterValue==None):
            return False
        if not ValidateData.funcIsValidString(strCorrectName):
            return False
        classType = type(parameterValue).__name__
        if(classType == strCorrectName):
            return True
        if(classType == 'instance'):
            if(parameterValue.__class__.__name__==strCorrectName):
                return True
            else:
                return False
        return False
