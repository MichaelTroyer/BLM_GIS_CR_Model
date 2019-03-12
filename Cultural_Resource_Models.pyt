# -*- coding: utf-8 -*-


"""
   __________  __  ___     __  ___          __     _____            
  / ____/ __ \/  |/  /    /  |/  /___  ____/ /__  / / (_)___  ____ _
 / /   / /_/ / /|_/ /    / /|_/ / __ \/ __  / _ \/ / / / __ \/ __ `/
/ /___/ _, _/ /  / /    / /  / / /_/ / /_/ /  __/ / / / / / / /_/ / 
\____/_/ |_/_/  /_/    /_/  /_/\____/\__,_/\___/_/_/_/_/ /_/\__, /  
                                                           /____/   

Author:
Michael Troyer

Date:
8/3/18


Summary:
Implments Chi-squared Goodness-of-Fit test for arbitrary input feature classes.
Tests each input feature class for significance and reweights the input classes
according to the standard residuals.

Computes site frequency across input analysis classes within surveyed areas as 'Actuals'.
Computes expected distribution under null hypothesis (no relationship) proportional
to the representation of each class within the surveyed parts of the analysis area.
Only calculates chi2 actuals and expected values for surveyed areas; the standard residuals
are then generalized to the study area in order to identify areas of lower/higher site potential. 


Inputs:
    Analysis Area:              polygon feature class
    Existing Sites:             polygon feature class       * Options to allow sub-selection
    Existing Surveys:           polygon feature class       * Options to allow sub-selection
    Input Analysis Features:    polygon feature class       * The layers to test
    Output Workspace:           Folder
    p-value:                    Double                      * Threshold for significance
    Enable Logging:             Boolean

Returns:
    DateTime-stamped database
        Inputs - Analysis Area: Feature Dataset
            * The entire study area
            * Reweighted feature classes according to inputs
        Inputs - Surveyed Area: Feature Dataset
            * Just the surveyed areas
            * Reweighted feature classes according to analysis inputs
            * These FCs also contain all the Chi2 calculations (acres, proprtion, expected, actuals)
        Analysis_Area
        Analysis_Area Sites
        Analysis_Area Surveys
        Analysis_Results
            * The interstection of the reweighted inputs
            * Includes each individual score
            * 'SUM_SCORE' is the sum of the standard residuals for each intersected polygon
        
    Report text file            Basic operational details/metadata
    Text logfile (optional)     Extended calculations and operation details


Usage:
    * Input Analysis Features must be categorical feature classes.
    * Analysis fields must be type text or integer.
    * Site counts are based on the input site polygon centroids.
    * Does not weight inputs - site size / intra-site details are not considered. Use the filters!
    * All tests and assessments of site distribution are based on surveyed areas only.
    * Does not consider sites within a study area that are outside a survey.
    * Do not use with linear features (roads, railroads, etc..) - centroids are not relevant.
    * Will not analyze non-significant inputs or those with any expected or observed values < 5.
    * Output projection for all files will match the input analysis area.


To Do:
- Check portability - works with local installs? [scipy may not be standard lib]

- Document tool help

- Delete gdb on hard exit

- Clear memory after loops for large datasets

- Sort actuals and expecteds for reporting

- Dissallow selection if no values supplied

- Validata workspace - no .gdb

- validation testing!

 - Quantify sample strength metric

- Change sites and survey thresholds to warnings?
    May yield too few sites - see above


*** Validation test set tool
*** Eliminate classes (longest shared boundary) tool

"""


######################################################################################################################
##
## IMPORTS
##
######################################################################################################################


from __future__ import division

import arcpy
import math
import datetime
import getpass
import os
import re
import sys
import traceback
import numpy as np
import scipy.stats as stats

from collections import defaultdict


######################################################################################################################
##
## GLOBALS
##
######################################################################################################################

##---Classes----------------------------------------------------------------------------------------------------------

class pyt_log(object):
    """A custom logging class that can simultaneously write to the console - AddMessage,
       write to an optional logfile, and/or a production report..."""
       
    def __init__(self, report_path, log_path, log_active=True):
        self.report_path = report_path
        self.log_path = log_path
        self.log_active = log_active
    
    def _write_arg(self, arg, path):
        txtfile = open(path, 'a')
        if hasattr(arg, '__iter__'):
            try:
                for k, v in arg.items():
                    txtfile.write('\t{}: {}\n'.format(k, v))
            except:
                for item in arg:
                    txtfile.write('\t{}\n'.format(item))
        else:
            txtfile.write((repr(arg)+'\n'))
        txtfile.close()

    def _writer(self, msg, path, *args):
        """A writer to write the msg, and unpacked variable"""
        if os.path.exists(path):
            write_type = 'a'
        else:
            write_type = 'w'
        with open(path, write_type) as txtfile:
            txtfile.write("\n"+msg+"\n")
            txtfile.close()
            if args:
                for arg in args:
                    self._write_arg(arg, path)

    def console(self, msg):
        """Print to console only"""
        arcpy.AddMessage(msg)

    def report(self, msg):
        """Write to report only"""
        self._writer(msg, path=self.report_path)

    def logfile(self, msg, *args):
        """Write to logfile only"""
        if self.log_active:
            path = self.log_path
            self._writer(msg, path, *args)
            
    def log_report(self, msg, *args):
        """Write to logfile and report"""
        self.report(msg)
        self.logfile(msg, *args)
        
    def log_all(self, msg, *args):
        """Write to all"""
        self.console(msg)
        self.report(msg)
        self.logfile(msg, *args)


# Exterior loop exceptions (will end the program)
        
class InsufficientSurveyCoverage (BaseException): pass
class InsufficientSiteSample     (BaseException): pass

# Interior loop exceptions (will exit the individual feature class loops)

class SingleFeatureError        (BaseException): pass
class InsufficientInputCoverage (BaseException): pass
class OverlappingPolygonsError  (BaseException): pass
class NullFeatureNameError      (BaseException): pass
class ObservedExpectedMismatch  (BaseException): pass
class NullpValue                (BaseException): pass
class NonSignificance           (BaseException): pass

                    
##---Functions--------------------------------------------------------------------------------------------------------

def deleteInMemory():
    """Delete in memory tables and feature classes
       reset to original worksapce when done"""

    # get the original workspace
    orig_workspace = arcpy.env.workspace

    # Set the workspace to in_memory
    arcpy.env.workspace = "in_memory"
    # Delete all in memory feature classes
    for fc in arcpy.ListFeatureClasses():
        try:
            arcpy.Delete_management(fc)
        except: pass
    # Delete all in memory tables
    for tbl in arcpy.ListTables():
        try:
            arcpy.Delete_management(tbl)
        except: pass
    # Reset the workspace
    arcpy.env.workspace = orig_workspace


def buildWhereClauseFromList(table, field, valueList):
    """Takes a list of values and constructs a SQL WHERE
       clause to select those values within a given field and table."""
    
    # Add DBMS-specific field delimiters
    fieldDelimited = arcpy.AddFieldDelimiters(arcpy.Describe(table).path, field)
    
    # Determine field type
    fieldType = arcpy.ListFields(table, field)[0].type
    
    # Add single-quotes for string field values
    if str(fieldType) == 'String':
        valueList = ["'%s'" % value for value in valueList]
        
    # Format WHERE clause in the form of an IN statement
    whereClause = "%s IN(%s)" % (fieldDelimited, ', '.join(map(str, valueList)))
    return whereClause


def get_acres(fc):
    """Check for an acres field in fc - create if doesn't exist or flag for calculation.
       Recalculate acres and return name of acre field"""
    
    # Add ACRES field to analysis area - check if exists
    field_list = [field.name for field in arcpy.ListFields(fc) if field.name.upper() == "ACRES"]
    
    # If ACRES/Acres/acres exists in table, flag for calculation instead
    if field_list:
        acre_field = field_list[0] # select the 'acres' variant
    else:
        arcpy.AddField_management(fc, "ACRES", "DOUBLE", 15, 2)
        acre_field = "ACRES"
    arcpy.CalculateField_management(fc, acre_field, "!shape.area@ACRES!", "PYTHON_9.3")
    acres_sum = sum([row[0] for row in arcpy.da.SearchCursor(fc, acre_field)])
    return acre_field, acres_sum


def standard_residuals(observed, expected):
    return [(obs - exp) / np.sqrt(exp) for obs, exp in zip(observed, expected)]


# Unused
def adjusted_standard_residuals(observed, expected):
    '''
    O - E / sqrt(nA * nB * (1 - nA/N) * (1 - nB/N) / N)
    where nA is the row total, nB is the column total, and N is total number of cases
    '''
    standard_adjusted_residuals = []
    N = len(expected)
    nB = sum(observed)
    for obs, exp in zip(observed, expected):
        if exp == 0:
            standard_adjusted_residuals.append(0)
        nA = obs + exp
        sar = (obs - exp) / math.sqrt((nA * nB * (1 - nA / N) * (1 - nB / N) / N))
        standard_adjusted_residuals.append(sar)
    return standard_adjusted_residuals


##---Variables--------------------------------------------------------------------------------------------------------

start_time = datetime.datetime.now()
user = getpass.getuser()
err_banner = "#" * 40


##---Settings---------------------------------------------------------------------------------------------------------

arcpy.env.addOutputsToMap = False
arcpy.env.overwriteOutput = True


######################################################################################################################
##
## EXECUTION
##
######################################################################################################################


class Toolbox(object):
    def __init__(self):
        self.label = "BLM_CO_Cultural_Resource_Modelling_Toolbox"
        self.alias = "BLM CO Cultural Resource Modelling Toolbox"

        # List of tool classes associated with this toolbox
        self.tools = [Chi2_Site_Distribution]


class Chi2_Site_Distribution(object):
    def __init__(self):
        self.label = "Chi2_Site_Distribution"
        self.description = "Chi-Squared Goodness-of-Fit Test for site distribution"
        self.canRunInBackground = True 

    def getParameterInfo(self):
        """Define parameter definitions"""

######################################################################################################################
##
## PARAMETERS
##
######################################################################################################################

        # Input Analysis Area
        param00=arcpy.Parameter(
            displayName="Input Modelling Analysis Area",
            name="Input_Boundary",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            )

        # Input Surveys feature class
        param01=arcpy.Parameter(
            displayName="Input Existing Surveys Feature Class",
            name="Input_Surveys",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            )
        
        # Surveys Subselection boolean        
        param02=arcpy.Parameter(
            displayName="Survey Selection Based on Case Value",
            name="Survey_Select_Boolean",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            enabled = "False",
            )

        # Surveys Subselection field
        param03=arcpy.Parameter(
            displayName="Select Survey Case Field",
            name="Survey_Select_Field",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            enabled = "False",
            )

        # Surveys Subselection value
        param04=arcpy.Parameter(
            displayName="Select Survey Case Value",
            name="Survey_Select_Value",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            enabled = "False",
            )

        # Input Sites feature class
        param05=arcpy.Parameter(
            displayName="Input Existing Sites Feature Class",
            name="Input_Sites",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            )

        # Sites Subselection boolean        
        param06=arcpy.Parameter(
            displayName="Sites Selection Based on Case Value",
            name="Sites_Select_Boolean",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            enabled = "False",
            )

        # Sites Subselection field
        param07=arcpy.Parameter(
            displayName="Select Sites Case Field",
            name="Sites_Select_Field",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            enabled = "False",
            )

        # Sites Subselection value
        param08=arcpy.Parameter(
            displayName="Select Sites Case Value",
            name="Sites_Select_Value",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            enabled = "False",
            )

        # Input Analysis Feature Classes
        param09=arcpy.Parameter(
            displayName="Input Analysis Feature Classes",
            name="Input_Analysis_Feature_Classes",
            datatype="GPValueTable",
            parameterType="Required",
            direction="Input",
            )
        param09.columns = [['DEFeatureClass', 'Feature Class'], ['GPString', 'Field']]
        param09.filters[1].type = 'ValueList'
        param09.filters[1].list = ['Field']
                
        # Output Location
        param10=arcpy.Parameter(
            displayName="Output Workspace",
            name="Out_fGDB",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input",
            )

        # p-value
        param11=arcpy.Parameter(
            displayName="p-value",
            name="p_value",
            datatype="Double",
            parameterType="Required",
            direction="Input",
            )
        
        # Optional logging
        param12=arcpy.Parameter(
            displayName="Enable Logging",
            name="Enable_Logging",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            )

        params = [param00, param01, param02, param03, param04, param05,
                  param06, param07, param08, param09, param10, param11, param12]

        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        params[0].filter.list = ["Polygon"]
        params[1].filter.list = ["Polygon"]
        params[5].filter.list = ["Polygon"]

        # Subselect Surveys
        if params[1].value:
            params[2].enabled = "True"
        else:
            params[2].enabled = "False"
            
        if params[2].value:
            params[3].enabled = "True"
            params[3].filter.type = "ValueList"
            params[3].filter.list = [f.name for f in arcpy.Describe(params[1].valueAsText).fields
                                     if f.type in ["String", "Integer", "SmallInteger"]]
        else:
            params[3].value = ""
            params[3].enabled = "False"
            
        if params[3].value:
            params[4].enabled = "True"
            field_select = params[3].value
            arcpy.Frequency_analysis(params[1].valueAsText, "in_memory\\freq", field_select)
            
            for field in arcpy.Describe(params[1].valueAsText).fields:
                if field.name == field_select:
                    if field.type in ("Integer", "SmallInteger"):
                        where = '"{}" IS NOT NULL'.format(field_select)
                    elif field.type == "String":
                        where = '"{0}" <> \'\' and "{0}" IS NOT NULL'.format(field_select)
            
            params[4].enabled = "True"
            params[4].filter.type = "ValueList"

            featurevalueList = [row[0] for row in arcpy.da.SearchCursor("in_memory\\freq", [field_select], where)]
            featurevalueList.sort()
            
            params[4].filter.list = featurevalueList
            
        else:
            params[4].value = ""
            params[4].enabled = "False"


        # Subselect Sites
        if params[5].value:
            params[6].enabled = "True"
        else:
            params[6].enabled = "False"
            
        if params[6].value:
            params[7].enabled = "True"
            params[7].filter.type = "ValueList"
            params[7].filter.list = [f.name for f in arcpy.Describe(params[5].valueAsText).fields
                                     if f.type in ["String", "Integer", "SmallInteger"]]
        else:
            params[7].value = ""
            params[7].enabled = "False"
            
        if params[7].value:
            params[8].enabled = "True"
            field_select = params[7].value
            arcpy.Frequency_analysis(params[5].valueAsText, "in_memory\\freq", field_select)
            
            for field in arcpy.Describe(params[5].valueAsText).fields:
                if field.name == field_select:
                    if field.type in ("Integer", "SmallInteger"):
                        where = '"{}" IS NOT NULL'.format(field_select)
                    elif field.type == "String":
                        where = '"{0}" <> \'\' and "{0}" IS NOT NULL'.format(field_select)
            
            params[8].enabled = "True"
            params[8].filter.type = "ValueList"

            featurevalueList = [row[0] for row in arcpy.da.SearchCursor("in_memory\\freq", [field_select], where)]
            featurevalueList.sort()
            
            params[8].filter.list = featurevalueList
            
        else:
            params[8].value = ""
            params[8].enabled = "False"


        # Filter analysis input fields
        if params[9].value and params[9].altered:
            fields = []
            for fc, _ in params[9].values:
                fields.extend([f.name for f in arcpy.Describe(fc).fields
                               if f.type in ["String", "Integer", "SmallInteger"]])
            fields = list(set(fields))
            fields.sort()
            params[9].filters[1].list = fields

###
# Defaults
###

        if not params[11].altered:
            params[11].value = 0.05
            
        if not params[12].altered:
            params[12].value = True

    
    def updateMessages(self, params):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if params[9].altered:
            for fc, val in params[9].values:
                if val not in [f.name for f in arcpy.Describe(fc).fields]:
                    error_msg = ("Field [{}] not found in feature class [{}]".format(val, os.path.basename(str(fc))))
                    params[9].setErrorMessage(error_msg)

            fc_paths = [str(path) for path, _ in params[9].values]
                
            if len(set(fc_paths)) != len(fc_paths):
                err_mg = ("Duplicate input layers: the same layer cannnot be used more than once.")
                params[9].setErrorMessage(err_mg)

        return


######################################################################################################################
##
## EXECUTION
##
######################################################################################################################

    def execute(self, params, messages):
        """The source code of the tool."""

        try:
            # Clear memory JIC
            deleteInMemory()

            # Date/time stamp for outputs
            dt_stamp = re.sub('[^0-9]', '', str(datetime.datetime.now())[:16])

            # p-limit for stats tests
            p_limit = params[11].value
                
            # Output directory
            output_path = params[10].valueAsText

            # Input analysis area name
            analysis_area_name = os.path.splitext(os.path.basename(params[0].valueAsText))[0]
            
            # Output fGDB name and full path
            gdb_name = "{}_Chi2_Analysis_{}".format(analysis_area_name, dt_stamp)
            gdb_path = os.path.join(output_path, gdb_name) + '.gdb'
           
            # Create a geodatabase
            arcpy.CreateFileGDB_management(output_path, gdb_name, "10.0")

            # Set workspace to fGDB
            arcpy.env.workspace = gdb_name
            
            # Get input analysis area spatial reference for outputs
            spatial_ref = arcpy.Describe(params[0].value).spatialReference
            
            # Create feature datasets
            inputs_analysis_area = os.path.join(gdb_path, 'Input_Variables')
            inputs_surveyed_area = os.path.join(gdb_path, 'Chi2_Analysis')
            arcpy.CreateFeatureDataset_management(gdb_path, 'Input_Variables', spatial_ref)
            arcpy.CreateFeatureDataset_management(gdb_path, 'Chi2_Analysis', spatial_ref)
        
            # Create the logger
            log_active = params[12].value
            report_path = os.path.join(
                output_path, "{}_Chi2_Analysis_{}_Report.txt".format(
                    analysis_area_name, dt_stamp))
            logfile_path = os.path.join(
                output_path, "{}_Chi2_Analysis_{}_Logfile.txt".format(
                    analysis_area_name, dt_stamp))
            logger = pyt_log(report_path, logfile_path, log_active)
            
            # Start logging
            logger.log_all("Chi2 Site Distribution "+str(datetime.datetime.now()))
            logger.log_report("_"*120+"\n")
            logger.log_all("Running environment: Python - {}\n".format(sys.version))
            logger.log_all("User: "+user+"\n")
            logger.log_all("Output Location:")
            logger.log_all('\t'+gdb_path+'\n')

            # Log input parameters
            logger.logfile('Input Parameters:\n')
            for param in params:
                logger.logfile(param.name + ': [{}]'.format(param.valueAsText))


######################################################################################################################
##
## MAIN PROGRAM
##
######################################################################################################################

### Get the analysis area
                
            # Copy analysis area and get acres
            analysis_area = os.path.join(gdb_path, 'Analysis_Area')
            arcpy.CopyFeatures_management(params[0].valueAsText, analysis_area)
            analysis_acres_field, analysis_acres_total = get_acres(analysis_area)
            logger.log_all('\nTotal acres within analysis area: {}\n'.format(round(analysis_acres_total, 4)))

            # Get the input surveys FC as a feature layer clipped to analysis area
            arcpy.Clip_analysis(params[1].valueAsText, params[0].valueAsText, "in_memory\\in_surveys")

### Get the surveys
            
            # Get the surveys subset if given a subselection and copy to fGDB
            if params[2].value:
                case_field = params[3].valueAsText
                case_value = params[4].valueAsText
                for field in arcpy.Describe("in_memory\\in_surveys").fields:
                    if field.name == case_field:                         
                        if field.type in ("Integer", "SmallInteger"):
                            where = '"{}" = {}'.format(case_field, case_value)
                        elif field.type == "String":
                            where = "{} = '{}'".format(case_field, case_value)
                        break
                logger.log_all("Surveys subselection: WHERE [{}] = '{}'".format(case_field, case_value))
                arcpy.MakeFeatureLayer_management("in_memory\\in_surveys", "in_memory\\surveys_raw", where)
            # If no sub-selection, keep everything
            else:
                logger.log_all("No surveys subselection - using all records")
                arcpy.MakeFeatureLayer_management("in_memory\\in_surveys", "in_memory\\surveys_raw")

            # Dissolve and get survey acreage
            arcpy.Dissolve_management("in_memory\\surveys_raw", "in_memory\\surveys")
            acres_field, survey_acres_total = get_acres("in_memory\\surveys")
            
            survey_coverage = survey_acres_total / analysis_acres_total          
            logger.log_all('Survey acres within analysis area: {}'.format(round(survey_acres_total, 2)))
            logger.log_all('Survey proportion within analysis area: {}\n'.format(round(survey_coverage, 3)))

            # Enforce minimum survey coverage for analysis
            if survey_coverage < 0.05:
                raise InsufficientSurveyCoverage
            
            arcpy.CopyFeatures_management("in_memory\\surveys", os.path.join(gdb_path, 'Analysis_Surveys'))

### Get the sites
            
            # Clip sites to survey coverage
            arcpy.Clip_analysis(params[5].valueAsText, "in_memory\\surveys", "in_memory\\in_sites")

            # Get the sites subset if given a subselection and copy to fGDB
            if params[6].value:
                case_field = params[7].valueAsText
                case_value = params[8].valueAsText
                for field in arcpy.Describe("in_memory\\in_sites").fields:
                    if field.name == case_field:                         
                        if field.type in ("Integer", "SmallInteger"):
                            where = '"{}" = {}'.format(case_field, case_value)
                        elif field.type == "String":
                            where = "{} = '{}'".format(case_field, case_value)
                        break
                logger.log_all("Sites subselection: where [{}] = '{}'".format(case_field, case_value))
                arcpy.MakeFeatureLayer_management("in_memory\\in_sites", "in_memory\\raw_sites", where)
            # If no sub-selection, keep everything
            else:
                logger.log_all("No sites subselection - using all records")
                arcpy.MakeFeatureLayer_management("in_memory\\in_sites", "in_memory\\raw_sites")
                
            site_count = int(arcpy.GetCount_management("in_memory\\raw_sites").getOutput(0))
            site_density = round(site_count/survey_acres_total, 4)
            acres_per_site = round(1/site_density, 2)
            logger.log_all('Sites identified for analysis: {}'.format(site_count))
            logger.log_all('Site density in surveyed areas (sites/acre): {}'.format(site_density))
            logger.log_all('Approximately 1 site every {} acres\n'.format(acres_per_site))
            
            if site_count < 30:
                raise InsufficientSiteSample
            
            arcpy.MakeFeatureLayer_management("in_memory\\raw_sites", "in_memory\\sites")
            arcpy.CopyFeatures_management("in_memory\\sites", os.path.join(gdb_path, 'Analysis_Sites'))

#### Main loop

            # Inputs 
            input_fc_paths_and_fields = params[9].value
               
            # Main data structure
            analysis_data = defaultdict(dict)
            
            for input_fc_path, field in input_fc_paths_and_fields:
                try:
                    in_name = os.path.basename(str(input_fc_path))
                    #drop '.shp' if input shapefile
                    in_name = os.path.splitext(in_name)[0]
                    logger.log_all('[{}]'.format(in_name))
                    logger.log_all(''.center(100, '-'))
                    path_to_clipped = os.path.join(inputs_surveyed_area, '{}_chi2'.format(in_name))
                    path_to_source = os.path.join(inputs_analysis_area, in_name)

                    analysis_data[in_name]['Source_Path'] = path_to_source
                    analysis_data[in_name]['Clipped Path'] = path_to_clipped
                    analysis_data[in_name]['Field'] = field

                    # Clip inputs by analysis area
                    arcpy.Clip_analysis(input_fc_path, analysis_area, "in_memory\\fc_clip")
                    arcpy.Dissolve_management("in_memory\\fc_clip", path_to_source, field)
                    
                    # Clip analysis area clips by existing survey coverage
                    arcpy.Clip_analysis(path_to_source, "in_memory\\surveys", "in_memory\\clip")
                    arcpy.Dissolve_management("in_memory\\clip", path_to_clipped, field)
                    
                    # Make sure both have updated acres - have same acres field
                    acres_field, analysis_acres = get_acres(path_to_source)
                    acres_field, surveyed_acres = get_acres(path_to_clipped)
                    
                    # Check more than one class represented
                    if int(arcpy.GetCount_management(path_to_clipped).getOutput(0)) <= 1:
                        raise SingleFeatureError
                    
                    # Make sure input feature class is fully (mostly) represented in surveyed area
                    # Must vary by less than 1/1000 of the surveyed acreage
                    survey_prop = surveyed_acres / survey_acres_total
                    logger.logfile('Feature Class Coverage in surveyed areas: [{}]'.format(survey_prop))
                    
                    if survey_prop < .999:  # Coverage must be incomplete
                        raise InsufficientInputCoverage
                    if survey_prop > 1.001:  # Surveyed area exceeds analysis area - must be double counting
                        raise OverlappingPolygonsError

                    # Each class is now a single row
                    # Calculate acres, proportion, expected sites count, and actual site count
                    class_acreages = {c[0]: c[1] for c in arcpy.da.SearchCursor(path_to_clipped, [field, acres_field])}
                                      
                    # Make sure no null class names
                    if not all(class_acreages.keys()):  # Have a weird value in dissolve field - not meaningful
                        raise NullFeatureNameError
                    
                    total_class_acres = sum(class_acreages.values())
                    class_proportions = {ca_name: ca / total_class_acres for ca_name, ca in class_acreages.items()}
                                          
                    logger.logfile('Total Class Acres - [{}]'.format(in_name), total_class_acres)
                    logger.logfile('Class Acreages - [{}]'.format(in_name), class_acreages)
                    logger.logfile('Class Proportions - [{}]'.format(in_name), class_proportions)
                    analysis_data[in_name]['Class_Acreages'] = class_acreages
                    analysis_data[in_name]['Class_Proportions'] = class_proportions
                                       
                    # Site count analysis
                    arcpy.FeatureToPoint_management("in_memory\\sites", "in_memory\\pts", "INSIDE")
                    arcpy.TabulateIntersection_analysis(path_to_clipped, field, "in_memory\\pts", "in_memory\\int")
                    n = int(arcpy.GetCount_management("in_memory\\sites").getOutput(0))

                    # Calculate the actual values
                    actual_values = {row[0]: row[1]
                                     for row in arcpy.da.SearchCursor("in_memory\\int", [field, "PNT_COUNT"])}
                    sum_actual = sum(actual_values.values())

                    # Calculate the expected values
                    expected_values = {val: int(round(n * (acres / total_class_acres)))
                                       for val, acres in class_acreages.items()}
                    sum_expected = sum(expected_values.values())
                    
                    # make sure all classes represented in actuals (calculated from intersect)
                    # Intesect may not have caught all the classes: i.e. where class = 0
                    for k, v in expected_values.items():
                        if k not in actual_values.keys():
                            actual_values[k] = 0
                            
                    logger.logfile('Actual values - [{}]'.format(in_name), actual_values)
                    analysis_data[in_name]['Actuals'] = actual_values
                    analysis_data[in_name]['Sum_Actuals'] = sum_actual                
                    logger.logfile('Expected values - [{}]'.format(in_name), expected_values)
                    analysis_data[in_name]['Expected'] = expected_values
                    analysis_data[in_name]['Sum_Expected'] = sum_expected
                    
                    logger.logfile('Actual sum {} : Expected sum {}'.format(sum_actual, sum_expected))
                    
                    # Make sure actual and expected vary by less than 1/1000 n - allowing a min of 2
                    site_count_threshold = .001 * n if .001 * n >= 2 else 2
                    if abs(sum_actual - sum_expected) > site_count_threshold:
                        raise ObservedExpectedMismatch

                    # Restructure the data into rows: name, acres, proprtion, expected value, actual value - sort by name
                    class_names = expected_values.keys()
                    data_table = [
                        (class_name,
                         analysis_data[in_name]['Class_Acreages'][class_name],
                         analysis_data[in_name]['Class_Proportions'][class_name],
                         analysis_data[in_name]['Expected'][class_name],
                         analysis_data[in_name]['Actuals'][class_name],
                         )
                        for class_name in class_names]

                    # Restructure the columns as rows
                    chi2_names, class_acres, class_proportions, expects, actuals = zip(*data_table)

                    # any exp == 0 or obs == 0?
                    ev_zeros = [(k,v) for k, v in expected_values.items() if v < 5]
                    av_zeros = [(k,v) for k, v in actual_values.items() if v < 5]
                    
                    if any(ev_zeros) or any(av_zeros):
                        # This is the only error scenario we want to continue past so that we can
                        # update the clipped FC table with proportion and count data so the user doesn't have
                        # to review the log file to figure out which classes need to be collapsed.
                        msg = ('[{}] Yields expected or observed values of less than 5.\n'
                               'The Chi-squared test statistic is not appropriate in this case.\n'
                               'Consider collapsing categories:\n'
                               '===============================\n'
                               '{}\n'
                               '==============================='
                               ''.format(in_name, '\n'.join(
                                   ['Expected: {} == {}'.format(k, v) for k, v in ev_zeros] +
                                   ['Actual: {} == {}'.format(k, v) for k, v in av_zeros]
                                   )))
                        logger.log_all(msg)
                        # All dummy data
                        p, std_res_raw = None, [None] * len(chi2_names)
                        std_res = {c_name: None for c_name in chi2_names}
                        analysis_data[in_name]['p_value'] = None
                        analysis_data[in_name]['Std_residuals'] = std_res
                        analysis_data[in_name]['Reweighted'] = False
                        
                    else:                
                        # Calculate the Chi-quare test statistic
                        chi_test, p = stats.chisquare(actuals, expects)
                        p = p if p > 0.00001 else 0.0  # don't bother with really small floating points
                        std_res_raw = standard_residuals(actuals, expects)
                        std_res = {c_name: res for c_name, res in zip(chi2_names, std_res_raw)}
                        
                        logger.log_all('[{}] Chi-squared test statistic: {}'.format(in_name, chi_test))
                        logger.log_all('[{}] Chi-squared p-value: {}'.format(in_name, p))
                        logger.log_all('[{}] Chi-squared df: {}'.format(in_name, len(actuals) - 1))
                        logger.logfile('\nStandard residuals: {}'.format(in_name), std_res)
                        analysis_data[in_name]['Chi2'] = chi_test
                        analysis_data[in_name]['p_value'] = p
                        analysis_data[in_name]['Std_residuals'] = std_res

                    results_rows = zip(chi2_names, class_acres, class_proportions, expects, actuals, std_res_raw)

                    add_fields = [
                                  ('Feature_Acres', 'Double'),
                                  ('Feature_Proportion', 'Double'),
                                  ('Expected_Site_Count', 'Integer'),
                                  ('Observed_Site_Count', 'Integer'),
                                  ('Standard_Residuals', 'Double'),
                                  ]
                    for field_name, field_type in add_fields:
                        arcpy.AddField_management(path_to_clipped, field_name, field_type)

                    # Write to table
                    with arcpy.da.UpdateCursor(
                        path_to_clipped, [field] + [af[0] for af in add_fields]) as cursor:
                        for row in cursor:
                            cl_name = row[0]
                            row[1] = analysis_data[in_name]['Class_Acreages'][cl_name]
                            row[2] = analysis_data[in_name]['Class_Proportions'][cl_name]
                            row[3] = analysis_data[in_name]['Expected'][cl_name]
                            row[4] = analysis_data[in_name]['Actuals'][cl_name]
                            row[5] = analysis_data[in_name]['Std_residuals'][cl_name]
                            cursor.updateRow(row)

                    # Filter out non-significant feature classes
                    if p is None:  # If there were zeros in actuals or expected - now exit the loop
                        raise NullpValue
                    if p > p_limit:
                        raise NonSignificance

                    # Add weight fields
                    weight_field = "Weight_{}".format(in_name)
                    analysis_data[in_name]['Weight_Field'] = weight_field
                    
                    # Add weight field and weight original and clipped FCs
                    for path_to_weight in [path_to_clipped, path_to_source]:
                        arcpy.AddField_management(path_to_weight, weight_field, "Double")
                        with arcpy.da.UpdateCursor(path_to_weight, [field, weight_field]) as cursor:
                            for row in cursor:
                                try:
                                    row[1] = analysis_data[in_name]['Std_residuals'][row[0]]
                                except KeyError:
                                    row[1] = 0  
                                cursor.updateRow(row)
                    analysis_data[in_name]['Reweighted'] = True

#### Interior loop exceptions

                except SingleFeatureError:
                    err_msg = (
                        '[{}] Cannot be analyzed.\n'
                        'Dissolving on field [{}] yields a single variable.\n'
                        'The Chi-squared test statistic is not appropriate in this case.\n'
                        'Consider expanding your categories.'
                        ''.format(in_name, field))
                    logger.log_all(err_msg)
                    continue

                except InsufficientInputCoverage:
                    err_msg = (
                        '[{}] Cannot be analyzed due to incomplete coverage.\n'
                        'Verify that your inputs cover the entire analysis area.'
                        ''.format(in_name))
                    logger.log_all(err_msg)
                    continue

                except OverlappingPolygonsError:
                    err_msg = (
                        '[{}] Cannot be analyzed due to excess acre calculations.\n'
                        'The sum area of the input exceeds the size of the anaysis area.\n'
                        'Verify that your inputs do not have overlapping polygons.'
                        ''.format(in_name))
                    logger.log_all(err_msg)
                    continue

                except NullFeatureNameError:
                    err_msg = (
                        '[{}] Cannot be analyzed due to a null value in the analysis field.\n'
                        'Verify that your inputs do not have empty or null values in the selected field.'
                        ''.format(in_name))
                    logger.log_all(err_msg)
                    continue

                except ObservedExpectedMismatch:
                    # Ideally these should all be caught by IncompleteCoverage and OerlappingInputs
                    err_msg = (
                        '[{}] Cannot be analyzed due to mismatch in the total observed and expected site counts.\n'
                        'The Chi-squared test statistic is not appropriate in this case.\n'
                        'This is likely due to an area calculation error in the input.\n'
                        'Verify that your inputs cover the entire analysis area and do not have overlapping polygons.'
                        ''.format(in_name))
                    logger.log_all(err_msg)
                    continue

                except NullpValue:
                    # Do nothing - result was previously documented
                    continue

                except NonSignificance:
                    err_msg = (
                        'Chi-squared (p={}) is not significant at the {} level.\n'
                        '{} Features will not be considered in final output.'
                        ''.format(p, p_limit, in_name))
                    logger.log_all(err_msg)
                    continue
                    
                except Exception as e:
                    # Something bad.. hard quit, handle in main block
                    raise e

                finally:
                    logger.log_all(''.center(100, '-'))
                    logger.log_all('\n')

### End main loop
### Compile results

            logger.console('Compiling Results..\n')
            
            # Intersect all reweighted source layers
            final_output = os.path.join(gdb_path, 'Analysis_Results')

            # List of dictionaries - name, data
            reweighted = [d for n, d in analysis_data.items() if d.get('Reweighted', False)]
            
            if len(reweighted) > 1:
                weight_fields = [item['Weight_Field'] for item in reweighted]
                source_fields = [item['Field'] for item in reweighted]
                intersect_paths = [item['Source_Path'] for item in reweighted]
                arcpy.Intersect_analysis(intersect_paths, final_output)
                
                # Remove unneeded fields
                keep_fields = weight_fields + source_fields
                drop_fields = [f for f in arcpy.ListFields(final_output, "*") if f.name not in keep_fields]
                for drop_field in drop_fields:
                    try:
                        arcpy.DeleteField_management(final_output, drop_field.name)
                    except:
                        pass

                # calculate acres
                get_acres(final_output)
                
                # Sum weights fields and calculate mean weight
                arcpy.AddField_management(final_output, "Sum_Std_Res", 'Double', 10, 4)
                arcpy.AddField_management(final_output, "Mean_Std_Res", 'Double', 10, 4)
                with arcpy.da.UpdateCursor(final_output, weight_fields + ["Sum_Std_Res", "Mean_Std_Res"]) as cur:
                    for row in cur:
                        row[-2] = sum(row[:-2])
                        row[-1] = sum(row[:-2]) / len(row[:-2])
                        cur.updateRow(row)
            else:
                logger.log_all('Not enough feature classes appropriate for summary analysis..\n')
            logger.log_all('Complete..\n')

######################################################################################################################
##
## EXCEPTIONS (Outer loop)
##
######################################################################################################################

        except InsufficientSurveyCoverage:
            msg = (
                "Insufficient survey coverage in analysis area.\n"
<<<<<<< HEAD
                "Model requires a minimum of 5 percent survey coverage for analysis.\n"
=======
                "Model requires a minimum of 10 percent survey coverage for analysis.\n"
>>>>>>> 3324790478fa542e9ca4323dda7b0c1f6a7a195c
                "[exiting program]")
            logger.log_all("{}\n\n{}\n\n{}".format(err_banner, msg, err_banner))
            arcpy.AddError(msg)
            
        except InsufficientSiteSample:
            msg = (
                "Too few sites in analysis area.\n"
                "Model requires a minimum of 30 sites for analysis.\n"
                "[exiting program]")
            logger.log_all("{}\n\n{}\n\n{}".format(err_banner, msg, err_banner))
            arcpy.AddError(msg)
            
        except:
            # [Rails]  <-- -->  [Car]
            err = str(traceback.format_exc())
            try:
                logger.logfile(err)
            except:
                pass
            arcpy.AddError(err)


######################################################################################################################
##
## CLEAN-UP
##
######################################################################################################################

        finally:
            deleteInMemory()


######################################################################################################################