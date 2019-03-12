# -*- coding: utf-8 -*-


"""
Author:
Michael Troyer

Date:
8/3/18


Summary:
Clips input feature classes and rasters to analysis area and copies to datetime
stamped geodatabase. Dissolves feature classes on field_name.


Inputs:
    Inut Clip Boundary:         path to polygon feature class
    Input Analysis Features:    multiValue List: path to feature class | field_name
    Input Analysis Rasters:     list of rasters
    Output Workspace:           folder

Returns:
    DateTime-stamped database
        clip_boundary
        Feature Class clips
        Raster clips

Usage:
    * Clip and process data for use with modeling tools.


To Do:

"""

import arcpy
import datetime
import os
import re
import sys
import traceback

                    
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


##---Variables--------------------------------------------------------------------------------------------------------

start_time = datetime.datetime.now()


##---Settings---------------------------------------------------------------------------------------------------------

arcpy.env.addOutputsToMap = False
arcpy.env.overwriteOutput = True


##--Execution---------------------------------------------------------------------------------------------------------

class Toolbox(object):
    def __init__(self):
        self.label = "Toolbox"
        self.alias = "Toolbox"
        self.tools = [Batch_Clip_Features_And_Rasters]


class Batch_Clip_Features_And_Rasters(object):
    def __init__(self):
        self.label = "Batch_Clip_Features_And_Rasters"
        self.description = "Batch Clip Features and Rasters"
        self.canRunInBackground = True 

    def getParameterInfo(self):
        """Define parameter definitions"""

        param00=arcpy.Parameter(
            displayName="Input Clip Boundary",
            name="Input_Boundary",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            )

        param01=arcpy.Parameter(
            displayName="Input Feature Classes",
            name="Input_Feature_Classes",
            datatype="GPValueTable",
            parameterType="Optional",
            direction="Input",
            )
        param01.columns = [['DEFeatureClass', 'Feature Class'], ['GPString', 'Field']]
        param01.filters[1].type = 'ValueList'
        param01.filters[1].list = ['Field']
            
        param02=arcpy.Parameter(
            displayName="Input Rasters",
            name="Input_Rasters",
            datatype="DERasterDataset",
            parameterType="Optional",
            direction="Input",
            multiValue=True,
            )
        # Output Location
        param03=arcpy.Parameter(
            displayName="Output Folder Location (will create a file geodatabase in this folder)",
            name="Out_fGDB",
            datatype="DEWOrkspace",
            parameterType="Required",
            direction="Input",
            )

        params = [param00, param01, param02, param03,]

        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True


    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        # Filter input data types
        params[0].filter.list = ["Polygon"]
        # params[1].filter.list = ["Polygon"]

        # Filter input field types
        if params[1].value and params[1].altered:
            fields = []
            for fc, _ in params[1].values:
                fields.extend([f.name for f in arcpy.Describe(fc).fields
                               if f.type in ["String", "Integer", "SmallInteger"]])
            fields = list(set(fields))
            fields.sort()
            params[1].filters[1].list = fields

   
    def updateMessages(self, params):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        
        if params[1].altered:
            for fc, val in params[1].values:
                if val not in [f.name for f in arcpy.Describe(fc).fields]:
                    error_msg = ("Field [{}] not found in feature class [{}]".format(val, os.path.basename(str(fc))))
                    params[1].setErrorMessage(error_msg)

            fc_paths = [str(path) for path, _ in params[1].values]
                
            if len(set(fc_paths)) != len(fc_paths):
                err_mg = ("Duplicate input layers: the same layer cannnot be used more than once.")
                params[9].setErrorMessage(err_mg)

        return


    def execute(self, params, messages):
        """The source code of the tool."""

        arcpy.AddMessage("Running environment: Python - {}\n".format(sys.version))
        try:
            # Clear memory JIC
            deleteInMemory()

            # Date/time stamp for outputs
            dt_stamp = re.sub('[^0-9]', '', str(datetime.datetime.now())[:16])

            # Output directory
            output_path = params[3].valueAsText
           
            # Output fGDB name and full path
            gdb_name = "DataPrep_BatchClips_{}".format(dt_stamp)
            gdb_path = os.path.join(output_path, gdb_name) + '.gdb'
           
            # Create a geodatabase
            arcpy.CreateFileGDB_management(output_path, gdb_name, "10.0")

            # Set workspace to fGDB
            arcpy.env.workspace = gdb_name
            
            arcpy.AddMessage("Output Geodatabase Location:")
            arcpy.AddMessage('\t'+gdb_path+'\n')
                
            # Copy clip boundary
            clip_boundary = os.path.join(gdb_path, 'Clip_Boundary')
            arcpy.CopyFeatures_management(params[0].valueAsText, clip_boundary)

            # Clip inputs feature classes
            input_fc_paths_and_fields = params[1].value
            
            for input_fc_path, field in input_fc_paths_and_fields:
                try:
                    fc_name = os.path.splitext(os.path.basename(str(input_fc_path)))[0]
                    out_fc_path = os.path.join(gdb_path, fc_name)

                    # Clip inputs by analysis area and dissolve on field
                    arcpy.Clip_analysis(input_fc_path, clip_boundary, "in_memory\\fc_clip")
                    arcpy.Dissolve_management("in_memory\\fc_clip", out_fc_path, field)
                    arcpy.AddMessage('[+] Successfully processed [{}]'.format(fc_name))
                except:
                    arcpy.AddMessage('[-] Could not process [{}]'.format(fc_name))

            arcpy.AddMessage(params[2].valueAsText)
            rasters = params[2].valueAsText.split(';')
            arcpy.AddMessage(str(rasters))
            for raster in rasters:
                try:
                    raster_name = os.path.basename(raster)
                    out_raster = os.path.join(gdb_path, 'raster_' + raster_name)

                    arcpy.Clip_management(
                        in_raster=raster,
                        out_raster=out_raster,
                        in_template_dataset=clip_boundary,
                        clipping_geometry="ClippingGeometry",
                        )

                    # raster_name = os.path.basename(raster)
                    # arcpy.AddMessage(raster_name)
                    # out_raster_path = os.path.join(gdb_path, raster_name)
                    # arcpy.AddMessage(out_raster_path)
                    # arcpy.Clip_management(
                    #     in_raster=raster,
                    #     out_raster=out_raster_path,
                    #     in_template_dataset=clip_boundary,
                    #     clipping_geometry="ClippingGeometry",
                    #     )
                    arcpy.AddMessage('[+] Successfully processed [{}]'.format(raster_name))
                except:
                    arcpy.AddMessage('[-] Could not process [{}]'.format(raster_name))
                    arcpy.AddMessage(str(traceback.format_exc()))

        except:  #                             _(shit)
            # [Rails]  <-- -->  [Car]  >--<o _/
            err = str(traceback.format_exc())
            arcpy.AddError(err)

        finally:
            deleteInMemory()