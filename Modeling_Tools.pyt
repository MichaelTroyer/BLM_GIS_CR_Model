import arcpy
import traceback


def buildWhereClauseFromList(table, field, valueList):
    """Takes a list of values and constructs a SQL WHERE
       clause to select those values within a given field and table."""
    
    # arcpy.AddMessage('table: ' + table)
    arcpy.AddMessage('field: ' + field)
    arcpy.AddMessage('vaues: ' + repr(valueList))
    
    # Add DBMS-specific field delimiters
    fieldDelimited = arcpy.AddFieldDelimiters(arcpy.Describe(table).path, field)

    # Determine field type
    fieldType = arcpy.ListFields(table, field)[0].type

    # Add single-quotes for string field values
    if str(fieldType) == 'String':
        valueList = ["'%s'" % value for value in valueList]
    
    # Format WHERE clause in the form of an IN statement
    whereClause = "%s IN (%s)" % (fieldDelimited, ', '.join(map(str, valueList)))

    arcpy.AddMessage('where = ' + whereClause)
    
    return whereClause


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Modelling Tools"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Collapse_Categories_into_Largest_Neighbors]


class Collapse_Categories_into_Largest_Neighbors(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Collapse_Categories_into_Largest_Neighbors"
        self.description = ""
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param00=arcpy.Parameter(
            displayName="Input Categorical Feature Class",
            name="Input_FC",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            )
        
        param01=arcpy.Parameter(
            displayName="Select Analysis Field",
            name="Select_Field",
            datatype="String",
            parameterType="Required",
            direction="Input",
            )

        param02=arcpy.Parameter(
            displayName="Select Values to Collapse",
            name="Select_Values",
            datatype="String",
            parameterType="Required",
            multiValue="True",
            direction="Input",
            enabled = "False",
            )

        param03=arcpy.Parameter(
            displayName="Output Feature Class",
            name="Output_FC",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output",
            )

        params = [param00, param01, param02, param03]
        return params


    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, params):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        
        in_featrs = params[0]
        sel_field = params[1]
        sel_value = params[2]

        sel_field.filter.type = "ValueList"
        sel_field.filter.list = [f.name for f in arcpy.Describe(in_featrs.valueAsText).fields
                                 if f.type in ["String", "Integer"]]

        if sel_field.value:
            sel_value.enabled = "True"
            sel_value.filter.type = "ValueList"
            arcpy.Frequency_analysis(in_featrs.valueAsText, "in_memory\\freq", sel_field.value)
            featrs = list(set([row[0] for row
                        in arcpy.da.SearchCursor("in_memory\\freq", sel_field.value) if row[0]]))
            featrs.sort()
            sel_value.filter.list = featrs

        else:
            sel_value.enabled = "False"

        return

    def updateMessages(self, params):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, params, messages):
        """The source code of the tool."""
        try:
            in_features  = params[0].valueAsText
            case_field   = params[1].valueAsText
            case_values  = [s.strip("'") for s in params[2].valueAsText.split(';')]
            out_features = params[3].valueAsText

            arcpy.MultipartToSinglepart_management(in_features, 'in_memory\\singlepart')
            arcpy.MakeFeatureLayer_management('in_memory\\singlepart', 'in_memory\\layer')
            where = buildWhereClauseFromList('in_memory\\layer', case_field, case_values)
            arcpy.SelectLayerByAttribute_management("in_memory\\layer", 'NEW_SELECTION', where)
            
            # Eliminate features by largest area
            arcpy.Eliminate_management("in_memory\\layer", out_features, "AREA")

        except:
            # [Rails]  <-- -->  [Car]
            err = str(traceback.format_exc())
            arcpy.AddError(err)