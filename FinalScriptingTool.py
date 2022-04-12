import arcpy
import numpy as np
import os

input_gdb = arcpy.GetParameterAsText(0)
blkgrp = arcpy.GetParameterAsText(1)
cntbnd = arcpy.GetParameterAsText(2)
mjr_rds = arcpy.GetParameterAsText(3)
schools = arcpy.GetParameterAsText(4)
hospitals = arcpy.GetParameterAsText(5)
crashes = arcpy.GetParameterAsText(6)
landuse = arcpy.GetParameterAsText(7)

residential_SQL = arcpy.GetParameterAsText(8)
commercial_SQL = arcpy.GetParameterAsText(9)

raster_ext_mask = arcpy.GetParameterAsText(10)
cell_size = arcpy.GetParameterAsText(11)

res_weight = arcpy.GetParameterAsText(12)
comm_weight = arcpy.GetParameterAsText(13)
school_weight = arcpy.GetParameterAsText(14)
hospital_weight = arcpy.GetParameterAsText(15)

output_fc = arcpy.GetParameterAsText(16)
output_gdb = arcpy.GetParameterAsText(17)

arcpy.env.workspace = input_gdb

SCPT_Alachua_blkgrps = output_gdb +"//SCPT_Alachua_blkgrps"

#Select the Blockgroups that are in Alachua County
Alachua_blkgrps = arcpy.management.SelectLayerByLocation(in_layer=blkgrp, overlap_type= "WITHIN",
                                    select_features=cntbnd)
arcpy.conversion.FeatureClassToFeatureClass(Alachua_blkgrps, output_gdb, "SCPT_Alachua_blkgrps")

SCPT_mjr_rd_blkgrp_intersection = output_gdb + "//SCPT_mjr_rd_blkgrp_intersection"
#Find the Intersection of Major Roads and the Alachua Block Groups
mjr_rd_intersection = arcpy.analysis.Intersect(in_features=[SCPT_Alachua_blkgrps, mjr_rds],
                        out_feature_class = SCPT_mjr_rd_blkgrp_intersection)

#Calculate the Summary Statistics for the Intersected Major Roads
mjr_rd_stats = arcpy.analysis.Statistics(in_table = SCPT_mjr_rd_blkgrp_intersection, 
            out_table = output_gdb + "//"+ "SCPT_mjr_rd_blkgrp_stats",
            statistics_fields = [["Shape_Length", "SUM"]],
            case_field = "GEOID10")

#Buffer the Major Roads by 300 ft
mjr_rd_buffer300ft = arcpy.Buffer_analysis(in_features = mjr_rds,
                        buffer_distance_or_field= "300 feet")

#Select Crashes Within the 300 ft buffer   
mjr_rd_crashes = arcpy.management.SelectLayerByLocation(in_layer = crashes,
                                       overlap_type= "WITHIN",
                                       select_features= mjr_rd_buffer300ft)
arcpy.conversion.FeatureClassToFeatureClass(in_features= mjr_rd_crashes, 
                                                out_path= output_gdb,
                                                out_name=  "SCPT_mjr_rd_crashes")      

SCPT_Alachua_blkgrp_crashes_join = output_gdb + "//"+ "SCPT_Alachua_blkgrp_crashes_join"

#Spatial Join the Crashes to the Alachua Block Group
Alachua_blkgrp_crashes_join = arcpy.analysis.SpatialJoin(target_features = Alachua_blkgrps, join_features = mjr_rd_crashes, 
                    out_feature_class = SCPT_Alachua_blkgrp_crashes_join, match_option = 'CONTAINS')       

#Join the Summary Table to the Alachua Block Group    
arcpy.management.JoinField(in_data = SCPT_Alachua_blkgrp_crashes_join, in_field = "GEOID10", 
                               join_table = output_gdb + "//"+ "SCPT_mjr_rd_blkgrp_stats", join_field = "GEOID10")

#Calculate the Crashes per linear mile
arcpy.AddField_management(in_table=SCPT_Alachua_blkgrp_crashes_join, field_name="Linear_mi", field_type = "DOUBLE")
arcpy.CalculateField_management(in_table=SCPT_Alachua_blkgrp_crashes_join, field="Linear_mi", 
                                    expression = "!SUM_Shape_Length!/1609", expression_type = "PYTHON3")

arcpy.AddField_management(in_table=SCPT_Alachua_blkgrp_crashes_join, field_name="Crash_rate", field_type = "DOUBLE")
arcpy.CalculateField_management(in_table=SCPT_Alachua_blkgrp_crashes_join, field="Crash_rate", 
                                expression = "!Join_Count!/!Linear_mi!", expression_type = "PYTHON3")

##############Classification of the Class Rate############################
input_arr = arcpy.da.TableToNumPyArray(in_table=SCPT_Alachua_blkgrp_crashes_join, field_names="Crash_rate", null_value = 0)

min_val = np.min(input_arr['Crash_rate'])
max_val = np.max(input_arr['Crash_rate'])
mean_val = np.mean(input_arr['Crash_rate'])

arcpy.management.AddField(in_table=SCPT_Alachua_blkgrp_crashes_join, field_name = "Crash_class", field_type ="TEXT")

codeblock = """
def classify(x, meanval, minval, maxval):
    if x is None:
        return "Low"
    elif x < (meanval + minval)/2:
        return "Low"
    elif (meanval + minval)/2 <= x < (meanval + maxval)/2:
        return "Medium"
    else:
        return "High"
"""

expression = "classify(!{}!, {}, {}, {})".format("Crash_rate", mean_val,
                            min_val, max_val)
arcpy.management.CalculateField(in_table=SCPT_Alachua_blkgrp_crashes_join, field = "Crash_class", expression = expression, 
                expression_type ="PYTHON3", code_block = codeblock)

#Part II Suitability Analysis
#Calculate the Euclidean Distance for each of the specified land uses
arcpy.env.mask = cntbnd
arcpy.env.extent = cntbnd
arcpy.env.snapRaster = raster_ext_mask
arcpy.env.cellSize = 30

#Schools
school_dist_rast = arcpy.sa.EucDistance(in_source_data=schools, cell_size=30)   # cell_size using default setup in env parameter
school_dist_rast.save(os.path.join(output_gdb, 'SCPT_EucDist_schools'))
# arcpy.conversion.FeatureClassToFeatureClass(in_features= school_dist_rast, 
#                                             out_path= output_gdb,
#                                             out_name= "SCPT_EucDist_schools")

arcpy.env.snapRaster = school_dist_rast

school_dist_1_9 = arcpy.sa.Slice(
    school_dist_rast,
    number_zones=9, 
    slice_type="NATURAL_BREAKS",
    base_output_zone=1
)

school_dist_1_9.save(output_gdb + "//" + "SCPT_school_dist_1_9")
# arcpy.conversion.FeatureClassToFeatureClass(in_features= school_dist_1_9, 
#                                             out_path= output_gdb,
#                                             out_name= "school_dist_1_9")

#Residential Landuse
residentialLU = arcpy.management.SelectLayerByAttribute(in_layer_or_view = landuse, where_clause= residential_SQL)

arcpy.conversion.FeatureClassToFeatureClass(in_features= residentialLU, 
                                                    out_path= output_gdb,
                                                    out_name= "SCPT_residential_LU")

SCPT_residential_LU = os.path.join(output_gdb, "SCPT_residential_LU")
residentialLU_rast = arcpy.sa.EucDistance(in_source_data = SCPT_residential_LU, cell_size=cell_size)   # cell_size using default setup in env parameter
residentialLU_rast.save(os.path.join(output_gdb, 'SCPT_EucDist_res_LU'))

# arcpy.conversion.FeatureClassToFeatureClass(in_features= residentialLU_rast, 
#                                             out_path= output_gdb,
#                                             out_name= "SCPT_residential_LU")

arcpy.env.snapRaster = residentialLU_rast

residentialLU_dist_1_9 = arcpy.sa.Slice(
    residentialLU_rast,
    number_zones=9, 
    slice_type="NATURAL_BREAKS",
    base_output_zone=1
)

residentialLU_dist_1_9.save(output_gdb + "//" + "SCPT_residentialLU_dist_1_9")
# arcpy.conversion.FeatureClassToFeatureClass(in_features= residentialLU_dist_1_9, 
#                                             out_path= output_gdb,
#                                             out_name= "SCPT_residentialLU_dist_1_9")

#Commercial Landuse
commercialLU = arcpy.management.SelectLayerByAttribute(in_layer_or_view = landuse, where_clause= commercial_SQL)

arcpy.conversion.FeatureClassToFeatureClass(in_features= commercialLU, 
                                                out_path= output_gdb,
                                                out_name= "SCPT_commercial_LU")

# Py_commercial_LU = os.path.join(output_gdb, "Py_commercial_LU")
commercialLU_rast = arcpy.sa.EucDistance(in_source_data= commercialLU, cell_size=cell_size)   # cell_size using default setup in env parameter
commercialLU_rast.save(os.path.join(output_gdb, 'SCPT_EucDist_comm_LU'))

# arcpy.conversion.FeatureClassToFeatureClass(in_features= commercialLU_rast, 
#                                             out_path= output_gdb,
#                                             out_name= "SCPT_EucDist_comm_LU")

arcpy.env.snapRaster = commercialLU_rast

commercialLU_dist_1_9 = arcpy.sa.Slice(
    commercialLU_rast,
    number_zones=9, 
    slice_type="NATURAL_BREAKS",
    base_output_zone=1
)             


#Hospitals
hospitals_dist_rast = arcpy.sa.EucDistance(in_source_data= hospitals, cell_size=cell_size)   # cell_size using default setup in env parameter
hospitals_dist_rast.save(os.path.join(output_gdb, 'SCPT_EucDist_hospitals'))
# arcpy.conversion.FeatureClassToFeatureClass(in_features= hospitals_dist_rast, 
#                                                 out_path= output_gdb,
#                                                 out_name= "SCPT_EucDist_hospitals")

arcpy.env.snapRaster = hospitals_dist_rast

hospitals_dist_1_9 = arcpy.sa.Slice(
    hospitals_dist_rast,
    number_zones=9, 
    slice_type="NATURAL_BREAKS",
    base_output_zone=1
)

# arcpy.conversion.FeatureClassToFeatureClass(in_features= hospitals_dist_1_9, 
#                                             out_path= output_gdb,
#                                             out_name= "Py_hospitals_dist_1_9")

dist_remap = arcpy.sa.RemapValue(
    [[1, 9], [2, 8], [3, 7], [4, 6],
     [6, 4], [7, 3], [8, 2], [9, 1]])

hospitals_dist_9_1 = arcpy.sa.Reclassify(hospitals_dist_1_9, "Value", dist_remap)

    # arcpy.conversion.FeatureClassToFeatureClass(in_features= hospitals_dist_9_1, 
    #                                             out_path= output_gdb,
    #                                             out_name= "Py_hospitals_dist_9_1")

#Take the Weighted Sum of the Slices 
# school_rast = output_gdb + "\\" + "Py_school_dist_1_9"
# hospital_rast = output_gdb + "\\" + "Py_hospitals_dist_9_1"
# residential_rast = output_gdb + "\\" + "Py_residentialLU_dist_1_9"
# commercial_rast = output_gdb + "\\" + "Py_commercialLU_dist_1_9"

WSumTable = arcpy.sa.WSTable([[school_dist_1_9, "Value", school_weight], [hospitals_dist_9_1, "Value", hospital_weight],
                        [residentialLU_dist_1_9, "Value", res_weight], [commercialLU_dist_1_9, "Value", comm_weight]])


outWeightedSum = arcpy.sa.WeightedSum(WSumTable)
# outWeightedSum.save(output_gdb+ "\\" + "Py_WeightedSum")

Py_WeightedSum_dist_1_9 = arcpy.sa.Slice(
        in_raster= outWeightedSum,
        number_zones=9, 
        slice_type="NATURAL_BREAKS",
    )

# Py_WeightedSum_dist_1_9.save(output_gdb + "\\" + "Py_WeightedSum_dist_1_9")

#Reclassify the Weights with 9 being the most desirable
Py_WeightedSum_dist_9_1 = arcpy.sa.Reclassify(Py_WeightedSum_dist_1_9, "Value", dist_remap)
# Py_WeightedSum_dist_9_1.save(output_gdb +"//"+ "Py_WeightedSum_reclass")


#Raster to Polygon
arcpy.conversion.RasterToPolygon(in_raster=Py_WeightedSum_dist_9_1,
            out_polygon_features= output_gdb + "//"+ "SCPT_WeightedSum_reclass_poly")


#Selection of all the 9s
Py_WSum_9 = arcpy.management.SelectLayerByAttribute(in_layer_or_view = output_gdb + "//"+ "SCPT_WeightedSum_reclass_poly", 
                                                    where_clause= """"gridcode" = 9""")


#Final selection
#Intersection of the High Crash Rate and the Weighted Sum files
Py_Final_Intersect = arcpy.analysis.Intersect(in_features=[SCPT_Alachua_blkgrp_crashes_join, Py_WSum_9], 
                                                   out_feature_class= output_gdb +"//" +"SCPT_Final_Intersect")           

#Add "Shape_Area" of Intersection so it's easier to distinguish before joining: just add a new name with same value
arcpy.AddField_management(in_table=Py_Final_Intersect, field_name="WSum_Area", field_type = "DOUBLE")
arcpy.CalculateField_management(in_table=Py_Final_Intersect, field="WSum_Area",
                                expression = "!Shape_Area!/4047", expression_type = "PYTHON3")

#Join the intersection to the original blk group layer
arcpy.management.JoinField(in_data = SCPT_Alachua_blkgrp_crashes_join, in_field = "GEOID10", 
                            join_table = Py_Final_Intersect, join_field = "GEOID10", fields="WSum_Area")    

arcpy.AddField_management(in_table=SCPT_Alachua_blkgrp_crashes_join, field_name="Area_ratio", field_type = "DOUBLE")
arcpy.CalculateField_management(in_table= SCPT_Alachua_blkgrp_crashes_join, field="Area_ratio", 
                                expression = "!WSum_Area! / !ACRES!", expression_type = "PYTHON3")       

Py_FinalSelection = arcpy.management.SelectLayerByAttribute(in_layer_or_view = SCPT_Alachua_blkgrp_crashes_join, 
                                                        where_clause= """Crash_class = 'High'""")        

arcpy.conversion.FeatureClassToFeatureClass(in_features= Py_FinalSelection, 
                                                out_path= output_gdb,
                                                out_name= "SCPT_FinalSelection")      

Py_FinalSelection2 = arcpy.management.SelectLayerByAttribute(in_layer_or_view = output_gdb + "//SCPT_FinalSelection", 
                                                        where_clause= '"Area_ratio" > 0.5')

arcpy.conversion.FeatureClassToFeatureClass(in_features= output_gdb + "//SCPT_FinalSelection", 
                                                out_path= output_gdb,
                                                out_name= output_fc)                                                                                  