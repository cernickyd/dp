# -*- coding: utf-8 -*-

import arcpy
from math import sqrt, acos, pi

def editLines (fc_bacis,fc_output,dem,inExField):
    #   Edits lines and prepares them for next step
    fc_extra = fc_output + "/extra"
    fc_intra = fc_output + "/intra"
    fc_merge = fc_output + "/merge"
    fc_dissolve = fc_output + "/dissolve"
    fc_SpatialJoin = fc_output + "/SpatialJoin"
    fc_input = fc_output + "/ready_lines"
    expression1 = inExField + "=1"
    expression2 = inExField + "=0"
    arcpy.FeatureClassToFeatureClass_conversion(fc_bacis, fc_output, "intra", expression1)
    arcpy.FeatureClassToFeatureClass_conversion(fc_bacis, fc_output, "extra", expression2)
    arcpy.Dissolve_management(fc_extra, fc_dissolve, "", "", "SINGLE_PART", "UNSPLIT_LINES")
    arcpy.SpatialJoin_analysis(fc_dissolve, fc_bacis, fc_SpatialJoin, "JOIN_ONE_TO_ONE", "", "", "INTERSECT")
    arcpy.Merge_management([fc_SpatialJoin, fc_intra], fc_merge)
    arcpy.InterpolateShape_3d(dem, fc_merge, fc_input)

    return fc_input

def splitPoints( exp, fc_input, fc_output, speedField, pnts_fc):
    splitPnts = []
    with arcpy.da.SearchCursor(fc_input, [speedField, "SHAPE@"], exp) as cursor:
        #   Reads the shape of feature class
        for row in cursor:
            polyLine = row[1]
            maxSpd = True
            #   Opens every line in feature class
            for line in polyLine:
                #print("next")
                total = len(line)
                start = 0
                end = 0
                #   Selects lines with less than 3 vertexes
                if total < 3:
                    v = row[0]
                    delta = abs(v - row[0])
                    if delta > 10 and maxSpd == True:
                        splitPnts.append(pnt1)
                        maxSpd = False

                    if v == row[0] and maxSpd == False:
                        splitPnts.append(pnt1)
                        maxSpd = True

                    #print(row[0], v)
                #   Selects lines with more than 3 vertexes
                else:
                    #   Calculates curvature of line. Takes 3 vertexes and calculates angle between them.
                    for i in range(0, total - 2):
                        pnt1 = line[i]
                        pnt2 = line[i + 1]
                        pnt3 = line[i + 2]
                        a = sqrt((pnt2.X - pnt1.X) ** 2 + (pnt2.Y - pnt1.Y) ** 2)
                        b = sqrt((pnt3.X - pnt2.X) ** 2 + (pnt3.Y - pnt2.Y) ** 2)
                        c = sqrt((pnt1.X - pnt3.X) ** 2 + (pnt1.Y - pnt3.Y) ** 2)
                        kosin = ((a * a + b * b - c * c) / (2 * a * b))
                        if kosin < -1:
                            kosin = -1
                        #   Gamma = angle between 3 vertexes
                        gamma = acos(kosin) * 180 / pi
                        odchylka = 180 - gamma
                        l = (a + b) / 1000
                        #   Curvature
                        K = odchylka / l
                        #   If all 3 vertexes lays on one straight line
                        if K == 0:
                            o = 100000
                        #   If not, calculates speed limit
                        else:
                            o = (360 / K) * l
                        r = o / 2 * pi
                        v0 = 3.6 * sqrt(9.81 * (r * 1000) * (0.25 + 0.01 * 2.5))
                        #   If the speed in the curve is higher than design speed, sets design speed.
                        if v0 > row[0]:
                            v = row[0]
                        #   Adds influence of slope
                        else:
                            v = v0

                        delta = abs(v - row[0])
                        if delta > 10 and maxSpd == True and i > (end + 2):
                            splitPnts.append(pnt2)
                            maxSpd = False
                            start = i

                        if v == row[0] and maxSpd == False and i > (start + 2):
                            splitPnts.append(pnt2)
                            maxSpd = True
                            end = i

                        #print(row[0], v)
    print(len(splitPnts))

    #   Creating split points

    arcpy.CreateFeatureclass_management(fc_output, 'splitPnts', 'POINT', spatial_reference=fc_input)
    with arcpy.da.InsertCursor(pnts_fc, ["SHAPE@XY"]) as cursor:
        for pnt in splitPnts:
            xy = (pnt.X, pnt.Y)
            cursor.insertRow([xy])

def addAttributes(finalLines):
    arcpy.AddField_management(finalLines, "ftSpeed", "FLOAT")
    arcpy.AddField_management(finalLines, "tfSpeed", "FLOAT")
    arcpy.AddField_management(finalLines, "ftTime", "FLOAT")
    arcpy.AddField_management(finalLines, "tfTime", "FLOAT")
    arcpy.AddField_management(finalLines, "avgTime", "FLOAT")
    
def addSpeeds(koef, finalLines):
    with arcpy.da.UpdateCursor(finalLines, ["ftSpeed", "tfSpeed", "rychlost", "intravilan", "SHAPE@"]) as cursor:
        #   Reads the shape of feature class
        for row in cursor:
            if row[3] == 1:
                row[0] = row[2]
                row[1] = row[2]
                cursor.updateRow(row)
                continue

            polyLine = row[4]
            # vs = row[0] # avarage speed of previous segment
            vArrayFT = []
            vArrayTF = []
            lenArray = []
            #   Opens every line in feature class
            for line in polyLine:
                print("next")
                total = len(line)

                #   Selects lines with less than 3 vertexes
                if total < 3:
                    pnt1 = line[0]
                    pnt2 = line[1]
                    v = row[2]
                    a = sqrt((pnt2.X - pnt1.X) ** 2 + (pnt2.Y - pnt1.Y) ** 2)

                    deltaZ = pnt2.Z - pnt1.Z
                    #   Calculates slope
                    slopeFT = (deltaZ / a) * 100
                    slopeTF = (-deltaZ / a) * 100
                    if slopeFT > 3:
                        vft = v * (1 - koef * slopeFT)
                    else:
                        vft = v

                    if slopeTF > 3:
                        vtf = v * (1 - koef * slopeTF)
                    else:
                        vtf = v

                    row[0] = vft
                    row[1] = vtf
                    cursor.updateRow(row)

                    print(row[2], vft, vtf)
                #   Selects lines with more than 3 vertexes
                else:
                    #   Calculates curvature of line. Takes 3 vertexes and calculates angle between them.
                    for i in range(0, total - 2):
                        pnt1 = line[i]
                        pnt2 = line[i + 1]
                        pnt3 = line[i + 2]
                        a = sqrt((pnt2.X - pnt1.X) ** 2 + (pnt2.Y - pnt1.Y) ** 2)
                        b = sqrt((pnt3.X - pnt2.X) ** 2 + (pnt3.Y - pnt2.Y) ** 2)
                        c = sqrt((pnt1.X - pnt3.X) ** 2 + (pnt1.Y - pnt3.Y) ** 2)
                        kosin = ((a * a + b * b - c * c) / (2 * a * b))
                        if kosin < -1:
                            kosin = -1
                        #   Gamma = angle between 3 vertexes
                        gamma = acos(kosin) * 180 / pi
                        odchylka = 180 - gamma
                        l = (a + b) / 1000
                        #   Curvature
                        K = odchylka / l
                        #   If all 3 vertexes lays on one straight line
                        if K == 0:
                            o = 100000
                        #   If not, calculates speed limit
                        else:
                            o = (360 / K) * l
                        r = o / 2 * pi
                        v0 = 3.6 * sqrt(9.81 * (r * 1000) * (0.25 + 0.01 * 2.5))
                        #   If the speed in the curve is higher than design speed, sets design speed.
                        if v0 > row[2]:
                            v = row[2]
                        else:
                            v = v0

                        #   Calculates slope
                        l = l * 1000
                        deltaZ = pnt3.Z - pnt1.Z
                        slopeFT = (deltaZ / l) * 100
                        slopeTF = -(deltaZ / l) * 100

                        #   Adds influence of slope
                        if slopeFT > 3:
                            vft = v * (1 - koef * slopeFT)
                        else:
                            vft = v

                        if slopeTF > 3:
                            vtf = v * (1 - koef * slopeTF)
                        else:
                            vtf = v

                        lenArray.append(l)
                        vArrayTF.append(vtf)
                        vArrayFT.append(vft)

                    vMeanFT = meanWght(lenArray, vArrayFT)
                    vMeanTF = meanWght(lenArray, vArrayTF)

                    row[0] = vMeanFT
                    row[1] = vMeanTF
                    cursor.updateRow(row)

                vArrayTF.clear()
                vArrayFT.clear()
                lenArray.clear()

def meanWght(lenArray, vArray):
    wghts = []
    mean = 0
    suma = sum(lenArray)
    for i in range(0,len(lenArray)-1):
        wghts.append(lenArray[i]/suma)
        mean = mean + vArray[i]*wghts[i]
    return mean

# -*- coding: utf-8 -*-
#   Imports
import arcpy, functions

#   To allow overwriting the outputs change the overwrite option to true.
arcpy.env.overwriteOutput = True

#   Inputs
fc_bacis = 'C:/Users/cerni/Desktop/Univerzita Karlova/Geografie/Diplomka/data/silnicni_sit.gdb/silnicni_sit_test_multiple'
fc_output = 'C:/Users/cerni/Desktop/Univerzita Karlova/Geografie/Diplomka/data/silnicni_sit.gdb'
dem = 'C:/Users/cerni/Desktop/Univerzita Karlova/Geografie/Diplomka/data/silnicni_sit.gdb/dem'
speedField = 'rychlost'
inExField = 'intravilan'
pnts_fc = fc_output + '/splitPnts'

#   Defining variables
exp = inExField + "= 0"
koef = 0.01813

# fc_bacis = arcpy.GetParameterAsText(0)
# fc_output = arcpy.GetParameterAsText(1)
# dem = arcpy.GetParameterAsText(2)
# speedField = arcpy.GetParameterAsText(3)
# inExField = arcpy.GetParameterAsText(4)

#   Cursor for reading rows in line input
fc_input = functions.editLines(fc_bacis,fc_output,dem)
functions.splitPoints(exp, fc_input, fc_output, speedField, pnts_fc)

finalLines = fc_output + '/finalLines'
arcpy.SplitLineAtPoint_management(fc_input,pnts_fc,finalLines,"0.25 Meters")
arcpy.AddField_management(finalLines, "ftSpeed", "FLOAT")
arcpy.AddField_management(finalLines, "tfSpeed", "FLOAT")

functions.addSpeeds(koef,finalLines)
