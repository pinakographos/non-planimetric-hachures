# This script runs in QGIS, and uses QGIS things. But it's not really
# doing anything that's inherently spatial. Just analyzing images and
# drawing lines.

#============================USER PARAMETERS============================
# These two parameters below are in image pixel units. So choosing 6 for
# the max_hachure_spacing means the script aims to make hachures 6 px
# apart when the shading raster is at its lightest
min_hachure_spacing = 3
max_hachure_spacing = 9

# When done, filter out all hachure lines shorter than this (in pixels)
stubFilter = 1 * max_hachure_spacing

# Enter the name of the raster layers we'll be using.  We need one to
# control the density of the hachure shading, and another with upslope
# data to control the direction of the lines. These are Blender images.

shading_layer_name = 'Shading'
upslope_layer_name = 'Upslope'

# Optionally, you can include a set of existing lines that show major
# breaks in the terrain (e.g., a line defining the top of a ridge or
# hill). The script will avoid running hachures over these lines.

# Set to None if you do not plan to use this.

breaks_layer_name = None

#============================PREPATORY WORK=============================
#--------STEP 0: Import various modules and such that are needed--------

import random
import time
import math
import os
import statistics

from collections import defaultdict

from qgis.PyQt.QtCore import (
    QVariant
)
from qgis.utils import iface
from qgis.core import (
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsField,
    QgsMemoryProviderUtils,
    QgsProcessingFeatureSourceDefinition,
    QgsPointXY,
    QgsGeometry,
    QgsFeature,
    QgsWkbTypes,
    edit
)
from qgis import processing

# This function reports back any errors found later on

def warn_user(error_type):
    # Here are our various error messages and levels
    # The format is ErrorNumber: (Text,Level)
    
    error_dict = {
        0: ('Done! Enjoy your freshly baked hachures!',
            Qgis.Success),
        1: ('One or both raster layer names incorrect',
            Qgis.Critical),
        2: ('min_hachure_spacing must not be more than '
            'max_hachure_spacing.',
            Qgis.Critical),
        3: ('min_hachure_spacing must be greater than 0',
            Qgis.Critical),
        4: ('max_hachure_spacing must be greater than 0',
            Qgis.Critical),
        5: ('Breaks layer not found; set to "None" if not being used',
            Qgis.Critical)
    }
    
    err = error_dict[error_type]
    
    iface.messageBar().pushMessage('Hachure Script',*err)
    
    if err[1] == Qgis.Critical:
        raise Exception(err[0])


#-----------------STEP 1: Set up some  useful variables-----------------
spacing_range = max_hachure_spacing - min_hachure_spacing
instance = QgsProject.instance()
crs = instance.crs()

phase = 0
# The script operates in two phases, starting with Phase 0.
# In Phase 0, we draw mostly-vertical hachures.
# In Phase 1, we draw mostly-horizontal hachures.


# If there's a breaks layer, it's treated as an existing set of lines to
# retain and work around.

if breaks_layer_name:
    try:
        start_hachures = list(
            instance.mapLayersByName(breaks_layer_name)[0].getFeatures()
            )
        breaks = start_hachures.copy()
        current_hachures = start_hachures.copy()
    except:
        warn_user(5)
else:
    start_hachures = None
    current_hachures = None


#------------------STEP 2: Handling Basic Input Errors------------------

checks = [
    (spacing_range < 0,2),
    (min_hachure_spacing <= 0,3),
    (max_hachure_spacing <= 0,4)
]

for condition,code in checks:
    if condition:
        warn_user(code)
        break

#------------------------STEP 3: Prep raster layers---------------------

try:
    upslope = instance.mapLayersByName(upslope_layer_name)[0]
    shading = instance.mapLayersByName(shading_layer_name)[0]
except:
    warn_user(1)
    
if (
    upslope.type() != QgsMapLayer.RasterLayer or
    shading.type() != QgsMapLayer.RasterLayer
):
    warn_user(1)

# We have the rasters we need; now set them up to be read later on.

provider = upslope.dataProvider()
extent = provider.extent()
rows = upslope.height()
cols = upslope.width()

# Our upslope image is 3 bands, giving the x/y/z components; we ignore
# the z component, though, as it's not needed.

upslope_blocks = [provider.block(x,extent,cols,rows) for x in [1,2]]

shading_block = shading.dataProvider().block(1, extent, cols, rows)

cell_width = extent.width() / cols
cell_height = extent.height() / rows

average_pixel_size = 0.5 * (upslope.rasterUnitsPerPixelX() +
                  upslope.rasterUnitsPerPixelY())

# This parameter can be tweaked if you find your lines are a bit jagged.

jump_distance = average_pixel_size * 3

#===========================CLASS DEFINITIONS===========================
#-------Check lines are used to check the spacing of the hachures-------

class Check_Line:
    def __init__(self,geometry):
        self.geometry = geometry
        self.length = geometry.length()
        
    def split_by_hachures(self):
        #Split this check_line according to our current list of hachures
        all_segments = []

        intersection_points = []
        for hachure_feature in current_hachures:
            hachure_geometry = hachure_feature.geometry()
            point = self.geometry.intersection(hachure_geometry)
            if point.wkbType() == QgsWkbTypes.MultiPoint:
                intersection_points += [CutPoint(
                    QgsGeometry.fromPointXY(p),hachure_feature)
                    for p in point.asMultiPoint()]
            elif point.wkbType() == QgsWkbTypes.Point:
                intersection_points += [CutPoint(point, hachure_feature)]
            # The intersection can return Empty or (rarely) 
            # a geometryCollection. We can safely skip over these
        
        for point in intersection_points:
            # This tells us where along the line to cut
            point.cut_location = self.geometry.lineLocatePoint(
                                     point.geometry)
                
        if len(intersection_points) > 0:
            # If we found intersections, use them to cut the line
            check_line_segments = cutpoint_splitter(self.geometry,
                                            intersection_points)
            all_segments += check_line_segments
        else:
            # If not, we should still return the unbroken line
            line_feature = QgsFeature()
            line_feature.setGeometry(self.geometry)
            all_segments.append(Segment(line_feature))
            
        return all_segments

#---Segments are check line pieces used to space or generate hachures---
class Segment:
    def __init__(self,segFeature):
        self.geometry = segFeature.geometry()
        self.length = self.geometry.length()
        self.feature = segFeature
        self.hachures = []

        # We'll want this set of vertices covering all along this line
        # for use in a couple of functions later
        
        densified_line = self.geometry.densifyByDistance(average_pixel_size)
        self.vertices = [(vertex.x(), vertex.y())
                    for vertex in densified_line.vertices()]

        # And here they are: these two functions use the densified set
        self.avg_xComp,self.avg_yComp = self.avg_xyComp()
        self.spacing = self.ideal_spacing()

        #:::::SET STATUS:::::
        self.status = None
        
        # Status stores info on how this segment should affect hachures
        # These values are used later in the subsequent_line function
        
        # The s factor (below) improves spacing as lines tilt away from
        # our check axis (horizontal or vertical). It forces highly-
        # angled lines to be spaced farther apart, by pretending they're
        # closer together than they really are. The s factor will always
        # be in the range from 0.7 to 1.
        # It's based on the value out of avg_xyComp.

        if phase == 0:
            component = self.avg_yComp
        else:
            component = self.avg_xComp
        
        factorMin = 0.7
 
        s = (1 - factorMin) * component + factorMin

        
        if self.length * s < (self.spacing * 0.9):
            self.status = 1
        elif self.length * s > (self.spacing * 2.2):
            self.status = 2

        # The 0.9 and 2.2 above are thermostat controls. Instead of a
        # line being "too short" when it exactly falls below its ideal
        # spacing, we let it get a little tighter to avoid near-parallel
        # hachures cycling on/off rapidly.
            
    def avg_xyComp(self):
        # What's the average x&y component of the hachure angle in this
        # area (under our Segment)? Helps us improve spacing as lines
        # tilt away from our check axis (horizontal or vertical).
        
        yComponents = []
        xComponents = []

        # Could likely do this more efficiently, but we'll just pretend
        # we're going to calculate a jump from each of these points.
        # That function gives us the hachure angle in the area.
        
        for coords in self.vertices:
            jump = jump_calculator(coords)
    
            if jump == -1:
                continue
                
            yComp = jump[1] / jump_distance
            xComp = jump[0] / jump_distance
            yComponents.append(yComp)
            xComponents.append(xComp)
        
        output = [0,0]
        
        if yComponents != []:
            output[1] = abs(statistics.fmean(yComponents))
        if xComponents != []:
            output[0] = abs(statistics.fmean(xComponents))
            
        return output    
        
    def ideal_spacing(self):

        # Examine the shading value under this line segment, and use it
        # to determine the ideal spacing of hachures.
        # A low value = a dark shading raster = denser hachures
        
        row_col_coords = [xy_to_rc(c) for c in self.vertices]
        
        samples = [sample_raster(c,1) for c in row_col_coords]

        pct = statistics.fmean(samples) / 65535

        # The 65535 is because we assume a 16-bit shading image. We are
        # essentially rescaling the shading value to a 0–1 range.

        spacing = min_hachure_spacing + (spacing_range * pct)
    
        return spacing
        
#-------------CutPoints mark where a check line is to be cut------------
class CutPoint:
    def __init__(self,point_geometry,hachure_feature):
        self.geometry = point_geometry
        self.hachure = hachure_feature
        self.cut_location = None

#=========================FUNCTION DEFINITIONS-=========================        
#--------Converts x/y coords to row/col for sampling the rasters--------
def xy_to_rc(location):
    x,y = location
        
    col = round((x - extent.xMinimum()) / cell_width - 0.5)
    row = round((extent.yMaximum() - y) / cell_height - 0.5)
    
    return (row,col)

#-------Samples the uphill (type = 0) or shading (type = 1) raster------

def sample_raster(location, type = 0):
    row,col = location
    
    if row >= rows or col >= cols or row < 0 or col < 0:
        # i.e., if we're out of bounds
        return 0
    
    if type == 0:
        return [block.value(row,col) / (65535.0 / 2) -1
                for block in upslope_blocks]
        # returns a list containing the x and y values from the uphill
        # raster, rescaled to -1.0 to 1.0

    else:
        return shading_block.value(row,col)


#--Take Segments & turn them into dashed lines based on ideal spacing---
def dash_maker(check_line_segment_list):
    
    output_segments = []
    
    for check_line_segment in check_line_segment_list:
        
        spacing = check_line_segment.spacing

        #We tune the spacing value based on the segment length to ensure
        #an integer number of dashes. This is rather like the automatic
        #dash/gap spacing in Adobe Illustrator

        #Our goal here is to split a segment into dashes & gaps, thusly:
        #  ----    ----    ----    ----    ----    ----    ----
        #Each dash length = spacing, surrounded by gaps half that width
        #Thus one unit looks like this: |  ----  |


        total_length = spacing * 2 #the length of a gap + dash + gap
        total_units = round(check_line_segment.length / total_length)
        
        if total_units == 0:
            #Just in case we round down to the point of having 0 dashes
            continue
        
        dash_gap_length = check_line_segment.length / total_units

        dash_width = dash_gap_length / 2
        #half of our gap-dash-gap is the dash

        gap_width = dash_width / 2
        start_point = gap_width
        end_point = dash_width + gap_width

        geometry = check_line_segment.geometry

        while True:
            substring_feature = QgsFeature()
            line_substring = geometry.constGet().curveSubstring(
                start_point, end_point)
            substring_feature.setGeometry(line_substring)

            output_segments.append(Segment(substring_feature))

            start_point += dash_gap_length
            end_point += dash_gap_length

            if end_point > check_line_segment.length:
               break

    if len(output_segments) > 0:       
        return output_segments
        
    else:
        return None 


#-------Turns list of tuples of xy coodinates into a line feature-----
def make_lines(coord_list):
    points = [QgsPointXY(x, y) for x, y in coord_list]
    polyline = QgsGeometry.fromPolylineXY(points)
    feature = QgsFeature()
    feature.setGeometry(polyline)
    
    return feature

#---------------------Cartesian distance calculator---------------------    
def dist(one,two):
    x1,y1 = one
    x2,y2 = two
    
    return math.sqrt((x1-x2)**2 + (y1-y2)**2)


#---Checks along a line to see where hachures need to be trimmed/begun--
def subsequent_line(check_line):
    global current_hachures

    # First we split the line according to the existing hachures
    
    split_line = check_line.split_by_hachures()
    
    too_short = []
    too_long = []

    for segment in split_line:
    
        if segment.status == 1:
            too_short.append(segment)
        elif segment.status == 2:
            too_long.append(segment)

    # too_short: this segment spans 2 hachures that are too close
    # too_long: segment's 2 hachures are too far apart

    # We first find which hachures must be clipped off
    
    to_clip = []

    for seg in too_short:
        hachures = seg.hachures
        if len(hachures) == 2:
            # Some segments won't touch enough hachures,
            # so we need to check there are 2
            
            if start_hachures:
                thinned = [h for h in hachures if h not in start_hachures]
            else:
                thinned = hachures
            # This confines our activity to only the hachures we're
            # actively generating right now. If, for example, we also
            # brought in a breaks layer, we want to keep those untouched

            # Later when we're in Phase 1, it will also prevent us from
            # trimming any Phase 0 hachures.


            # Below: if we had 2 hachures before, and were left with 1,
            # that means the other was one we *have* to keep as it was a
            # start_hachure. So, we need to clip the other.

            if len(thinned) == 1:
                to_clip.extend(thinned)
            elif len(thinned) == 0:
                pass
            else: #clip one randomly
                random.shuffle(hachures)
                to_clip.append(hachures[0])
                
    # to_clip can have duplicates. A hachure may have too_short segments
    # on each side, and both of them choose that particular hachure as
    # the one that needs to be clipped off. So we remove duplicates:
    
    to_clip = list(set(to_clip))
    
    # Remove those to be clipped from the current hachures
    current_hachures = [f for f in current_hachures if f not in to_clip]

    # Clip them, then put them back
    clipped_hachures = haircut(check_line,to_clip)
    current_hachures += clipped_hachures
    
    #Let's next deal with adding new hachures to the too_long segments
    
    made_additions = False
    if len(too_long) > 0:
        
        dashes = dash_maker(too_long)
  
        if dashes: #this could come back with None so we must check
            made_additions = True
            additions = hachure_generator(dashes)
    
    if made_additions:
        current_hachures += additions

#---Clips off hachures that need to stop at this particular check_line--
def haircut(check_line,hachure_list):
    if phase == 0:
        value = check_line.geometry.asPolyline()[0].y()
    else:
        value = check_line.geometry.asPolyline()[0].x()


    # Lines are clipped by reading out their points and keeping only
    # those that are on the correct side of our check line. Check lines
    # can be horizontal or vertical, so this also checks the phase.
    
    clipped = []
    for hachure in hachure_list:
        hachure_geo = hachure.geometry()

        points = hachure_geo.asPolyline()
        
        new_points = []
        
        for pt in points:
            if phase == 0:
                checkVal = pt.y()
            else:
                checkVal = pt.x()
                
            if checkVal <= value:
                new_points.append(pt)
            else:
                break

        new_geom = QgsGeometry.fromPolylineXY(new_points)
        
        feat = QgsFeature()
        feat.setGeometry(new_geom)
        clipped.append(feat)
  
    return clipped

#-------------------Starts our first set of hachures--------------------
def first_line(check_line):
    global current_hachures
    
    # First we split our check_line into even segments
    spacing = max_hachure_spacing * 3 
    output_segments = []

    length = check_line.length
    start_point = 0
    end_point = spacing

    i = spacing
    cut_locations = []
    while i < length:
        cut_locations.append(i)
        i += spacing
            
    output_segments.extend(master_splitter(check_line,cut_locations))

    # Then turn them into dashes
    dashes = dash_maker(output_segments)
    
    if dashes:
        current_hachures = hachure_generator(dashes)

#---Takes a single line geometry and splits it at a list of locations---        
def master_splitter(check_line,cut_locations):
    start_point = 0
    cut_locations.append(check_line.length)
    cut_locations.sort()
    
    segment_list = []
    
    for cut_spot in cut_locations:
        
        line_substring = check_line.geometry.constGet().curveSubstring(
                             start_point,cut_spot)
        new_feature = QgsFeature()
        new_feature.setGeometry(line_substring)
        segment_list.append(Segment(new_feature))
        start_point = cut_spot
        
    return segment_list

#---Like master_splitter, but uses CutPoints instead of cut locations---
def cutpoint_splitter(line_geometry,CutPoint_list):
    CutPoint_list.sort(key = lambda x: x.cut_location)
    
    # CutPoints hold info on what hachure generated them; we want to add
    # that info to the subsequent segments
    
    segment_list = []
    
    # Add first segment
    line_substring = line_geometry.constGet().curveSubstring(
                         0,CutPoint_list[0].cut_location)
    new_feature = QgsFeature()
    new_feature.setGeometry(line_substring)
    segment_list.append(Segment(new_feature))

    # Then do all the middle cuts & append hachure data to the Segments
    for i in range(0,len(CutPoint_list)):
        start_point = CutPoint_list[i]
        start_location = start_point.cut_location
        if i == len(CutPoint_list) - 1:
            # Checks if we're at end of the list & handles final segment
            end_location = line_geometry.length()
        else:
            end_point = CutPoint_list[i+1]
            end_location = end_point.cut_location
        line_substring = line_geometry.constGet().curveSubstring(
                             start_location,end_location)
        new_feature = QgsFeature()
        new_feature.setGeometry(line_substring)
        new_segment = Segment(new_feature)
        segment_list.append(new_segment)
        if i != len(CutPoint_list) - 1:
            new_segment.hachures = [start_point.hachure,end_point.hachure]
            
    return segment_list


#---------------Figure out our next point on the hachure----------------

# Given x/y coordinates for a point on a hachure line, figure out where
# the next point on the line should be. This requires sampling the
# upslope raster and doing a bit of math.

def jump_calculator(coords):
    x,y = coords
    rc = xy_to_rc(coords)
    values = sample_raster(rc) # Get the x and y upslope component values
    
    if values == 0: #if we go out of bounds, stop this line
        return -1
        
    delta_x,delta_y = values #our uphill components
        
    # now normalize this to a vector of jump_distance length
    
    leng = jump_distance / math.hypot(delta_x, delta_y)
    
    delta_x = delta_x * leng
    delta_y = delta_y * leng

    # This bit ensures we're drawing hachures in the correct direction
    
    if (delta_x < 0 and phase == 1) or (delta_y < 0 and phase == 0):
        delta_x *= -1
        delta_y *= -1
    
    return delta_x,delta_y
    
#--Generates new hachures starting at the middle of any given segment---
def hachure_generator(segment_list):

    #First we need the midpoint in each line, to begin our hachure from  
    start_points = []
    
    for segment in segment_list:
        
        midpoint = segment.length / 2
        
        midpoint = segment.geometry.interpolate(midpoint)        
        
        start_points.append(midpoint.asPoint())
    
    #Next loop through the start_points & make hachures
    
    feature_list = []
    
    for coords in start_points:
        
        line_coords = [coords]
        
        for i in range(0,150):
            # this loop is a failsafe in case other checks below fail
            # to stop the hachure when they should
            x,y = line_coords[-1]

            # Figure out the direction to jump next, with -1 == failure
            jump = jump_calculator(line_coords[-1])
            if jump == -1:
                break
            else:
                delta_x,delta_y = jump


            # This confines us to just mostly-horizontal or mostly-
            # vertical hachures, depending on the phase
            if phase == 0 and abs(delta_y) < abs(delta_x):
                break
            elif phase == 1 and abs(delta_y) >= abs(delta_x):
                break
            
            new_y = y + delta_y
            new_x = x + delta_x
      
            # Hachures can bounce back and forth occasionally &
            # should stop. If lines are zig-zagging, every other point
            # will be separated by only a small distance

            if (len(line_coords) > 3 and
                dist(line_coords[-1], line_coords[-3])
                < (jump_distance * 1.5)):
                
            # Snip off the last couple points if we've gone bad:
                del line_coords[-2:]
                break

            line_coords += [(new_x,new_y)]
            
        if len(line_coords) > 1:
            # if we stopped before we even got 2 points, don't bother
            # connecting them into a line.
            feature_list.append(make_lines(line_coords))
    
    return feature_list
    
#==============FUNCTIONS OVER; BEGIN CHECK LINE PREPARATION=============

# Generate horizontal lines covering the raster. The number of lines is
# equal to the raster height.

line_interval = extent.height() / rows
current_y = extent.yMinimum()
x_min = extent.xMinimum()
x_max = extent.xMaximum()

line_geometries = []

for i in range(rows):
    line = QgsLineString([QgsPoint(x_min,current_y),QgsPoint(x_max,current_y)])
    line_geometries.append(QgsGeometry(line))
    current_y += line_interval
    
check_lines_h = [Check_Line(a) for a in line_geometries]

# Repeat the process, this time making *vertical* lines to cover the
# raster width, one per raster pixel

line_interval = extent.width() / cols
current_x = extent.xMinimum()
y_min = extent.yMinimum()
y_max = extent.yMaximum()

line_geometries = []

for i in range(cols):
    line = QgsLineString([QgsPoint(current_x,y_min),QgsPoint(current_x,y_max)])
    line_geometries.append(QgsGeometry(line))
    current_x += line_interval
    
check_lines_v = [Check_Line(a) for a in line_geometries]


#=======BEGIN MAIN LOOP: USE EACH LINE TO GENERATE/CHECK HACHURES=======


# We begin with the first_line function, and check to see if it produced
# any hachures. If not, we go to the next check line and try again.
# Otherwise, it calls subsequent_line() for the rest of the check lines.

for line in check_lines_h:
     if current_hachures:
         subsequent_line(line)
     else:
         first_line(line)

# We sometimes pick up errant duplicates, so let's clean the final list
current_hachures = list(set(current_hachures))

# And filter out any stubs 
current_hachures = [feat for feat in current_hachures
                    if feat.geometry().length() > stubFilter]

# We're now done with phase 0: we've drawn vertical hachures and checked
# them with horizontal lines. We now must repeat this process, drawing
# horizontal hachures that are checked with vertical lines (Phase 1).

# First, we copy our current hachure set into our start_hachures
# This is important because if a hachure is in start_hachures, it won't
# ever be trimmed. It's locked in its current form, and all other
# hachures will be trimmed to give it appropriate space.

start_hachures = current_hachures.copy()

# start_hachures now contains our first set of hachures, plus any
# optional breakline features we included at the start

# now move on to the vertical-basis
phase = 1

# And loop through the hachure generation steps.

for line in check_lines_v:
     if current_hachures:
         subsequent_line(line)
     else:
         first_line(line)

# Clear duplicates again, remove stubs
current_hachures = list(set(current_hachures))
stub_filtered = [feat for feat in current_hachures
                 if feat.geometry().length() > stubFilter]
if breaks_layer_name:
    filtered = [feat for feat in stub_filtered if feat not in breaks]
else:
    filtered = stub_filtered
    
# Add our hachures to the map, plus a length attribute
# the length attribute can be useful for further stub filtering

hachureLayer = QgsVectorLayer('linestring','Hachures','memory')
hachureLayer.setCrs(crs)

field = QgsField('Length', QVariant.Double)
hachureLayer.dataProvider().addAttributes([field])
hachureLayer.updateFields()

for feature in filtered:
    feature.setAttributes([feature.geometry().length()])

with edit(hachureLayer):
    hachureLayer.dataProvider().addFeatures(filtered)
    
instance.addMapLayer(hachureLayer)

# Besides outputting the raw hachures, we can also improve them with
# some quick QGIS tools

# Firstly, we smooth them out:

params = {'LINES_IN': hachureLayer,
    'LINES_OUT':'TEMPORARY_OUTPUT',
    'METHOD':2,
    'SENSITIVITY':3,
    'ITERATIONS':100,
    'PRESERVATION':10,
    'SIGMA':5}

first_smooth = processing.run("sagang:linesmoothing",params)['LINES_OUT']
params['SIGMA'] = 1
params['LINES_IN'] = first_smooth

second_smooth = processing.run("sagang:linesmoothing",params)['LINES_OUT']

# Then, if we had a breaks layer to begin with, I find that it's good to
# help those breaks stand out a bit by clearing their neighborhood.
# We simply buffer them, then delete hachures nearby.

if breaks_layer_name:
    
    params = {'INPUT': instance.mapLayersByName(breaks_layer_name)[0],
        'DISTANCE':5,
        'SEGMENTS':5,
        'END_CAP_STYLE':0,
        'JOIN_STYLE':0,
        'MITER_LIMIT':2,
        'DISSOLVE':False,
        'SEPARATE_DISJOINT':False,
        'OUTPUT':'TEMPORARY_OUTPUT'}
    buffered_breaks = processing.run("native:buffer", params)['OUTPUT']


    params = {'INPUT':QgsVectorLayer(second_smooth),
        'OVERLAY':buffered_breaks,
        'OUTPUT':'TEMPORARY_OUTPUT',
        'GRID_SIZE':None
        }
    differenced = processing.run("native:difference", params)['OUTPUT']
    differenced.setName('Smoothed Hachures')
    instance.addMapLayer(differenced)

else:
    instance.addMapLayer(QgsVectorLayer(second_smooth,'Smoothed Hachures'))

warn_user(0)
