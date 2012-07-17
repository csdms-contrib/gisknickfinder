#Name: GISKnickFinder.  
#This python code can be used to find knickpoints and extract information about streams, it utilizes built-in functions of ArcGIS.
#Copyright (C) 2012  Francis Rengers
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; version 2

#You should have received a copy of the GNU General Public License
#along with this program; if not, write to the Free Software
#Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


import arcpy, math, numpy as np
from arcpy import env
from arcpy.sa import *
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import ArcFuncs
reload(ArcFuncs)

env.workspace = "C:/GIS/Folder"
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = 1

#########################################################
#Input
Watershed = "Path to DEM Watershed (this must be larger than the size of your watershed shapefile"
WatershedMask = env.workspace+"Watershed.shp"
Streamtops=env.workspace+"StreamTop.shp" #This is a point file for the top of each of the streams that you are interested in
SpatialR="C:/Program Files (x86)/ArcGIS/Coordinate Systems/Projected Coordinate Systems/UTM/NAD 1983/NAD 1983 UTM Zone 13N.prj" #This is your projection (I left mine in here just so you could see what one looks like)
CurvThres=-15 #This is the curvature value that highlights your knickpoint. If this value is too high you will have more knickpoints than actually exist. If it is too low then you won't find all of the knickpoints. So you should play around with it.
FlowAccThres=5000 #This is the drainage area threshold. If it is very small, you will highlight knickpoints that are not in your channel. If it is too high then you will miss knickpoints in your channel.
##############################################
#Output
Curvpro=env.workspace+"/curvpropy"
Rastpoly=env.workspace+"/Rastpoly.shp"
Dissolvepoly=env.workspace+"/Dissolvpoly.shp"
HCpnts=env.workspace+"/HCpnts.shp"
FlowdirArc=env.workspace+ "/FlowdirArc"
Costpath = env.workspace+"/costtohcs.shp"
Fillgrid=env.workspace+"/FillWS"
FlowAcc=env.workspace+"/FlowAccArc"
dxRaster=env.workspace+"/dxRas"
HClocations=env.workspace +"/HCLocations"
ncostpath=env.workspace+"/outcost"
HCraster=env.workspace+"/HCraster"
DAreaatHC=env.workspace+"/DAreaatHC"
StrPathPnts=env.workspace+"/StrPathPnts"
indstrcostpath=env.workspace+"/strpath"

try:
    ########################################################################
    #ARCGIS OPERATIONS
    #First I'll clip the watershed, fill the DEM, create the flow direction, and flow accumulation grids
    Watershedstuff=ArcFuncs.WatShedfunc(Watershed, WatershedMask, Fillgrid, FlowdirArc, FlowAcc)

    #Next I take the flow direction grid and recode it to reflect distance 
    #dxgrid=ArcFuncs.dxfunc(Newflowdir, dxRaster)
    dxgrid=ArcFuncs.dxfunc(Raster(FlowdirArc), dxRaster)

    #In this section I will obtain the locations of the headcuts and I will simplify them to make sure that I just have one point
    #at each headcut
    #Allhcs=ArcFuncs.findhcs(outFill, Curvpro, outFlowAccumulation, HClocations, Rastpoly, Dissolvepoly, HCpnts)
    Allhcs=ArcFuncs.findhcs(Raster(Fillgrid), Curvpro, CurvThres, FlowAccThres, Raster(FlowAcc), HClocations, Rastpoly, Dissolvepoly, HCpnts)
    
    #Now I will try to calculate the least cost path to all of the headcuts.
    #Make a search cursor on the shapefile with all of the headcut locations
    #LCostPth=ArcFuncs.leastcostfunc(HCs, SpatialR, Fillgrid, Newflowdir, ncostpath)
    CountLCostPth=ArcFuncs.leastcostfunc(HCpnts, SpatialR, Fillgrid, Raster(FlowdirArc), ncostpath)
    #Count=np.savetxt('Count.txt',LCostPth)

    #I need to calculate the least cost path to the tops of all of the streams too. I will use this to
    #check the K and m values that I obtain.
    CountStTop=ArcFuncs.leastcostfunc(Streamtops, SpatialR, Fillgrid, Raster(FlowdirArc),indstrcostpath) 
    
    ######################################################################
    #PYTHON OPERATIONS
    #Now I create 3 python lists, one for the least cost path
    #this is an array of all ones, the second array is the drainage area corresponding
    #to each cell in the path, and the third array is the dx array showing the distance to the next cell.

    #First I do this just for the headcuts
    PathArrays=ArcFuncs.hcpathvectors2(FlowAcc, dxRaster, ncostpath, CountLCostPth, Fillgrid)
    areaarray=PathArrays[0]
    dxarray=PathArrays[1]
    for i in range(0,len(areaarray)):
        savearea=np.asarray(areaarray[i],dtype='float')
        savedx=np.asarray(dxarray[i],dtype='float')
        np.savetxt(env.workspace+'/AreaArray'+str(i)+'.txt',savearea)
        np.savetxt(env.workspace+'/dxArray'+str(i)+'.txt',savedx)        

    #Next I do this just for the stream paths
    StrPathArrays=ArcFuncs.hcpathvectors2(FlowAcc, dxRaster, indstrcostpath, CountStTop, Fillgrid)
    strareaarray=StrPathArrays[0]
    strdxarray=StrPathArrays[1]
    strelevarray=StrPathArrays[2]
    for i in range(0,len(strareaarray)):
        savestrarea=np.asarray(strareaarray[i],dtype='float')
        savestrdx=np.asarray(strdxarray[i],dtype='float')
        savestrelev=np.asarray(strelevarray[i],dtype='float')
        np.savetxt(env.workspace+'/StrAreaArray'+str(i)+'.txt',savestrarea)
        np.savetxt(env.workspace+'/StrdxArray'+str(i)+'.txt',savestrdx)
        np.savetxt(env.workspace+'/StrElevArray'+str(i)+'.txt',savestrelev)

    #Now I want to grab the X, Y coordinates of each cell in the raster.
    StrXY=ArcFuncs.getxyofstr(CountStTop, indstrcostpath, StrPathPnts)
    for i in range(0,len(StrXY)):
        savestrXY=np.asarray(StrXY[i],dtype='float')
        np.savetxt(env.workspace+'/StrXY'+str(i)+'.txt',savestrXY)
        
    print "It worked!"
except:
    print "Wah Wah"
    #print arcpy.GetMessage()
    
