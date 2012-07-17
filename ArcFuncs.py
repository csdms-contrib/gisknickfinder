#This is a function that I will use to do all of my Arc work
#Started April 11, 2012
import arcpy, math, numpy as np, operator
from arcpy import env
from arcpy.sa import *
env.workspace = "C:/GIS/Bijou/Python/Hcpyproj"
arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = 1

def WatShedfunc(Watershed, WatershedMask, Fillgrid, FlowdirArc, FlowAcc):
    print "Starting Watershed Module"
    ###################################################
    #I will begin by extracting my DEM to the watershed and creating some of the watershed characteristics
    ExtractWS = ExtractByMask(Watershed, WatershedMask)
    #Next I create a flow direction grid in arc using the hydrology tools
    outFill = Fill(ExtractWS)
    outFill.save(Fillgrid)
    Newflowdir=FlowDirection(outFill)
    Newflowdir.save(FlowdirArc)
    outFlowAccumulation = FlowAccumulation(Newflowdir)
    outFlowAccumulation.save(FlowAcc)
    return ExtractWS, outFill, Newflowdir, outFlowAccumulation
    print "Watershed Module has a problem"

def dxfunc(Newflowdir, dxRaster):
    #Now I'm going to reclassify the flow direction grid to get a value of dx
    #dxRas = Con(((Newflowdir == 1) | (Newflowdir == 3) | (Newflowdir == 5) | (Newflowdir == 7)), 1, math.sqrt(2))
    dxRas = Con(((Newflowdir == 1) | (Newflowdir == 64) | (Newflowdir == 16) | (Newflowdir == 4)), 1, math.sqrt(2))
    dxRas.save(dxRaster)
    return dxRas

def findhcs(outFill, Curvpro, CurvThres, FlowAccThres, outFlowAccumulation, HClocations, Rastpoly, Dissolvepoly, HCpnts):
    outCurve = Curvature(outFill, 1, Curvpro)
    Raster(Curvpro).save(env.workspace + "/curvtest")
    #Identifying headcuts below
    #HcLocations=((Raster(Curvpro)<-15) & (outFlowAccumulation>5000))
    HcLocations=((Raster(Curvpro)<CurvThres) & (outFlowAccumulation>FlowAccThres))
    HcLocations.save(HClocations)
    #Now I have to clean up the pixels so I only have one point at each headcut
    #First use region group function
    RegGroup=RegionGroup(HcLocations, "EIGHT", "WITHIN", "ADD_LINK")
    #Next I need to use Raster to polygon
    Rast2Poly=arcpy.RasterToPolygon_conversion(RegGroup, Rastpoly, "NO_SIMPLIFY", "Value")
    #Now I have to Dissolve the polygon
    Dissolve=arcpy.Dissolve_management(Rastpoly, Dissolvepoly,"GRIDCODE", "", "MULTI_PART")
    #I think I found an error here. I noticed that if I only have one headcut, when I use the dissolve function it will give me
    #a point for my headcut, and a point that represents the rest of the area of the watershed. I obviously want to delete the
    #point for the rest of the watershed area, and just keep the headcut point.
    #I add a column to the dissolve feature class called area, calculate the area for each field in the dissolve polygon, then
    #delete any field that is > 5 m^2 (this should be a large enough buffer that I don't accidently delete my headcut. 
    #The last thing I have to do is create a final shapefile for the headcut points
    arcpy.AddField_management(Dissolve, "Area", "LONG", 9, "", "", "", "NULLABLE")
    arcpy.CalculateField_management(Dissolve, "Area",'!shape.area!', "PYTHON")
    #Now I use the update cursor to delete any really large areas, so that I'm only left with my headcut
    rows = arcpy.UpdateCursor(Dissolve)
    for row in rows:
        if row.Area > 5:
            rows.deleteRow(row)
    HCs=arcpy.FeatureToPoint_management(Dissolve, HCpnts)
    return outCurve, HcLocations, RegGroup, Rast2Poly, Dissolve, HCs

def leastcostfunc(HCs, SpatialR, Fillgrid, Newflowdir, ncostpath):
    rows = arcpy.SearchCursor(HCs)
    #Now I'm going to have to loop through the shapefile, take one point out at a time, and turn that into a new shapefile
    count=1
    for row in rows:
        #I need to make a shapefile and add just one of the headcut points to that.
        pnt=row.shape.getPart(0)
        point = arcpy.Point()
        #Getting the X & Y coordiantes of each point in the shape file 
        point.X=pnt.X
        point.Y=pnt.Y
        spatial_reference = arcpy.SpatialReference(SpatialR)
        template = HCs
        Dumshape=arcpy.CreateFeatureclass_management(env.workspace, "dummypoint"+str(count)+".shp", "Point", template, "DISABLED", "DISABLED", spatial_reference)
        #Now I need to make an insert cursor to add the point to my new shapefile
        rowz = arcpy.InsertCursor(Dumshape)
        newpnt = rowz.newRow()
        newpnt.shape=point
        rowz.insertRow(newpnt)
        del rowz, newpnt
        #Now I'm going to create the least cost path to the point in the new shape file that I just created
        outCostPath = CostPath(Dumshape, Fillgrid, Newflowdir, "EACH_CELL","FID")
        #The least cost path has 2 vaules, I want to make it all ones, so I do that below using the setnull method
        Newcostpath = SetNull(outCostPath, 1, "VALUE < 0")
        Newcostpath.save(ncostpath+str(count))
        count=count+1
    del row, rows
    return count

def AreaDistance(HCpnts, FlowAcc, DAreaatHC):
    HeadcutRast=arcpy.FeatureToRaster_conversion(HCpnts, "GRIDCODE", HCraster, 1)
    DAatHC = SetNull(HeadcutRast, FlowAcc, "VALUE > 0")
    DAatHC.save(DAreaatHC)

def hcpathvectors(FlowAcc, dxRas, ncostpath, count,Fillgrid):
    print "Getting the Area and dx Arrays"
    #I am going to convert the rasters to numpy grids
    AreaArray=arcpy.RasterToNumPyArray(FlowAcc)
    AreaArray[AreaArray <= 0] = 0
    dxArray=arcpy.RasterToNumPyArray(dxRas)
    dxArray[dxArray <= 0] = 0
    ElevArray=arcpy.RasterToNumPyArray(Fillgrid)
    ElevArray[ElevArray <= 0] = 0
    #Below I get the rows and column numbers of the array (they should all be the same size)
    [r,c]=AreaArray.shape
    #Now I initialize lots of python lists
    CostVec=[]
    AreaVec=[]
    dxVec=[]
    elevVec=[]
    CostArray2=[]
    AreaArray2=[]
    dxArray2=[]
    elevArray2=[]
    #Below I march through each row and column, and if my cost vector
    #is equal to 1, then I grab the drainage area array and I grab the
    #flow direction grid value
    for i in range(1,count):
        Kostpath=ncostpath+str(i)
        CostArray = arcpy.RasterToNumPyArray(Kostpath)
        test=np.sum(CostArray)
        print "Cost Array sum is: "+str(test)
        #I use the CostArray to zero out values outside of the path to the headcut
        AreaVec=AreaArray[CostArray == 1]
        dxVec=dxArray[CostArray == 1]
        elevVec=ElevArray[CostArray == 1]
        lendx=len(dxVec)
        print "dx length is: "+str(lendx)
        AreaArray2.append(AreaVec)
        dxArray2.append(dxVec)
        elevArray2.append(elevVec)
        CostVec=[]
        AreaVec=[]
        dxVec=[]
        elevVec=[]
    return AreaArray2, dxArray2, elevArray2

def hcpathvectors2(FlowAcc, dxRas, ncostpath, count,Fillgrid):
    print "Getting the Area and dx Arrays"
    #I am going to convert the rasters to numpy grids
    AreaArray=arcpy.RasterToNumPyArray(FlowAcc)
    AreaArray[AreaArray <= 0] = 0
    dxArray=arcpy.RasterToNumPyArray(dxRas)
    dxArray[dxArray <= 0] = 0
    ElevArray=arcpy.RasterToNumPyArray(Fillgrid)
    ElevArray[ElevArray <= 0] = 0
    #Below I get the rows and column numbers of the array (they should all be the same size)
    [r,c]=AreaArray.shape
    #Now I initialize lots of python lists
    CostVec=[]
    AreaVec=[]
    dxVec=[]
    elevVec=[]
    CostArray2=[]
    AreaArray2=[]
    dxArray2=[]
    elevArray2=[]
    #Below I march through each row and column, and if my cost vector
    #is equal to 1, then I grab the drainage area array and I grab the
    #flow direction grid value
    for i in range(1,count):
        Kostpath=ncostpath+str(i)
        CostArray = arcpy.RasterToNumPyArray(Kostpath)
        test=np.sum(CostArray)
        print "Cost Array sum is: "+str(test)
        #I use the CostArray to zero out values outside of the path to the headcut
        AreaVec=AreaArray[CostArray == 1]
        dxVec=dxArray[CostArray == 1]
        elevVec=ElevArray[CostArray == 1]
        RealAreaVec=list(AreaVec)
        #Areaind=range(len(RealAreaVec))#These should be the current indices of Area
        #Areaind.sort(lambda x,y:RealAreaVec[x]-RealAreaVec[y])#This is where I sort 
        Areaind=[ p for (p,q) in sorted(enumerate(RealAreaVec), key=operator.itemgetter(1))]
        AreaVec=[ AreaVec[j] for j in Areaind ]
        dxVec=[ dxVec[j] for j in Areaind ]
        elevVec=[ elevVec[j] for j in Areaind ]
        lendx=len(dxVec)
        print "dx length is: "+str(lendx)
        AreaArray2.append(AreaVec)
        dxArray2.append(dxVec)
        elevArray2.append(elevVec)
        CostVec=[]
        AreaVec=[]
        dxVec=[]
        elevVec=[]
    return AreaArray2, dxArray2, elevArray2

def getxyofstr(CountStTop, indstrcostpath, StrPathPnts):
    AllStreamXYs=[]
    for i in range(1,CountStTop):
        StreamXYs=[]
        inputrast=indstrcostpath+str(i)
        arcpy.RasterToPoint_conversion(inputrast, StrPathPnts+str(i)+".shp", "VALUE")
        Poly=StrPathPnts+str(i)+".shp"
        #arcpy.Copy_management(in_data, in_features)
        arcpy.AddXY_management(Poly)
        #Now I'm going to grab the X,Y values from the shape file and put them in a python list
        rows = arcpy.SearchCursor(Poly)
        for row in rows:
            LocalXY=[]
            pnt=row.shape.getPart(0)
            point = arcpy.Point()
            #Getting the X & Y coordiantes of each point in the shape file 
            point.X=pnt.X
            LocalXY.append(point.X)
            point.Y=pnt.Y
            LocalXY.append(point.Y)
            StreamXYs.append(LocalXY)
        AllStreamXYs.append(StreamXYs)
    return AllStreamXYs
            