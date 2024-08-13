# trace generated using paraview version 5.12.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
from numpy import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

n = 167
# create a new 'AMReX/BoxLib Grid Reader'
plot_filesseries = AMReXBoxLibGridReader(registrationName='plot_files.series', FileNames=['/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/plot_files.series'])


# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on plot_filesseries
plot_filesseries.CellArrayStatus = ['qc', 'qrain', 'rain_accum']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
plot_filesseriesDisplay = Show(plot_filesseries, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
plot_filesseriesDisplay.Representation = 'Outline'

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=plot_filesseries)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Outline'

# hide data in view
Hide(plot_filesseries, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=cellDatatoPointData1)

# Properties modified on contour1
contour1.Isosurfaces = [1e-05]

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(cellDatatoPointData1, renderView1)

for i in arange(0, n, 1):

    print("Doing scene %d of %d"%(i,n))


    if(i<=9):
        filename = "/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/supercell_3d_00%d.stl"%i
    elif(i<=99):
        filename = "/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/supercell_3d_0%d.stl"%i
    elif(i<=999):
        filename = "/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/supercell_3d_%d.stl"%i

    print("filename is %s", filename)

# save data
    SaveData(filename, proxy=contour1)

    animationScene1.GoToNext()
