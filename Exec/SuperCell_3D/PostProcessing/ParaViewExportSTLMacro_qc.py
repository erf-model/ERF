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
plot_filesseries.CellArrayStatus = ['qc', 'rain_accum']

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

for i in arange(5, n, 1):

	print("Doing scene %d of %d"%(i,n)) 

	
	if(i<=9):
		filename = "/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/STLFiles/qr/supercell_3d_00%d.stl"%i
	elif(i<=99):
		filename = "/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/STLFiles/qr/supercell_3d_0%d.stl"%i
	elif(i<=999):
		filename = "/pscratch/sd/n/nataraj2/ERF/SuperCell_3D/SolnFiles/4th_Standard_alpha66_YOpenXPer/STLFiles/qr/supercell_3d_%d.stl"%i

	print("filename is %s", filename)

# save data
	SaveData(filename, proxy=contour1)

	animationScene1.GoToNext()

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
#layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
#layout1.SetSize(1918, 1074)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
#renderView1.CameraPosition = [43318.13587266245, -345293.42631722393, 60345.64466193015]
#renderView1.CameraFocalPoint = [0.0, 0.0, 12000.0]
#renderView1.CameraViewUp = [-0.10221585482631135, 0.12534642745465943, 0.9868334166142126]
#renderView1.CameraParallelScale = 90934.04203047394


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
