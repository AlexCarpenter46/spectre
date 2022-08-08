#!/usr/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging


def generate_trace(file_destination, output, camera_angle, shift, Ah_Surfaces):
    trace = """
# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10
from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
# create a new 'XDMF Reader'
volumeFilesxmf = XDMFReader(registrationName='VolumeFiles.xmf', FileNames=["""
    trace += "'" + file_destination + "VolumeFiles.xmf'])"
    trace += """
volumeFilesxmf.PointArrayStatus = ['Lapse', 'Shift', 'SpatialRicciScalar']
volumeFilesxmf.GridStatus = ['Evolution']
# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# get the material library
materialLibrary1 = GetMaterialLibrary()
# get display properties
volumeFilesxmfDisplay = GetDisplayProperties(volumeFilesxmf, view=renderView1)
# get color transfer function/color map for 'Lapse'
lapseLUT = GetColorTransferFunction('Lapse')
lapseLUT.RGBPoints = [0.47452330589294434, 0.231373, 0.298039, 0.752941,
 0.7356098592281342, 0.865003, 0.865003, 0.865003, 0.996696412563324, 0.705882,
  0.0156863, 0.14902]
lapseLUT.ScalarRangeInitialized = 1.0
# get opacity transfer function/opacity map for 'Lapse'
lapsePWF = GetOpacityTransferFunction('Lapse')
lapsePWF.Points = [0.47452330589294434, 0.0, 0.5, 0.0, 0.996696412563324,
 1.0, 0.5, 0.0]
lapsePWF.ScalarRangeInitialized = 1
# get animation scene
animationScene1 = GetAnimationScene()
# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()
# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=volumeFilesxmf)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'Lapse']
slice1Display.LookupTable = lapseLUT
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleArray = 'Lapse'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 60.0
slice1Display.SelectScaleArray = 'Lapse'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'Lapse'
slice1Display.GaussianRadius = 3.0
slice1Display.SetScaleArray = ['POINTS', 'Lapse']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Lapse']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.878652036190033, 0.0, 0.5,
 0.0, 0.996696412563324, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.878652036190033, 0.0, 0.5,
 0.0, 0.996696412563324, 1.0, 0.5, 0.0]
# hide data in view
Hide(volumeFilesxmf, renderView1)
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]
# update the view to ensure updated data information
renderView1.Update()
# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)
# Properties modified on slice1
slice1.Triangulatetheslice = 0
# update the view to ensure updated data information
renderView1.Update()
# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=slice1)
calculator1.Function = ''
# Properties modified on calculator1
calculator1.Function = ''
# show data in view
calculator1Display = Show(calculator1, renderView1, 'GeometryRepresentation')
# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = ['POINTS', 'Lapse']
calculator1Display.LookupTable = lapseLUT
calculator1Display.SelectTCoordArray = 'None'
calculator1Display.SelectNormalArray = 'None'
calculator1Display.SelectTangentArray = 'None'
calculator1Display.OSPRayScaleArray = 'Lapse'
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.ScaleFactor = 60.0
calculator1Display.SelectScaleArray = 'Lapse'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.GlyphTableIndexArray = 'Lapse'
calculator1Display.GaussianRadius = 3.0
calculator1Display.SetScaleArray = ['POINTS', 'Lapse']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = ['POINTS', 'Lapse']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
calculator1Display.DataAxesGrid = 'GridAxesRepresentation'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [0.47452330589294434,
 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [0.47452330589294434,
 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0]
# hide data in view
Hide(slice1, renderView1)
# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on calculator1
calculator1.Function = 'abs(SpatialRicciScalar)'
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on calculator1
calculator1.ResultArrayName = 'AbsValueSpatialRicci'
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on lapseLUT
lapseLUT.Discretize = 0
# Apply a preset using its name. Note this may not work as
# expected when presets have duplicate names.
lapseLUT.ApplyPreset('Rainbow Uniform', True)
# invert the transfer function
lapseLUT.InvertTransferFunction()
# create a new 'Warp By Scalar'
warpByScalar1 = WarpByScalar(registrationName='WarpByScalar1',
 Input=calculator1)
warpByScalar1.Scalars = ['POINTS', 'AbsValueSpatialRicci']
# show data in view
warpByScalar1Display = Show(warpByScalar1, renderView1,
 'GeometryRepresentation')
# get color transfer function/color map for 'AbsValueSpatialRicci'
absValueSpatialRicciLUT = GetColorTransferFunction('AbsValueSpatialRicci')
absValueSpatialRicciLUT.RGBPoints = [1.3087053964190598e-15, 0.231373, 0.298039,
 0.752941, 0.5333755612373359, 0.865003, 0.865003, 0.865003, 1.0667511224746704,
  0.705882, 0.0156863, 0.14902]
absValueSpatialRicciLUT.ScalarRangeInitialized = 1.0
# trace defaults for the display properties.
warpByScalar1Display.Representation = 'Surface'
warpByScalar1Display.ColorArrayName = ['POINTS', 'AbsValueSpatialRicci']
warpByScalar1Display.LookupTable = absValueSpatialRicciLUT
warpByScalar1Display.SelectTCoordArray = 'None'
warpByScalar1Display.SelectNormalArray = 'None'
warpByScalar1Display.SelectTangentArray = 'None'
warpByScalar1Display.OSPRayScaleArray = 'AbsValueSpatialRicci'
warpByScalar1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByScalar1Display.SelectOrientationVectors = 'None'
warpByScalar1Display.ScaleFactor = 60.0
warpByScalar1Display.SelectScaleArray = 'AbsValueSpatialRicci'
warpByScalar1Display.GlyphType = 'Arrow'
warpByScalar1Display.GlyphTableIndexArray = 'AbsValueSpatialRicci'
warpByScalar1Display.GaussianRadius = 3.0
warpByScalar1Display.SetScaleArray = ['POINTS', 'AbsValueSpatialRicci']
warpByScalar1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByScalar1Display.OpacityArray = ['POINTS', 'AbsValueSpatialRicci']
warpByScalar1Display.OpacityTransferFunction = 'PiecewiseFunction'
warpByScalar1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByScalar1Display.PolarAxes = 'PolarAxesRepresentation'
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
warpByScalar1Display.ScaleTransferFunction.Points = [1.3087053964190598e-15,
 0.0, 0.5, 0.0, 1.0667511224746704, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
warpByScalar1Display.OpacityTransferFunction.Points = [1.3087053964190598e-15,
 0.0, 0.5, 0.0, 1.0667511224746704, 1.0, 0.5, 0.0]
# hide data in view
Hide(calculator1, renderView1)
# show color bar/color legend
warpByScalar1Display.SetScalarBarVisibility(renderView1, True)
# update the view to ensure updated data information
renderView1.Update()
# get opacity transfer function/opacity map for 'AbsValueSpatialRicci'
absValueSpatialRicciPWF = GetOpacityTransferFunction('AbsValueSpatialRicci')
absValueSpatialRicciPWF.Points = [1.3087053964190598e-15, 0.0, 0.5, 0.0,
 1.0667511224746704, 1.0, 0.5, 0.0]
absValueSpatialRicciPWF.ScalarRangeInitialized = 1
# Properties modified on warpByScalar1
warpByScalar1.Normal = [0.0, 0.0, -1.0]
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on warpByScalar1
warpByScalar1.ScaleFactor = 1.5
# update the view to ensure updated data information
renderView1.Update()
# set scalar coloring
ColorBy(warpByScalar1Display, ('POINTS', 'Lapse'))
# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(absValueSpatialRicciLUT, renderView1)
# rescale color and/or opacity maps used to include current data range
warpByScalar1Display.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
warpByScalar1Display.SetScalarBarVisibility(renderView1, True)
# Properties modified on warpByScalar1Display
warpByScalar1Display.Opacity = 0.8
# Hide orientation axes
renderView1.OrientationAxesVisibility = 0
# hide color bar/color legend
warpByScalar1Display.SetScalarBarVisibility(renderView1, False)
LoadPalette(paletteName='GradientBackground')"""
    if (Ah_Surfaces == 0):
        trace += """
# create a new 'XDMF Reader'
ahAxmf = XDMFReader(registrationName='AhA.xmf', FileNames=["""
        trace += "'" + file_destination + "AhA.xmf'])"
        trace += """
ahAxmf.PointArrayStatus = ['RicciScalar']
ahAxmf.GridStatus = ['Evolution']
# show data in view
ahAxmfDisplay = Show(ahAxmf, renderView1, 'UnstructuredGridRepresentation')
# get color transfer function/color map for 'RicciScalar'
ricciScalarLUT = GetColorTransferFunction('RicciScalar')
ricciScalarLUT.RGBPoints = [1.9988615925761155, 0.231373, 0.298039, 0.752941,
1.999732833173388, 0.865003, 0.865003, 0.865003, 2.0006040737706603,
0.705882, 0.0156863, 0.14902]
ricciScalarLUT.ScalarRangeInitialized = 1.0
# get opacity transfer function/opacity map for 'RicciScalar'
ricciScalarPWF = GetOpacityTransferFunction('RicciScalar')
ricciScalarPWF.Points = [1.9988615925761155, 0.0, 0.5, 0.0, 2.0006040737706603,
1.0, 0.5, 0.0]
ricciScalarPWF.ScalarRangeInitialized = 1
# trace defaults for the display properties.
ahAxmfDisplay.Representation = 'Surface'
ahAxmfDisplay.ColorArrayName = ['POINTS', 'RicciScalar']
ahAxmfDisplay.LookupTable = ricciScalarLUT
ahAxmfDisplay.SelectTCoordArray = 'None'
ahAxmfDisplay.SelectNormalArray = 'None'
ahAxmfDisplay.SelectTangentArray = 'None'
ahAxmfDisplay.OSPRayScaleArray = 'RicciScalar'
ahAxmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ahAxmfDisplay.SelectOrientationVectors = 'None'
ahAxmfDisplay.ScaleFactor = 0.09697864643297437
ahAxmfDisplay.SelectScaleArray = 'RicciScalar'
ahAxmfDisplay.GlyphType = 'Arrow'
ahAxmfDisplay.GlyphTableIndexArray = 'RicciScalar'
ahAxmfDisplay.GaussianRadius = 0.004848932321648718
ahAxmfDisplay.SetScaleArray = ['POINTS', 'RicciScalar']
ahAxmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ahAxmfDisplay.OpacityArray = ['POINTS', 'RicciScalar']
ahAxmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ahAxmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
ahAxmfDisplay.PolarAxes = 'PolarAxesRepresentation'
ahAxmfDisplay.ScalarOpacityFunction = ricciScalarPWF
ahAxmfDisplay.ScalarOpacityUnitDistance = 0.23102038037530326
ahAxmfDisplay.OpacityArrayName = ['POINTS', 'RicciScalar']
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
ahAxmfDisplay.ScaleTransferFunction.Points = [1.9988615925761155, 0.0, 0.5, 0.0,
 2.0006040737706603, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
ahAxmfDisplay.OpacityTransferFunction.Points = [1.9988615925761155, 0.0, 0.5,
 0.0, 2.0006040737706603, 1.0, 0.5, 0.0]
# show color bar/color legend
ahAxmfDisplay.SetScalarBarVisibility(renderView1, True)
# update the view to ensure updated data information
renderView1.Update()
# create a new 'Transform'
transform1 = Transform(registrationName='Transform1', Input=ahAxmf)
transform1.Transform = 'Transform'
# show data in view
transform1Display = Show(transform1, renderView1,
 'UnstructuredGridRepresentation')
# trace defaults for the display properties.
transform1Display.Representation = 'Surface'
transform1Display.ColorArrayName = ['POINTS', 'RicciScalar']
transform1Display.LookupTable = ricciScalarLUT
transform1Display.SelectTCoordArray = 'None'
transform1Display.SelectNormalArray = 'None'
transform1Display.SelectTangentArray = 'None'
transform1Display.OSPRayScaleArray = 'RicciScalar'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 0.09697864643297437
transform1Display.SelectScaleArray = 'RicciScalar'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'RicciScalar'
transform1Display.GaussianRadius = 0.004848932321648718
transform1Display.SetScaleArray = ['POINTS', 'RicciScalar']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', 'RicciScalar']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.ScalarOpacityFunction = ricciScalarPWF
transform1Display.ScalarOpacityUnitDistance = 0.23102038037530326
transform1Display.OpacityArrayName = ['POINTS', 'RicciScalar']
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform1Display.ScaleTransferFunction.Points = [1.9988615925761155,
 0.0, 0.5, 0.0, 2.0006040737706603, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform1Display.OpacityTransferFunction.Points = [1.9988615925761155,
 0.0, 0.5, 0.0, 2.0006040737706603, 1.0, 0.5, 0.0]
# hide data in view
Hide(ahAxmf, renderView1)
# show color bar/color legend
transform1Display.SetScalarBarVisibility(renderView1, True)
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on transform1.Transform
transform1.Transform.Translate = [0.0, 0.0, 2.5]
# update the view to ensure updated data information
renderView1.Update()
# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=transform1.Transform)
# turn off scalar coloring
ColorBy(transform1Display, None)
# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(ricciScalarLUT, renderView1)
# change solid color
transform1Display.AmbientColor = [0.0, 0.0, 0.0]
transform1Display.DiffuseColor = [0.0, 0.0, 0.0]
# create a new 'XDMF Reader'
ahBxmf = XDMFReader(registrationName='AhB.xmf', FileNames=["""
        trace += "'" + file_destination + "AhB.xmf'])"
        trace += """
ahBxmf.PointArrayStatus = ['RicciScalar']
ahBxmf.GridStatus = ['Evolution']
# show data in view
ahBxmfDisplay = Show(ahBxmf, renderView1, 'UnstructuredGridRepresentation')
# trace defaults for the display properties.
ahBxmfDisplay.Representation = 'Surface'
ahBxmfDisplay.ColorArrayName = ['POINTS', 'RicciScalar']
ahBxmfDisplay.LookupTable = ricciScalarLUT
ahBxmfDisplay.SelectTCoordArray = 'None'
ahBxmfDisplay.SelectNormalArray = 'None'
ahBxmfDisplay.SelectTangentArray = 'None'
ahBxmfDisplay.OSPRayScaleArray = 'RicciScalar'
ahBxmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ahBxmfDisplay.SelectOrientationVectors = 'None'
ahBxmfDisplay.ScaleFactor = 0.09697864697486241
ahBxmfDisplay.SelectScaleArray = 'RicciScalar'
ahBxmfDisplay.GlyphType = 'Arrow'
ahBxmfDisplay.GlyphTableIndexArray = 'RicciScalar'
ahBxmfDisplay.GaussianRadius = 0.0048489323487431206
ahBxmfDisplay.SetScaleArray = ['POINTS', 'RicciScalar']
ahBxmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ahBxmfDisplay.OpacityArray = ['POINTS', 'RicciScalar']
ahBxmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ahBxmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
ahBxmfDisplay.PolarAxes = 'PolarAxesRepresentation'
ahBxmfDisplay.ScalarOpacityFunction = ricciScalarPWF
ahBxmfDisplay.ScalarOpacityUnitDistance = 0.23102038188159646
ahBxmfDisplay.OpacityArrayName = ['POINTS', 'RicciScalar']
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
ahBxmfDisplay.ScaleTransferFunction.Points = [1.998871756148613, 0.0, 0.5, 0.0,
 2.000595342018353, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
ahBxmfDisplay.OpacityTransferFunction.Points = [1.998871756148613, 0.0, 0.5,
 0.0, 2.000595342018353, 1.0, 0.5, 0.0]
# show color bar/color legend
ahBxmfDisplay.SetScalarBarVisibility(renderView1, True)
# update the view to ensure updated data information
renderView1.Update()
# create a new 'Transform'
transform2 = Transform(registrationName='Transform2', Input=ahBxmf)
transform2.Transform = 'Transform'
# show data in view
transform2Display = Show(transform2, renderView1,
 'UnstructuredGridRepresentation')
# trace defaults for the display properties.
transform2Display.Representation = 'Surface'
transform2Display.ColorArrayName = ['POINTS', 'RicciScalar']
transform2Display.LookupTable = ricciScalarLUT
transform2Display.SelectTCoordArray = 'None'
transform2Display.SelectNormalArray = 'None'
transform2Display.SelectTangentArray = 'None'
transform2Display.OSPRayScaleArray = 'RicciScalar'
transform2Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform2Display.SelectOrientationVectors = 'None'
transform2Display.ScaleFactor = 0.09697864697486241
transform2Display.SelectScaleArray = 'RicciScalar'
transform2Display.GlyphType = 'Arrow'
transform2Display.GlyphTableIndexArray = 'RicciScalar'
transform2Display.GaussianRadius = 0.0048489323487431206
transform2Display.SetScaleArray = ['POINTS', 'RicciScalar']
transform2Display.ScaleTransferFunction = 'PiecewiseFunction'
transform2Display.OpacityArray = ['POINTS', 'RicciScalar']
transform2Display.OpacityTransferFunction = 'PiecewiseFunction'
transform2Display.DataAxesGrid = 'GridAxesRepresentation'
transform2Display.PolarAxes = 'PolarAxesRepresentation'
transform2Display.ScalarOpacityFunction = ricciScalarPWF
transform2Display.ScalarOpacityUnitDistance = 0.23102038188159646
transform2Display.OpacityArrayName = ['POINTS', 'RicciScalar']
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
transform2Display.ScaleTransferFunction.Points = [1.998871756148613, 0.0, 0.5,
 0.0, 2.000595342018353, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
transform2Display.OpacityTransferFunction.Points = [1.998871756148613, 0.0, 0.5,
 0.0, 2.000595342018353, 1.0, 0.5, 0.0]
# hide data in view
Hide(ahBxmf, renderView1)
# show color bar/color legend
transform2Display.SetScalarBarVisibility(renderView1, True)
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on transform2.Transform
transform2.Transform.Translate = [0.0, 0.0, 2.5]
# update the view to ensure updated data information
renderView1.Update()
# turn off scalar coloring
ColorBy(transform2Display, None)
# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(ricciScalarLUT, renderView1)
# change solid color
transform2Display.AmbientColor = [0.0, 0.0, 0.0]
transform2Display.DiffuseColor = [0.0, 0.0, 0.0]
# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=transform2.Transform)"""
    trace += """
# set active source
SetActiveSource(warpByScalar1)
# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(registrationName='AnnotateTimeFilter1',
 Input=warpByScalar1)
# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1,
 'TextSourceRepresentation')
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = 'Time: {time:} M'
# update the view to ensure updated data information
renderView1.Update()
# Properties modified on annotateTimeFilter1Display
annotateTimeFilter1Display.FontSize = 30
# update the view to ensure updated data information
renderView1.Update()"""
    if (shift == 0):
        trace += """
# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=warpByScalar1,
GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'ShiftVector']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.5
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'Every Nth Point'
glyph1.Stride = 60
# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'AbsValueSpatialRicci']
glyph1Display.LookupTable = absValueSpatialRicciLUT
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'AbsValueSpatialRicci'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'ShiftVector'
glyph1Display.ScaleFactor = 57.59131774902344
glyph1Display.SelectScaleArray = 'AbsValueSpatialRicci'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'AbsValueSpatialRicci'
glyph1Display.GaussianRadius = 2.879565887451172
glyph1Display.SetScaleArray = ['POINTS', 'AbsValueSpatialRicci']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'AbsValueSpatialRicci']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [1.7175290703552726e-11, 0.0, 0.5,
 0.0, 1.058470606803894, 1.0, 0.5, 0.0]
# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [1.7175290703552726e-11, 0.0,
 0.5, 0.0, 1.058470606803894, 1.0, 0.5, 0.0]
# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)
# update the view to ensure updated data information
renderView1.Update()
# turn off scalar coloring
ColorBy(glyph1Display, None)
# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(absValueSpatialRicciLUT, renderView1)
# get the time-keeper
timeKeeper1 = GetTimeKeeper()
# set active source
SetActiveSource(warpByScalar1)
# set active source
SetActiveSource(glyph1)"""
    trace += """
#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================
# get layout
layout1 = GetLayout()
#--------------------------------
# saving layout sizes for layouts
#layout/tab size in pixels
layout1.SetSize(1810, 1156)
#-----------------------------------
# camera placement settings"""
    if (camera_angle == 1):
        trace += """
renderView1.CameraPosition = [0.0, 0.0, 36.90869716569761]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.6894899550131899]
renderView1.CameraViewUp = [0, 1, 0]
renderView1.CameraParallelScale = 424.27024700303446"""
    elif (camera_angle == 2):
        trace += """
renderView1.CameraPosition = [-30.093571816984163, -18.55482256667294,
 8.558827016411096]
renderView1.CameraFocalPoint = [6.698318283229103e-17, 4.073218385661855e-17,
 0.6894899550131903]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 424.27024700303446"""
    else:
        trace += """
renderView1.CameraPosition = [-16.82050334149194, -1.8559824362222395,
 2.6797079790683345]
renderView1.CameraFocalPoint = [0.0836179084847767, -1.0057516274428049,
 0.22771979385109334]
renderView1.CameraViewUp = [0, 0, 1]
renderView1.CameraParallelScale = 519.6152422706632"""
    with open(output + ".py", "w") as py_file:
        py_file.write(trace)


def parse_args():
    """
    Parse the command line arguments
    """
    import argparse as ap
    parser = ap.ArgumentParser(
        description="Generate trace file for visualizing SpECTRE data.",
        formatter_class=ap.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--file-destination',
        required=True,
        help="The common prefix of the H5 volume files to load, excluding "
        "the node number integer(s)")
    parser.add_argument('--output',
                        '-o',
                        required=True,
                        help="Output file name, a py extension will be added")
    parser.add_argument('--camera-angle',
                        default=0,
                        type=int,
                        help="Which camera angle you'd like")
    parser.add_argument('--shift',
                        '-s',
                        default=0,
                        type=int,
                        help="Disable output of shift arrows.")
    parser.add_argument('--Ah-Surfaces',
                        '-a',
                        default=0,
                        type=int,
                        help="Disable output of Ah Surfaces.")
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    input_args = parse_args()
    generate_trace(**vars(input_args))
