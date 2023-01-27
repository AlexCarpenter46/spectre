#!/usr/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import click
import glob
import logging
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


def generate_video(file_destination, volume_xmf, aha_xmf, ahb_xmf,
                   output_destination, camera_angle):
    """Generate Pictures from XMF files for BBH Visualizations

    Generate pictures from BBH runs using the XMF files generated using
    generate-xdmf. For files to be read properly, the XMF files being pointed to
    should be in the same directory as the Volume.h5 and Surfaces.h5 files.

    To splice all the pictures into a video, look into FFmpeg"""

    # create a new 'XDMF Reader'
    volulmeFilesxmf = XDMFReader(
        registrationName=volume_xmf,
        FileNames=[file_destination + '/' + volume_xmf])
    volulmeFilesxmf.PointArrayStatus = ['Lapse', 'SpatialRicciScalar']
    # get animation scene
    animationScene1 = GetAnimationScene()
    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # show data in view
    volulmeFilesxmfDisplay = Show(volulmeFilesxmf, renderView1,
                                  'UnstructuredGridRepresentation')
    # get color transfer function/color map for 'Lapse'
    lapseLUT = GetColorTransferFunction('Lapse')
    # get opacity transfer function/opacity map for 'Lapse'
    lapsePWF = GetOpacityTransferFunction('Lapse')
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    volulmeFilesxmfDisplay.ScaleTransferFunction.Points = [
        0.4745234251022339, 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0
    ]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    volulmeFilesxmfDisplay.OpacityTransferFunction.Points = [
        0.4745234251022339, 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0
    ]
    # reset view to fit data
    renderView1.ResetCamera(False)
    # get the material library
    materialLibrary1 = GetMaterialLibrary()
    # show color bar/color legend
    volulmeFilesxmfDisplay.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # create a new 'Slice'
    slice1 = Slice(registrationName='Slice1', Input=volulmeFilesxmf)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]
    # Properties modified on slice1.SliceType
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]
    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice1Display.ScaleTransferFunction.Points = [
        0.4745234251022339, 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0
    ]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice1Display.OpacityTransferFunction.Points = [
        0.4745234251022339, 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0
    ]
    # hide data in view
    Hide(volulmeFilesxmf, renderView1)
    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # Properties modified on lapseLUT
    lapseLUT.Discretize = 0
    # Properties modified on slice1
    slice1.Triangulatetheslice = 0
    # update the view to ensure updated data information
    renderView1.Update()
    # create a new 'Warp By Scalar'
    warpByScalar1 = WarpByScalar(registrationName='WarpByScalar1',
                                 Input=slice1)
    warpByScalar1.Scalars = ['POINTS', 'Lapse']
    # Properties modified on warpByScalar1
    warpByScalar1.Scalars = ['POINTS', 'SpatialRicciScalar']
    warpByScalar1.ScaleFactor = 2.5
    warpByScalar1.Normal = [0.0, 0.0, -1.0]
    # show data in view
    warpByScalar1Display = Show(warpByScalar1, renderView1,
                                'GeometryRepresentation')
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    warpByScalar1Display.ScaleTransferFunction.Points = [
        0.4745234251022339, 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0
    ]
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    warpByScalar1Display.OpacityTransferFunction.Points = [
        0.4745234251022339, 0.0, 0.5, 0.0, 0.996696412563324, 1.0, 0.5, 0.0
    ]
    # hide data in view
    Hide(slice1, renderView1)
    # show color bar/color legend
    warpByScalar1Display.SetScalarBarVisibility(renderView1, True)
    # update the view to ensure updated data information
    renderView1.Update()
    # Apply a preset using its name.
    lapseLUT.ApplyPreset('Rainbow Uniform', True)
    # invert the transfer function
    lapseLUT.InvertTransferFunction()
    if (aha_xmf != ''):
        # create a new 'XDMF Reader'
        ahAxmf = XDMFReader(registrationName=aha_xmf,
                            FileNames=[file_destination + '/' + aha_xmf])
        ahAxmf.PointArrayStatus = ['RicciScalar']
        # show data in view
        ahAxmfDisplay = Show(ahAxmf, renderView1,
                             'UnstructuredGridRepresentation')
        # get color transfer function/color map for 'RicciScalar'
        ricciScalarLUT = GetColorTransferFunction('RicciScalar')
        # get opacity transfer function/opacity map for 'RicciScalar'
        ricciScalarPWF = GetOpacityTransferFunction('RicciScalar')
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        ahAxmfDisplay.ScaleTransferFunction.Points = [
            1.9988708955637904, 0.0, 0.5, 0.0, 2.000606718406889, 1.0, 0.5, 0.0
        ]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        ahAxmfDisplay.OpacityTransferFunction.Points = [
            1.9988708955637904, 0.0, 0.5, 0.0, 2.000606718406889, 1.0, 0.5, 0.0
        ]
        # show color bar/color legend
        ahAxmfDisplay.SetScalarBarVisibility(renderView1, True)
        # update the view to ensure updated data information
        renderView1.Update()
        # create a new 'Transform'
        transform1 = Transform(registrationName='Transform1', Input=ahAxmf)
        transform1.Transform = 'Transform'
        # Properties modified on transform1.Transform
        transform1.Transform.Translate = [0.0, 0.0, 2.5]
        # show data in view
        transform1Display = Show(transform1, renderView1,
                                 'UnstructuredGridRepresentation')
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        transform1Display.ScaleTransferFunction.Points = [
            1.9988708955637904, 0.0, 0.5, 0.0, 2.000606718406889, 1.0, 0.5, 0.0
        ]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        transform1Display.OpacityTransferFunction.Points = [
            1.9988708955637904, 0.0, 0.5, 0.0, 2.000606718406889, 1.0, 0.5, 0.0
        ]
        # hide data in view
        Hide(ahAxmf, renderView1)
        # show color bar/color legend
        transform1Display.SetScalarBarVisibility(renderView1, True)
        # update the view to ensure updated data information
        renderView1.Update()
        # turn off scalar coloring
        ColorBy(transform1Display, None)
        # Hide the scalar bar for this color map if no visible data.
        HideScalarBarIfNotNeeded(ricciScalarLUT, renderView1)
        # change solid color
        transform1Display.AmbientColor = [0.0, 0.0, 0.0]
        transform1Display.DiffuseColor = [0.0, 0.0, 0.0]
        # toggle 3D widget visibility (only when running from the GUI)
        Hide3DWidgets(proxy=transform1.Transform)
    if (ahb_xmf != ''):
        # create a new 'XDMF Reader'
        ahBxmf = XDMFReader(registrationName=ahb_xmf,
                            FileNames=[file_destination + '/' + ahb_xmf])
        ahBxmf.PointArrayStatus = ['RicciScalar']
        # show data in view
        ahBxmfDisplay = Show(ahBxmf, renderView1,
                             'UnstructuredGridRepresentation')
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        ahBxmfDisplay.ScaleTransferFunction.Points = [
            1.9988745472336247, 0.0, 0.5, 0.0, 2.000593573487867, 1.0, 0.5, 0.0
        ]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        ahBxmfDisplay.OpacityTransferFunction.Points = [
            1.9988745472336247, 0.0, 0.5, 0.0, 2.000593573487867, 1.0, 0.5, 0.0
        ]
        # show color bar/color legend
        ahBxmfDisplay.SetScalarBarVisibility(renderView1, True)
        # update the view to ensure updated data information
        renderView1.Update()
        # create a new 'Transform'
        transform2 = Transform(registrationName='Transform2', Input=ahBxmf)
        transform2.Transform = 'Transform'
        # Properties modified on transform2.Transform
        transform2.Transform.Translate = [0.0, 0.0, 2.5]
        # show data in view
        transform2Display = Show(transform2, renderView1,
                                 'UnstructuredGridRepresentation')
        # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        transform2Display.ScaleTransferFunction.Points = [
            1.9988745472336247, 0.0, 0.5, 0.0, 2.000593573487867, 1.0, 0.5, 0.0
        ]
        # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        transform2Display.OpacityTransferFunction.Points = [
            1.9988745472336247, 0.0, 0.5, 0.0, 2.000593573487867, 1.0, 0.5, 0.0
        ]
        # hide data in view
        Hide(ahBxmf, renderView1)
        # show color bar/color legend
        transform2Display.SetScalarBarVisibility(renderView1, True)
        # update the view to ensure updated data information
        renderView1.Update()
        # turn off scalar coloring
        ColorBy(transform2Display, None)
        # Hide the scalar bar for this color map if no visible data.
        HideScalarBarIfNotNeeded(ricciScalarLUT, renderView1)
        # change solid color
        transform2Display.AmbientColor = [0.0, 0.0, 0.0]
        transform2Display.DiffuseColor = [0.0, 0.0, 0.0]
        # toggle 3D widget visibility (only when running from the GUI)
        Hide3DWidgets(proxy=transform2.Transform)
        # set active source
        SetActiveSource(transform1)

    LoadPalette(paletteName='GradientBackground')
    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0
    # set active source
    SetActiveSource(warpByScalar1)
    # hide color bar/color legend
    warpByScalar1Display.SetScalarBarVisibility(renderView1, False)
    # Properties modified on warpByScalar1Display
    warpByScalar1Display.Opacity = 0.8
    # get layout
    layout1 = GetLayout()
    # layout/tab size in pixels
    layout1.SetSize(1920, 1080)

    #-----------------------------------
    # saving camera placements for views
    # Top down view
    if (camera_angle == 1):
        renderView1.CameraPosition = [0.0, 0.0, 36.90869716569761]
        renderView1.CameraFocalPoint = [0.0, 0.0, 0.6894899550131899]
        renderView1.CameraViewUp = [0, 1, 0]
        renderView1.CameraParallelScale = 424.27024700303446
    # Wide/Inbetween View
    elif (camera_angle == 2):
        renderView1.CameraPosition = [
            -30.093571816984163, -18.55482256667294, 8.558827016411096
        ]
        renderView1.CameraFocalPoint = [
            6.698318283229103e-17, 4.073218385661855e-17, 0.6894899550131903
        ]
        renderView1.CameraViewUp = [0.0, 0.0, 1.0]
        renderView1.CameraParallelScale = 424.27024700303446
    # Side View
    else:
        # current camera placement for renderView1
        renderView1.CameraPosition = [
            -29.944619336722987, -3.666072157343372, 2.895224044348878
        ]
        renderView1.CameraFocalPoint = [
            -0.13267040638072278, 0.6356115665206243, -0.37352608789235847
        ]
        renderView1.CameraViewUp = [0.0, 0.0, 1.0]
        renderView1.CameraParallelScale = 519.6152422706632

    # save animation
    SaveAnimation(output_destination + '/Binary_Pic.png',
                  renderView1,
                  ImageResolution=[1920, 1080],
                  FrameRate=60)


@click.command(help=generate_video.__doc__)
@click.option('--file-destination',
              '-f',
              required=True,
              help="Path to directory containing the volume files xmf and"
              "h5 files")
@click.option('--output-destination',
              '-o',
              required=True,
              help="Output directory where pictures will be placed")
@click.option('--volume-xmf',
              '-v',
              required=True,
              help="Xmf file for VolumeData")
@click.option(
    '--aha-xmf',
    '-a',
    help="Optional xmf file for AhA for apparent horizon visualization")
@click.option(
    '--ahb-xmf',
    '-b',
    help="Optional xmf file for AhB for apparent horizon visualization")
@click.option(
    '--camera-angle',
    default=0,
    type=int,
    help="Determines which camera angle to use: Default 0 is a side view,"
    "1 is a top down, 2 is further out but inbetween side and top down view")
def generate_video_command(**kwargs):
    _rich_traceback_guard = True
    generate_video(**kwargs)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    generate_video_command(help_option_names=["-h", "--help"])
