#!/usr/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import os

import click
import rich.columns


def ah_vis(ah_xmf, render_view):
    import paraview.simple as pv

    Ah_xmf = pv.XDMFReader(registrationName=ah_xmf, FileNames=[ah_xmf])
    Ah_xmf.PointArrayStatus = ["RicciScalar"]
    ricciScalarLUT = pv.GetColorTransferFunction("RicciScalar")
    ricciScalarPWF = pv.GetOpacityTransferFunction("RicciScalar")
    Ah_xmf_Display = pv.Show(
        Ah_xmf, render_view, "UnstructuredGridRepresentation"
    )
    Ah_xmf_Display.SetScalarBarVisibility(render_view, False)
    render_view.Update()
    transform1 = pv.Transform(registrationName="Transform1", Input=Ah_xmf)
    transform1.Transform = "Transform"
    transform1.Transform.Translate = [0.0, 0.0, 2.5]
    transform1_display = pv.Show(
        transform1, render_view, "UnstructuredGridRepresentation"
    )
    pv.Hide(Ah_xmf, render_view)
    transform1_display.SetScalarBarVisibility(render_view, False)
    render_view.Update()
    pv.ColorBy(transform1_display, None)
    pv.HideScalarBarIfNotNeeded(ricciScalarLUT, render_view)
    transform1_display.AmbientColor = [0.0, 0.0, 0.0]
    transform1_display.DiffuseColor = [0.0, 0.0, 0.0]
    pv.Hide3DWidgets(proxy=transform1.Transform)


def render_bbh(
    output_dir,
    volume_xmf,
    aha_xmf,
    ahb_xmf,
    camera_angle,
    color_map,
    show_grid,
    stride,
):
    """Generate Pictures from XMF files for BBH Visualizations

    Generate pictures from BBH runs using the XMF files generated using
    generate-xdmf. For files to be read properly, the XMF files being pointed to
    should be in the same directory as the VolumeData.h5 and SurfacesData.h5
    files.

    To splice all the pictures into a video, try using FFmpeg"""
    import paraview.simple as pv

    version = pv.GetParaViewVersion()
    if version < (5, 11) or version > (5, 11):
        print(
            "WARNING: Your Paraview version is not 5.11, "
            "the script may not work correctly."
        )

    # Volume Data Visualization
    volume_files_xmf = pv.XDMFReader(
        registrationName=volume_xmf, FileNames=[volume_xmf]
    )

    # Check for Lapse and SpatialRicciScalar
    variables = volume_files_xmf.PointData.keys()
    lapse_flag = False
    ricci_scalar_flag = False
    for x in variables:
        if x == "Lapse":
            lapse_flag = True
        if x == "SpatialRicciScalar":
            ricci_scalar_flag = True
    if not (lapse_flag and ricci_scalar_flag):
        print(
            "The volume data doesn't have the Lapse or SpatialRicciScalar "
            "output, the script will not work correctly without those."
        )
        return

    render_view = pv.GetActiveViewOrCreate("RenderView")

    # Color the grid
    color_transfer_function = pv.GetColorTransferFunction("Lapse")
    opacity_function = pv.GetOpacityTransferFunction("Lapse")
    color_transfer_function.Discretize = 0
    color_transfer_function.ApplyPreset(color_map, True)
    color_transfer_function.InvertTransferFunction()

    # Slice volume data
    slice = pv.Slice(registrationName="slice", Input=volume_files_xmf)
    slice.SliceType = "Plane"
    slice.HyperTreeGridSlicer = "Plane"
    slice.SliceOffsetValues = [0.0]
    slice.SliceType.Normal = [0.0, 0.0, 1.0]
    slice.Triangulatetheslice = 0

    # Warp grid by spatial ricci scalar
    warp_by_scalar = pv.WarpByScalar(
        registrationName="WarpByScalar", Input=slice
    )
    warp_by_scalar.Scalars = ["POINTS", "SpatialRicciScalar"]
    warp_by_scalar.ScaleFactor = 2.5
    warp_by_scalar.Normal = [0.0, 0.0, -1.0]
    warp_by_scalar_display = pv.Show(
        warp_by_scalar, render_view, "GeometryRepresentation"
    )
    warp_by_scalar_display.SetScalarBarVisibility(render_view, False)

    # Apparent Horizon Visualization
    if aha_xmf:
        ah_vis(aha_xmf, render_view)
    if ahb_xmf:
        ah_vis(ahb_xmf, render_view)

    if show_grid:
        warp_by_scalar_display.Representation = "Surface With Edges"

    pv.LoadPalette(paletteName="GradientBackground")
    render_view.OrientationAxesVisibility = 0
    pv.SetActiveSource(warp_by_scalar)
    warp_by_scalar_display.Opacity = 0.8
    pv.ColorBy(warp_by_scalar_display, ("POINTS", "Lapse"))
    layout = pv.GetLayout()
    layout.SetSize(1920, 1080)

    # Camera placements
    # Top down view
    if camera_angle == "1":
        render_view.CameraPosition = [0.0, 0.0, 36.90869716569761]
        render_view.CameraFocalPoint = [0.0, 0.0, 0.6894899550131899]
        render_view.CameraViewUp = [0, 1, 0]
        render_view.CameraParallelScale = 424.27024700303446
    # Wide/Inbetween View
    elif camera_angle == "2":
        render_view.CameraPosition = [
            -89.0,
            -17.0,
            25.0,
        ]
        render_view.CameraFocalPoint = [
            -0.3921962951264054,
            1.6346750682876983,
            -0.34522248814953405,
        ]
        render_view.CameraViewUp = [
            0.0,
            0.0,
            1.0,
        ]
    # Side View
    else:
        render_view.CameraPosition = [
            -29.944619336722987,
            -3.666072157343372,
            2.895224044348878,
        ]
        render_view.CameraFocalPoint = [
            -0.13267040638072278,
            0.6356115665206243,
            -0.37352608789235847,
        ]
        render_view.CameraViewUp = [0.0, 0.0, 1.0]
        render_view.CameraParallelScale = 519.6152422706632

    # save animation
    animation_scene = pv.GetAnimationScene()
    animation_scene.Stride = stride
    pv.SaveAnimation(
        os.path.join(output_dir, "Binary_Pic.png"),
        render_view,
        ImageResolution=[1920, 1080],
        FrameRate=60,
    )


@click.command(help=render_bbh.__doc__)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Output directory where pictures will be placed",
)
@click.option(
    "--volume-xmf",
    "-v",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    required=True,
    help="Xmf file for VolumeData",
)
@click.option(
    "--aha-xmf",
    "-a",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Optional xmf file for AhA visualization",
)
@click.option(
    "--ahb-xmf",
    "-b",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Optional xmf file for AhB visualization",
)
@click.option(
    "--camera-angle",
    "-c",
    default="0",
    type=click.Choice(["0", "1", "2"]),
    help=(
        "Determines which camera angle to use: Default 0 is a side view,"
        "1 is a top down, 2 is further out but inbetween side and top down view"
    ),
)
@click.option(
    "--color-map",
    "-m",
    default="Rainbow Uniform",
    help=(
        'Determines how to color the domain, common color maps are "Inferno'
        ' (matplotlib)", "Viridis (matplotlib)"'
    ),
)
@click.option(
    "--show-grid",
    is_flag=True,
    help="Show grid lines",
)
@click.option(
    "--stride",
    default=1,
    help=(
        "Stride determines how many timesteps to move forward before taking"
        " another picture"
    ),
)
def render_bbh_command(**kwargs):
    _rich_traceback_guard = True
    render_bbh(**kwargs)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    render_bbh_command(help_option_names=["-h", "--help"])
