 \cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond
# Running a Single Black Hole and Visualizing it In Paraview

Welcome to the Running and Visualizing a Single Black Hole Tutorial. This
tutorial is designed to make you more comfortable with running a single black
hole and give you some of the key concepts of visualizing the data produced in
SpECTRE simulations. The concepts used for visualization can be applied to
binary simulations as well!

## Prerequisits

__Important!__
Before doing any of the commands below, make sure you have built
EvolveGhSingleBlackHole within your build directory. If you're unfamiliar with
this process, cd into your build directory, source and load modules and then run
the command
```
make -j4 EvolveGhSingleBlackHole
```

## Running a Single Black Hole

To submit and run a job, we'll need 2 things; the submit script for the cluster
you're on and the KerrSchild.yaml file. To find the right submit script look in
__spectre/support/SubmitScripts/__ and for the yaml file you'll want to look in
__spectre/tests/InputFiles/GeneralizedHarmonic/__. I'm on the Ocean cluster so
my submit script is __OceanClang.sh__. You'll want to copy both of these files
into a new directory.

We'll need to make a few edits to these files before they're ready to run. Let's
start with the KerrSchild.yaml file. Start by changing the OuterRadius in the
DomainCreator to 50.0
```
DomainCreator:
  Shell:
    InnerRadius: 1.9
    OuterRadius: 50.0
```
It should look something like this. We want to change the outer radius so that
our visualization doesn't look like a tiny sad donut!

Next, we'll want to look at the EventsAndTriggers section. Here, find the
VolumeData subfile option, we want to edit the Interval above to be 100. SpECTRE
simulations can output 100's of Gb's of data if you're not careful. Changing the
interval allows us to control how much data we're writing to disk (This
simulation will output somewhere between 25-40Mb). Then, within
VariablesToObserve, we'll want to add __- Lapse__ and __- SpatialRicciScalar__.
Once you've done this, this section should look something like this.
```
  - - Slabs:
        EvenlySpaced:
          Interval: 100
          Offset: 0
    - - ObserveFields:
          SubfileName: VolumeData
          VariablesToObserve:
            - SpacetimeMetric
            - Pi
            - Phi
            - GaugeH
            - PointwiseL2Norm(GaugeConstraint)
            - PointwiseL2Norm(ThreeIndexConstraint)
            - PointwiseL2Norm(FourIndexConstraint)
            - Lapse
            - SpatialRicciScalar
          InterpolateToMesh: None
          CoordinatesFloatingPointType: Double
          FloatingPointTypes: [Double]
```
The last thing we need to do in EventsAndTriggers is replace these lines above
completion
```
  - - Slabs:
        Specified:
          Values: [3]
```
with these lines
```
  - - TimeCompares:
        Comparison: GreaterThanOrEqualTo
        Value: 10.0
```
This means our simulation will evolve the black hole to time 10M. This will give
us 11 'frames' of data in Paraview. The way we know this is by looking at the
Interval, InitialTimeStep and TimeCompares value options. At the top, our
InitialTimeStep is 0.01 and since we set the Interval for VolumeData is 100,
that means every 100 steps it'll write the data to disk. This means we'll get a
frame at time 0, 1, ..., 10.

Now, let's open the submit script. We'll start by editing some of the options at
the top of the file. The -J option controls the what job name will appear in the
queue, the --nodes option specifies how many nodes to use, we'll change this
option to 1 and finally, -t specifies how much time to allocate to this job,
this means that even if the job doesn't finish by this time, it will be
terminated by the server, this can be left alone for now. Next, we'll need to
change and specify 4-5 commands depending on what cluster you're on.
```
export SPECTRE_BUILD_DIR=${HOME}/spectre/build
export SPECTRE_MODULE_DIR=${HOME}/spectre_deps/modules/
export SPECTRE_RUN_DIR=${PWD}/Run
export SPECTRE_EXECUTABLE=${SPECTRE_BUILD_DIR}/bin/EvolveGhSingleBlackHole
export SPECTRE_INPUT_FILE=${PWD}/KerrSchild.yaml
```
If you're working on a cluster that doesn't require spectre a spectre
dependecies directory, then you should safely comment out or remove the
SPECTRE_MODULE_DIR line altogether.

With this, you're ready to submit a job! Now, while in the directory containing
both files, all you need to do is run this command
```
sbatch OceanClang.sh
```
Where you change out OceanClang.sh for the submit script you copied.

Your job should be running now! A way to check this is by looking at the
spectre.stdout and spectre.stderr files. Once your job is complete (this should
be rather quick) you should run these two commands
1. spectre generate-xdmf GhKerrSchildVolume*.h5 --output VolumeFiles
--subfile-name VolumeData
2. spectre generate-xdmf GhKerrSchildSurfaces.h5 --output AhA --subfile-name AhA
The first command specifies the files needed for visualizing the Volume data
that's been written to the GhKerrSchildVolume.h5 files and the second command
allows us to visualize the apparent horizon of the black hole!
(If you're visualizing a binary, and want the other black hole, run the second
command again, replacing AhA with AhB)

## Paraview Time
Now that we have generated the xmf files needed, you can either download the run
directory to your local machine using rsync (for mac or linux) or scp (windows)
or you can visualize the data using your cluster and connecting your local
Paraview client to the cluster using pvserver. Since this run was short and the
data is small, I'll be downloading it to my local.

Once you've downloaded the Run directory, we can open Paraview 5.10.1. Once it's
open, find the file icon and click on it or use the file bar and click on open.

\image html Paraview_Tutorial1.png

Now navigate to where you downloaded the Run directory and open the
VolumeFiles.xmf file. This will prompt you to choose an XDMF reader, choose the
"XDMF Reader" option, clicking on any other option will crash Paraview :)

Once you've opened the xmf file and click apply, you should get a picture of a
gumball!

\image html Paraview_Tutorial2.png

This gumball is the computational domain where all the action occurs. However,
we can't see any of the interesting stuff going on yet. For that, we'll need to
take a slice through the z-axis. This can be done by clicking the slice option
and specifying the Z normal.

\image html Paraview_Tutorial3.png

Once you click apply, you should get a picture like this. The empty section in
the middle is where our black hole lives! That region is called the excision
surface.

\image html Paraview_Tutorial4.png

Now that we have a slice, we can use some of the filters to visualize the
quantities we output. The first filter we're going to use is Warp by Scalar, you
can find it by navigating the filters tab and looking through alphabetical.
Once you've found and clicked on it, you'll see soemething like this.

\image html Paraview_Tutorial5.png

Now change the Scalars to SpatialRicciScalar and I encourage you to mess around
with the scale factor to see what happens and find something you like! Once
you've done that, you should get a picture like this.

\image html Paraview_Tutorial6.png

It's time to visualize the black hole! For this we'll need to open the other
xmf file we generated. Open the file AhA.xmf with the same "XDMF Reader" option
and click apply. This should put a small gumball resting in the well we've
created with the Warp by Scalar option. This small gumball is the apparent
horizon of the black hole!

\image html Paraview_Tutorial7.png

However, the apparent horizon is now shielding our eyes from seeing the
curvature it's imposing on our computational domain. To remedy this, we'll need
to move the apparent horizon up a bit. To do this, we'll need the Transform
filter. Once you find it, you should see a screen like this.

\image html Paraview_Tutorial8.png

I once again recommend you to mess around with the Translate options to see what
happens! For me, I usually settle on moving it -3 in the z direction so that it
floats a bit above the Warp by Scalar. Once you've found a translation you like,
you should see something like this.

\image html Paraview_Tutorial9.png

Finally, let's learn about how to change the colors in Paraview. Let's first
change the color of the Apparent Horizon. If you scroll down in the properties
where you changed the translation, you'll find a section for colors. Change the
coloring from RicciScalar to Solid color. Now click edit, you can choose to make
it any color you like. I'll be going with solid black for no particular reason.
Next, let's change the colors of the slice. To do this click on the
WarpByScalar1 section in the Pipeline Browser and scroll down through the
properties until you find coloring again. Change the coloring to Lapse and then
click edit. This should open the Color Map Editor on the right side of Paraview.
Right under mapping data you can find all the dfiferent color maps. A good
practice when using color maps is always using unifrom color maps since we're
representing data with color here. Some of the uniform color maps are Inferno,
Viridis and Uniform Rainbow. I also encourage you to mess around with the
buttons to the right of the color range because some of them are very useful.
Another very useful option in the Color Map Editor is the
__Use Log Scale When Mapping Data To Colors__ option. This can be very helping
when trying to see tiny differences. Another option to make the colors smoother
and nicer to look at is the Color Discretization option. Unchecking this option
smooths out the color mapping. Once you're satisfied with all the options you
should get a picture like this!

\image html Paraview_Tutorial10.png

Congratulations, you've just run a black hole simulation and visualized it using
Paraview! Now, all that's left to do is click play and see what happens! For a
single black hole, nothing too exciting should happen. However, now you should
have a good grasp of core concepts and ideas for visualizing binaries!
