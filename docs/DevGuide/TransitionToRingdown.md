\cond NEVER
Distributed under the MIT License.
See LICENSE.txt for details.
\endcond

# Transition to Ringdown {#dev_guide_transition_to_ringdown}

\tableofcontents

### Introduction

For simulating binary black holes (BBH) in SpECTRE, we have 2 main simulations.
The binary evolution and the single black hole evolution. The transition to
ringdown as the name intends, transitions our BBH simulations from the binary
domain with 2 excisions to the ringdown domain with one excision. This requires
careful placement and sizing/shaping of the excision in the ringdown domain.
This can only be done once we have found the common horizon multiple times in
the binary evolution. The way we place, shape, and size the excision is through
functions of time maps for rotation, expansion, translation, and shape and the
inner radius supplied to the domain.

\warning The transition to ringdown is still experimental, it has only been
tested on the following cases. Equal mass (q = 1) non-spinning, q = 2
non-spinning, and q = 1 non-spinning headon-collision.

In this Dev-Guide we'll walk through exactly what this script is doing.

### Setting the Functions of Time

There are 4 functions of time that need to be set in order to start the
ringdown. The expansion and rotation functions of time are the same as the ones
used in the inspiral, with the only difference being that they settle to
constant values during the ringdown. The translation and shape functions of time
need to be constructed using multiple common horizon finds. This is because we
don't have a control system updating functions of time that track the common
horizon's shape and position during the binary evolution.

To construct the translation and shape functions of time we start by
creating a sphere domain that has the settle to constant versions of the
expansion and rotation map, as well as a FromVolumeFile version of the
translation map used in the inspiral. The FromVolumeFile version reads in the
history of the translation map functions of time used at various times in the
inspiral. The inner radius of this sphere domain is set to 0.01 so that every
point on the common horizon strahlkorper we transform to this domain will be
able to be mapped to a block. We then take around 10 common horizon
strahlkorpers from the inspiral inertial frame at different times and transform
them to the ringdown distorted frame. After the strahlkorper has been
transformed it is recentered so that the expansion center is where the physical
center is. We then save the shape coefficients of the strahlkorper in the
distorted frame. These will be used to construct the shape function of time used
in the ringdown. We then take the expansion center point of the strahlkorper in
the ringdown distorted frame and transform it back to the inspiral inertial
frame. These center points are saved and will be used to construct the
translation function of time for the ringdown. Once we have the shape
coefficients and center points of the strahlkorper at multiple times, we can
then fit these coefficients and center points to a cubic polynomial. Once we
have the cubic polynomial, we can easily get the first and second time
derivatives of the center positions and shape coefficients. These functions and
their first two time derivatives are then used to initialize the translation and
shape maps at the match time of the ringdown.

\note We need somewhere around 10 common horizon finds at different times to get
good fit for the translation and shape functions of time. We also set the 00
coefficient of shape to be 0 since we control the size of the excision with the
inner radius option in the ringdown yamls.

### Choosing an Excision Radius

Once we have the functions of time that will be used in the ringdown, we can now
start the process of choosing a good excision radius. The way we choose an
excision radius is by taking the strahlkorpers that represent the excisions A/B
from the inspiral and transform this to a sphere domain with the correct
functions of time. We start with a sphere domain inner radius of 0.94 * the
average radius of the common horizon strahlkorper at the match time. We then
transform the excisions from the inspiral inertial frame to the ringdown grid
frame. We then loop over every point to ensure that it's inside the ringdown
excision. If a point or points is not inside the ringdown excision, then we set
the minimum excision radius to the point farthest from the center of the common
horizon so that all of the points fit. We then try an excision radius that is
3/4 of the way between the maximum value that should fit the excisions and the
minimum value that should fit the excisions. If all the points fit, then we
increase the L_max on the excision strahlkorpers to make sure multiple
resolutions fit inside this ringdown excision. If all the points fit again then
we have found an excision radius to use in the ringdown. This is all done using
2 main loops, the outer loop, and the inner loop. These are described in more
detail below.

#### Outer Loop

The outer loop controls the L_Max of the excisions and rescales the shape
coefficients by the current inner radius / average radius of common horizon. It
starts with an L_Max of 20 and increases the L_Max by 6 for each iteration.
After the L_Max is set and the strahlkorpers representing the excisions from the
inspiral inertial frame are made, it moves to the inner loop. If the inner loop
converges, and it is not the first iteration of the outer loop, then the L_Max
is increased by 6 and it goes through the inner loop again. If the absolute
value of the difference between the excision radius chosen in the last iteration
and the current iteration are within 2e-4 / q where q is the mass ratio of the
binary then the outer loop  has converged and we have found a suitable excision
radius.

#### Inner Loop

The inner loop controls the radius of the excision used in the ringdown. It
starts by creating a sphere domain with an inner radius of 0.94 * the average
radius of the common horizon at the match time and the settle to const functions
of time for expansion and rotation and the translation and rescaled shape
functions of time constructed for the ringdown. We then transform the excisions
from the inspiral inertial frame to the ringdown grid frame. We then loop over
every point checking that it's inside the ringdown excision. If a point or
points is not inside the ringdown excision, then we set the minimum excision
radius to the point farthest from the center of the common horizon so that all
of the points fit. We then try an excision radius that is 3/4 of the way between
the average radius of the common horizon at the match time and the minimum
excision radius that should fit the excisions. This loop iterates multiple times
until the difference between the last excision radius chosen and the current
excision radius chosen is less than 2e-4 / q where q is the mass ratio of the
binary. When this loop converges, the excision radius is sent to the outer loop
and then the outer loop checks if it has converged.

Not sure if I'll use this.
The reason we can't use the translation function of time from the inspiral is
because it is tracking the center of the apparent horizons A and B not the
center of the common horizon strahlkorper.
