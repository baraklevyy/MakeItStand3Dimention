# Make It Stand

In this program, we implemented a methodology to modify 3D objects in order to balance them
on defined balance points. The methodology is described in detail in the paper:
_Make it stand: balancing shapes for 3D fabrication_ [[1]](#1).


## Application
Please refer the link for more detailed information.
<a href="Barak_Make_it_Stand.pdf">This is a link</a>
 

In reality when implementing 3D the oblect may not be balance as we desire. My program implement a method that helps users create a balanced 3D object when the only argument necessary are: 1. The point of equallibirium. 2. The orientation of the object.
My program adjust the center of mass of the object and a consequence balance the entire object. 

First we have to choose the point of equilibrium using the menu bar of the program. The desired point of equilibrium will first appear red and when the program ready to balance it'll turn into green.


<img src="pointOfEqulibrium.gif" alt="drawing" width="600">
<br/>
<br/>

Later we want tp adjust the object position. The program let the user to user three rotational angels: Yaw, Roll and Pitch.
Point worth mentioning is the plane shown in the program, this plane represents the ground.

<img src="anglesMoving.gif" alt="drawing" width="600">
<br/>
<br/>


Now that we defined how our object will look statically we can begin the carving process
to make a balanced object. The program let the user to inspect at individual carving step (optional using the **carving for one iteration** button). 

<img src="oneIterationGif.gif" alt="drawing" width="600">
<br/>
<br/>

Beyond single step at a time, the program can carve the whole object until reaching the balanced critiria. 

<img src="carveGif.gif" alt="drawing" width="600">
<br/>
<br/>


## Method From The Bird's Eye View 
In order to balance the object we have to discretize it's shape.
The process to discretize the object known in the jargon as voxelization.
The classification of the object as voxels simplifies the task of object balancing by voxels carving.
Point worth mentioning is that within pre-defined distance from the border of the object we stop carving in order to keep the object of one piece. 

After the voxelization process we want to carve some voxels in order to adject the COM (center of mass).
In a way the COM is exactly above the POB (point of balance) - the object is stable and can stand. For the sake of voxel carving we drawing a plane that is parallel to the gravity direction, it's root is the POB and vertical to COM->POB vector. (Refer to the link above)

As one can infer, one side of the mentioned plane is contained the COM of the object while at the same time the other side is not.
The method is to erase voxels from the contained side, that is, removing voxels from that side potentially make the COM closer to the POB and therfore make the whole object steady.

#### Reference
https://doi.org/10.1145/2461912.2461957
