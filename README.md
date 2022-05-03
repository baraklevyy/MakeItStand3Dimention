# Make It Stand

In this program, we implemented a methodology to modify 3D objects in order to balance them
on defined balance points. The methodology is described in detail in the paper:
_Make it stand: balancing shapes for 3D fabrication_ [[1]](#1).


## Application
Please refer the link for detailed information.
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

As well as to single step at a time, the program can carve the whole object until reaching the balanced critiria. 

<img src="carveGif.gif" alt="drawing" width="600">
<br/>
<br/>


## Implementation of Balancing
The first step to balancing an object is to discretize the object into voxels — commonly known as
voxelization. Converting the object into voxels simplies the problem so that we can carve 
an object by removing a voxel. An important constraint to follow is that voxels within a certain
distance from the border of the mesh should not be removed. This maintains the object's 
integrity when printed.

Next, we want to carve the object to adjust its center of mass. Having a center of mass
directly above the point-of-balance would mean the object is balanced. To carve the object, we 
draw a plane on the point-of-balance. This plane is parallel to the direction of gravity,
and perpendicular to the point-of-balance to center-of-mass vector. An illustration is shown in the paper.

One side of this plane contains the object's current center of mass. This is the side we want to
remove voxels from because removing a voxel will potentially bring the center of mass closer to the plane 
and balance point, and make the object more balanced. Removing a voxel from the other 
side of the plane will always make the shape more unbalanced.

## Code Summary
`main.cpp`
* Displays mesh objects and balancing options
* Controls the workflow of balance properties selection

`inner_void_mesh.h` and `inner_void_mesh.cpp`
* `InnerVoidMesh` class holds the current state of carving; e.g. holds a list of voxels, whether they
are removed, etc.
* `Voxel` holds the definition for a voxel
* Voxalizes the mesh
* Decides which voxel to carve

`list3d.h` and `list3d.cpp`
* 3D vector wrapper. Used to hold the 3D list of voxels

`center_of_mass.h` and `center_of_mass.cpp`
* Computes the center of mass of an object

`grid_util.h` and `grid_util.cpp`
* Helper functions to help transform (translate, scale, rotate) and voxelize objects

## Next Step
This implementation of 3D balancing could be improved by introducing deformation 
to balance the object, allowing the use of points of balance, and using polygons
as bases of balance. These were features described in the paper and were left out
in consideration of time.

## Dependencies

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd libigl-example-project/
    git clone https://github.com/libigl/libigl.git

Make sure to use the `--recursive` flag if you're cloning this project. This will clone the submodules for this project — specifically, `LIBIGL`.


## References
<a id="1">[1]</a> 
Romain Prévost, Emily Whiting, Sylvain Lefebvre, and Olga Sorkine-Hornung. 2013. 
Make it stand: balancing shapes for 3D fabrication. 
ACM Trans. Graph. 32, 4, Article 81 (July 2013), 10 pages. 
DOI: https://doi.org/10.1145/2461912.2461957
