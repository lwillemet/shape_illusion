# shape_illusion
Simulation of the exploration of a physical and a virtual bump with a finite-element model
For a description of the mechanical model, see the project Finger-skin-model available on Github.

The two scripts (PhysicalBump and VirtualBump) solves stresses and strains when the finger is sliding on a surface. At the end of the execution, you will obtain two matrices UPhy and FPhy (resp. UVir and FVir) corresponding respectively to the displacements and the forces acting on each element. The first element corresponds to the bone and the other from the left one to right one of the chain. The length of the matrices is 2*(Ndx+1) where Ndx is the number of particles (default value : 501). The first half of the vector is the normal displacement (resp. force) and the second half, the tangential displacement (resp. force). The number of columns of matrices is equal to the length of the time vector.

![profile](https://user-images.githubusercontent.com/58175648/119663819-0582cd80-be33-11eb-9493-b8b1020d8b9f.png)

# PhysicalBump
The relief follows a gaussian curve with period defined to 4 mm, corresponding approximately to a third of the contact area (red curve on the graph above). 

# VirtualBump
The friction variation profile is sinusoidal with the same period (blue curve on the graph above).

# curspace
This function is called at the beginning of the script finger-skin-model-static and creates an evenly spaced points distribution along the external layer of the finger in 2D. Then, at the beginning of the stimulation each spring has the same length at rest.

# f_U
This function is called in the script finger-skin-model-static to compute the displacement of each particles knowing the previous displacement of these elements and the forces exerted on each. As Runge-Kutta is used at the 4th order, the function is called 4 times.
