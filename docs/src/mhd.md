# Magnetohydrodynamics Model

!!! note
    As of June 2021, KaTeX lacks the full support for equation numbering.

## Ideal MHD


The ideal MHD equations can be written in (near) conservative form as
```math
\begin{aligned}
\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho \mathbf{u}) &= 0, \\
\frac{\partial \rho \mathbf{u}}{\partial t} + \nabla\cdot \left[ \rho \mathbf{u} \mathbf{u} + \bar{\bar{I}} (p + \frac{1}{2} B^2)
- \mathbf{B}\mathbf{B} \right] &= -\mathbf{B}\nabla\cdot\mathbf{B}, \\
\frac{\partial \mathbf{B}}{\partial t} + \nabla\cdot (\mathbf{u}\mathbf{B} - \mathbf{B}\mathbf{u})
&= -\mathbf{u}\nabla\cdot\mathbf{B}, \\
\frac{\partial e}{\partial t} + \nabla\cdot\left[ \mathbf{u} ( e + p + \frac{1}{2} B^2 )
-\mathbf{u}\cdot\mathbf{B}\mathbf{B} \right] &= -\mathbf{u}\cdot\mathbf{B}\nabla\cdot\mathbf{B},
\end{aligned}
```
where ``\rho`` is mass density, ``\mathbf{u}`` is velocity, ``p`` is pressure,
``\mathbf{B}`` is the magnetic field and ``\bar{\bar{I}}`` is the identity matrix. 
The total energy density is
```math
  \tag{totalE}
  e = \frac{p}{\gamma-1} + \frac{\rho u^2}{2} + \frac{B^2}{2},
```
where ``\gamma`` is the adiabatic index.

!!! note
    We keep the divergence of the magnetic field on the RHS because of numerical accuracy. Even though physically valid by far in nature, it is not guaranteed to be zero in a numerical model.

We can also substitute the energy equation with pressure equation for a non-conservative form:
```math
  \tag{pressure}
  \frac{\partial p}{\partial t} + \nabla\cdot (p \mathbf{u}) = -(\gamma - 1) p \nabla \cdot \mathbf{u}, 
```
The source terms proportional to ``\nabla\cdot\mathbf{B}`` and ``\nabla\cdot\mathbf{u}`` on the right hand sides are not evaluated as a flux, but the ingredients (normal components of ``\mathbf{B}`` and ``\mathbf{u}`` are calculated and used later. In total, there are 8 primitive variables ``(\rho,\mathbf{u},\mathbf{b},p)``. Note that in BAT-S-RUS, the momentum ``\rho\mathbf{u}`` is actually stored instead of velocity ``\mathbf{u}``. 

In Cartesian coordinates, we can express the above equations in more detailed and compact form. If we use ``\mathbf{A}`` to represent all the state variables, we have
```math
\frac{\partial\mathbf{A}}{\partial t}+\nabla\cdot\vec{\mathbf{F}}(\mathbf{B})=\mathbf{S},
```
by defining the state and flux vectors as
```math
\begin{aligned}
\text{state: } \mathbf{A} = 
	\begin{bmatrix}
	\rho \\ \rho u_x \\ \rho u_y \\ \rho u_z \\ B_x \\ B_y \\ B_z \\ p \\ e \\
	\end{bmatrix},
\quad
\text{flux: } \vec{\mathbf{F}} = 
	&\begin{bmatrix}
	\rho u_x \\ \rho u_x u_x + (p+\frac{1}{2}B^2)- B_xB_x \\ \rho u_x u_y - B_x B_y \\ \rho u_x u_z - B_xB_z \\ 0 \\ u_x B_y - B_x u_y \\
	u_x B_z - B_x u_z \\ p u_x \\u_x (e+p+\frac{1}{2}B^2) - (u_xB_x+u_yB_y+u_zB_z)B_x
	\end{bmatrix}\widehat{x} + \\
	&\begin{bmatrix}
	\rho u_y \\ \rho u_y u_x - B_yB_x \\ \rho u_y u_y + (p+\frac{1}{2}B^2)- B_y B_y \\ \rho u_y u_z - B_yB_z \\ u_y B_x - B_y u_x \\ 0 \\
	u_y B_z - B_y u_z \\ p u_y \\ u_y (e+p+\frac{1}{2}B^2) - (u_xB_x+u_yB_y+u_zB_z)B_y
	\end{bmatrix}\widehat{y} +\\
	&\begin{bmatrix}
	\rho u_z \\ \rho u_z u_x - B_zB_x \\ \rho u_z u_y - B_z B_y \\ \rho u_z u_z + (p+\frac{1}{2}B^2) - B_zB_z \\ u_z B_x - B_z u_x \\  u_z B_y - B_z u_y \\ 0 \\ p u_z \\u_z (e+p+\frac{1}{2}B^2) - (u_xB_x+u_yB_y+u_zB_z)B_z
	\end{bmatrix}\widehat{z}, \\
\text{source: } \mathbf{S} = 
	&\begin{bmatrix}
	0 \\ 
	-B_x(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	-B_y(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	-B_z(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	-u_x(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	-u_y(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	-u_z(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	-(\gamma-1)p(\frac{\partial u_x}{\partial x}+\frac{\partial u_y}{\partial y}+\frac{\partial u_z}{\partial z}) \\
	-(u_xB_x+u_yB_y+u_zB_z)(\frac{\partial B_x}{\partial x}+\frac{\partial B_y}{\partial y}+\frac{\partial B_z}{\partial z}) \\
	\end{bmatrix}.
\end{aligned}
```
Note that the last two equations overspecify the system. We only need to choose one of them.

Up to this stage, all the equations above are exact. The next question is how do we solve this system numerically?
