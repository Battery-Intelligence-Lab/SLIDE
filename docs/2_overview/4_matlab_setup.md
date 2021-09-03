---
layout: default
title: MATLAB setup
nav_order: 4
---


# MATLAB setup

<p style='text-align: justify;'>
Chebyshev spectral methods are used to discretise the diffusion PDE in space. The result of 
the mathematical derivation are some matrices which are needed to calculate the derivatives, 
actual concentrations etc.
</p>

$$\frac{\partial z_{2:N}}{\partial t} = D_{s,i}(T) \tilde{A}_iz_{2:N}  + \tilde{B}_i j_i(t) $$

$$c_i = \begin{bmatrix} \tilde{C}_iz_{i,2:N} + \frac{\tilde{D}_ij_i(t)}{D_{s,i(T)}} \\ \tilde{C}_{i,c} \left( \tilde{C}_{i}z_{i,2:N} +  \frac{\tilde{D}_ij_i(t)}{D_{s,i(T)}} \right) + \frac{j_i(t)D_{s,i}(T)}{\tilde{D}_i} \end{bmatrix}$$


| Parameter             | Description |
|-----------------------|-------------|
| $$z_{2:N}$$           | Transformed concentration at the (positive) inner nodes in electrode $$i$$                                                            |
| $$D_{s,i}(T)$$        | Solid diffusion constant of electrode $$i$$ at temperature $$T$$                                                                      |
| $$\tilde{A}_i$$       | Diagonal matrix of the state space model for electrode $$i$$ (same matrix as $$\Lambda$$)                                             |
| $$\tilde{B}_i$$       | Column matrix of the state space model for electrode $$i$$                                                                            |
| $$j_i(t)$$            | Current density at time $$t$$ on electrode $$i$$                                                                                      |
| $$c_i$$               |  Concentration at the inner and centre nodes of electrode $$i$$ [mol / m$$^3$$]                                                       |
| $$\tilde{C}_i$$       | Matrix of the state space model linking the actual and transformed concentrations for electrode $$i$$                                 |
| $$\tilde{D}_i$$       | Matrix of the state space model linking the actual concentration and the current density for electrode $$i$$                          |
| $$\tilde{C}_{i,c}$$   | Matrix of the state space model linking the actual concentration at the centre node to the actual concentration of the inner nodes    |


These matrices are calculated by the MATLAB function ```modelSetup.m```, which also writes their values in *.```*.csv``` files. In the C++ code, the values of these matrices are read from these files and stored in a Struct.

However, the values of the matrices depend on 3 parameters that the user can change: the radius of each particle and the number of discretisation nodes. Of course, the same values must be used by the MATLAB code as by the C++ code (else you are spatially discretising for a radius r1, while you are calculating things for a radius r2, which obviously produces wrong results). When the C++ code reads the matrices, it checks that the parameters are identical, and throws an error if they are not. In that case, the user has to change the parameters in the MATLAB code and re-run it to re-calculate the matrices for the new parameters.


In the MATLAB code, these parameters are defined at the top of the function ```modelSetup.m```
- **nch:** the number of inner Chebyshev nodes in the positive domain. The value of this parameter must be the same value as the one used by the C++ code (which is defined on top of ```constants.hpp``` as ```nch```)
- **Rp:** The radius of the positive particle. The value must be the same as the value used by the C++ code, which is defined in the constructors of the child-classes of Cell: ```cell_fit.cpp```, ```cell_KokamNMC.cpp```, ```cell_LGChemNMC.cpp``` and ```cell_user.cpp```
- **Rn:** The radius of the negative particle. The value must be the same as the value used by the C++ code, which is defined in the constructors of the child-classes of Cell: ```cell_fit.cpp```, ```cell_KokamNMC.cpp```, ```cell_LGChemNMC.cpp``` and ```cell_user.cpp```.


The MATLAB code calls some MATLAB functions which implement the mathematical formulas described in the Spectral SPM and the eigenvalue transformation described above. It creates the following ```*.csv``` files:

- **Cheb_input:** contains the number of inner nodes, and the radius. This file is read by the C++ code to verify that MATLAB used the same values of ```nch```, ```Rp``` and ```Rn``` as the C++ code
- **Cheb_nodes:** contains the nondimentional location of the inner nodes of the positive Chebyshev domain. The file is read by the C++ code to fill in the variable ```xch``` in the Model-structure
- **Cheb_Ap:** diagonal elements of the matrix ```Ap``` (the other elements are 0)
- **Cheb_An:** diagonal elements of the matrix ```An``` (the other elements are 0)
- **Cheb_Bp:** elements of the column matrix ```Bp```
- **Cheb_Bn:** elements of the column matrix ```Bn```
- **Cheb_Cp:** elements of the matrix ```Cp```
- **Cheb_Cn:** elements of the matrix ```Cn```
- **Cheb_Cc:** elements of the column matrix to get the concentration at the centre node
- **Cheb_Dp:** elements of the column matrix ```Dp```
- **Cheb_Dn:** elements of the column matrix ```Dn```
- **Cheb_Q:** elements of the Chebyshev integration matrix ```Q```
- **Cheb_Vp:** inverse of the eigenvector of the positive electrode
- **Cheb_Vn:** inverse of the eigenvectors of the negative electrode


As long as the input parameters (```N``` or ```nch```, ```Rp``` and ```Rn```) are not changed, the initial MATLAB setup can be skipped.

The functions used to calculate the updated matrices come from:

- A MATLAB Differentiation Matrix Suite by JAC Weideman and SC Reddy. The details are published in JAC Weideman, SC Reddy, A MATLAB differentiation matrix suite, ACM Transactions of Mathematical Software, Vol 26, pp 465-519 (2000). (DOI 10.1145/365723.365727). (See [1](https://uk.mathworks.com/matlabcentral/fileexchange/29-dmsuite), [2](http://appliedmaths.sun.ac.za/~weideman/research/differ.html))

- [Chebfun V4](https://uk.mathworks.com/matlabcentral/fileexchange/23972-chebfun-v4-old-version-please-download-current-version-instead?focused=5592790&tab=function) by the Chebfun Team. 


```
Copyright (c) 2015, The Chancellor, Masters and Scholars of the University of Oxford, and the Chebfun Developers
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```