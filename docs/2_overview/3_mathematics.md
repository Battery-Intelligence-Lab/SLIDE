---
layout: default
title: Preliminaries
nav_order: 3
---


# Mathematical background

_Slide_ builds upon the spectral single particle model implemented by Adrien Bizeray, Jorn Reniers and David Howey, which is available on [GitHub](https://github.com/davidhowey/Spectral_li-ion_SPM). The reader is recommended to first familiarise with the spectral SPM before turning attention to Slide, which extends the spectral SPM with various degradation mechanisms.

There is one main mathematical difference. The spectral SPM used the *transformed concentration* $$u$$, which was obtained by multiplying the radius with the lithium concentration at that radius. The final equation of the state space model for solid diffusion of the spectral SPM was:

$$\frac{\partial u_{2:N}}{\partial t} = D_{s,i}(T)Au_{2:N}  + Bj_i(t)$$

Where $$u$$ is the transformed concentration at the inner nodes, $$D_{s,i}(T)$$ is the solid diffusion
coefficient of electrode $$i$$ at temperature $$T$$, $$j_i(t)$$ is the current density on 
electrode $$i$$ at time $$t$$ and $$A$$ and $$B$$ are the state space matrices. 
Although the accuracy of this model is very high even for a low number of discretisation nodes ($$N$$),
the matrix $$A$$ is full. This project adds one extra transformation to the eigenspace in order to 
obtain a sparse (diagonal) state space matrix.

Using the eigenvalue decomposition of $$A$$ we can write:

$$\frac{\partial u_{2:N}}{\partial t} = D_{s,i}(T)V^{-1}\Lambda Vu_{2:N}  + VBj_i(t)$$

Using $$z_{2:N}=Vu_{2:N}$$ and $$\tilde{B}=VB$$  we can write:

$$\frac{\partial z_{2:N}}{\partial t} = D_{s,i}(T)V^{-1}\Lambda Vz_{2:N}  + \tilde{B}j_i(t)$$

This equation has the same format as the original equation, but the matrix $$\Lambda$$ is diagonal, which increases calculation speed and reduces numerical errors.

This transformation has one additional advantage. One of the eigenvectors represents a uniform concentration. This means that if we start with a uniform concentration, only one value of $$z$$ will be non-zero. 
This also means that the corresponding eigenvalue in $$\Lambda$$ is 0 (the time derivative of a uniform
concentration must be 0 if there is no external current $$j_i(t)$$ such that the amount of lithium is conserved.
