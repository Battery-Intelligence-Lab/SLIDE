---
title: 'SLIDE: A C++ software for high-fidelity simulation of lithium-ion battery energy storage system degradation'
tags:
  - lithium-ion batteries
  - degradation
  - battery pack simulation
  - C++
authors:
  - name: Jorn M. Reniers
    orcid: 0000-0002-4186-8696
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Volkan Kumtepeli
    orcid: 0000-0003-2392-9771
    affiliation: 1
  - name: David A. Howey
    orcid: 0000-0002-0620-3955
    affiliation: 1
affiliations:
  - name: Department of Engineering Science, University of Oxford, OX1 3PJ, Oxford, UK
    index: 1
  - name: Brill Power, OX1 1HU, Oxford, UK
    index: 2
date: 18 May 2022
bibliography: paper.bib
---

# Summary

Disclaimer: paper writing is still ongoing; please do not use this version as a reference. 

Day by day, the effects of climate change become more and more obvious and ambitious net-zero goals set by countries drive the deployment of renewable energy generators futher. This, in turn, increases the need for grid balancing, voltage support, and other services that leads to a higher demand for grid-scale energy storage systems. Due to rapid cost declines of lithium-ion batteries, they are increasingly becoming an important part of grid infrastructure, participating in the various markets for
frequency control, system reserves, and wholesale energy trading [@reiners2022digital], [@reniers2019review].

``Slide`` (simulator for lithium-ion degradation) is a comprehensive modelling software for fast simulation of multi-cell lithium-ion battery energy storage system degradation. It aims to provide the efficient tools to implement a high-fidelity lithium-ion battery model, configure and perform simulation runs, and evaluate the resulting data. Code for this software project is mainly written in C++ to do fast simulations of degradation of lithium-ion batteries. Simulating 5000 1C CC cycles should take less than 1 minute; adding a CV phase doubles the calculation time to below 2 minutes. The underlying battery model is the Single Particle Model (SPM) with a coupled bulk thermal model. ``Slide`` adds various degradation models on top of the SPM. The equations were taken from literature and implemented in one large coupled model. Users can easily select which models they want to include in their simulations. They can set the values of the fitting parameters of those degradation models to fit their own data. The project uses object oriented programming in C++, see documentation for more details. 

Other projects for lithium-ion battery simulation and analysis may include: [@kumtepeli2020energy], [@naumann2017simses], [@moller2022simses], and [@tranter2022liionpack].

# Mathematical background

``Slide`` builds upon the spectral single particle model implemented by Adrien Bizeray, Jorn Reniers and David Howey, which is available on [GitHub](https://github.com/davidhowey/Spectral_li-ion_SPM). The reader is recommended to first familiarise with the spectral SPM before turning attention to Slide, which extends the spectral SPM with various degradation mechanisms.

There is one main mathematical difference. The spectral SPM used the *transformed concentration* $u$, which was obtained by multiplying the radius with the lithium concentration at that radius. The final equation of the state space model for solid diffusion of the spectral SPM was:

$$\frac{\partial u_{2:N}}{\partial t} = D_{s,i}(T)Au_{2:N}  + Bj_i(t)$$

Where $u$ is the transformed concentration at the inner nodes, $D_{s,i}(T)$ is the solid diffusion
coefficient of electrode $i$ at temperature $T$, $j_i(t)$ is the current density on 
electrode $i$ at time $t$ and $A$ and $B$ are the state space matrices. 
Although the accuracy of this model is very high even for a low number of discretisation nodes ($N$),
the matrix $A$ is full. This project adds one extra transformation to the eigenspace in order to 
obtain a sparse (diagonal) state space matrix.

Using the eigenvalue decomposition of $A$ we can write:

$$\frac{\partial u_{2:N}}{\partial t} = D_{s,i}(T)V^{-1}\Lambda Vu_{2:N}  + VBj_i(t)$$

Using $z_{2:N}=Vu_{2:N}$ and $\tilde{B}=VB$  we can write:

$$\frac{\partial z_{2:N}}{\partial t} = D_{s,i}(T)V^{-1}\Lambda Vz_{2:N}  + \tilde{B}j_i(t)$$

This equation has the same format as the original equation, but the matrix $\Lambda$ is diagonal, which increases calculation speed and reduces numerical errors.

This transformation has one additional advantage. One of the eigenvectors represents a uniform concentration. This means that if we start with a uniform concentration, only one value of $z$ will be non-zero. 
This also means that the corresponding eigenvalue in $\Lambda$ is 0 (the time derivative of a uniform
concentration must be 0 if there is no external current $j_i(t)$ such that the amount of lithium is conserved.

# Acknowledgements

We gratefully acknowledge the contributions by [Battery Intelligence Lab](https://howey.eng.ox.ac.uk) members. 


# References

