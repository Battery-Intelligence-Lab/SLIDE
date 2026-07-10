---
layout: default
title: Advanced
nav_order: 1
---


# Advanced changes to the code 

> The pages in this section describe legacy v3 internals. For v4, use [Add a cell model](../v4/extending-cell-model.html) and [Add an ageing mechanism](../v4/extending-ageing.html). v4 uses compiled batches and scalar-generic kernels, not runtime `Cell_user`/`DEG_ID` registration.
<p style='text-align: justify;'>
This document gives some explanation of how to do advanced changes to the code. It is highly recommended you do not attempt this unless you know how to program in C++ (including object-oriented programming) and you understand how the code in this project works.

All guidance in this document is only guidance, not a procedure which is guaranteed to work. The changed mentioned in this document haven’t been tested so the list of things to change might not be exhaustive (and errors will occur if there are things missing and you only do what is listed here).
</p>

[Adding new cell type](2_new_cell.md)\
[Adding new degradation](3_new_degradation.md)\
[Adding new state](4_new_state.md)
