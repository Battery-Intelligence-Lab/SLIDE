---
layout: default
title: v4 preview
nav_order: 0
---

# SLIDE v4 preview

The v4 API is built around four layers:

1. value-semantic cell, electrode, parameter, experiment, and pack descriptions;
2. cold compilation into validated model compositions and flat topology;
3. allocation-free batch and pack stepping over a structure-of-arrays state arena;
4. symmetric C++, Python, and MATLAB entry points.

Start with [installation](installation.html), then choose a language quickstart. Read [compatibility](compatibility.html) before assuming that a PyBaMM object or optional backend is supported.

The performance contract and detailed migration evidence live in
[`PLAN.md`](https://github.com/Battery-Intelligence-Lab/SLIDE/blob/master/PLAN.md).
The docs focus on supported user and contributor workflows.
