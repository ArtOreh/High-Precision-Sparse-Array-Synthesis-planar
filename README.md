# Continuous Sparse Planar Array Synthesis (Work in Progress)

This repository contains the results of an ongoing research on 2D sparse antenna arrays. The source code is currently private until the paper is submitted.

## Key Features of the Results
* **Continuous Coordinates:** Points are not restricted to a grid.
* **Uniform Amplitude:** All array elements have equal power excitation.
* **Circular Aperture:** Layouts are bounded within a unit circle.

## Folder Structure

* `results/` — Contains the output data.
  * `coordinates/` — Text files with precise `[x, y]` coordinates.
  * `plots/` — Array layouts and 2D pattern graphs. Note: Gray circles represent a safety zone of radius $d_{min}/2$ around each element to ensure physical spacing.
