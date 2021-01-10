# KarinMosimann
Image analysis codes used for performing Mercator map projection

This repository contains image analysis codes for performing Mercator map projection in images of zebrafish embryo development obtained from a multiview light-sheet microscope
The reference preprint for this work can be found here: https://www.biorxiv.org/content/10.1101/2020.11.11.355693v2

The repository contains the following image analysis codes:

1. Get_isosurface.m
This code obtains the isosurface from fused images and performs a sphere fit of the point cloud

2. Map_projection_Mercator.m
This code peforms Mercator map projection of fused spherical data

3. Map_projection_pixel_map.m
This code determines the pixel size in the projected map

4. Compare_embryo_fit_sphere_over_time.m
This code compares the sphere fit radius over time
