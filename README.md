# MidSurfer

![Teaser](images/teaser.png)

## Overview

MidSurfer is a novel parameter-free method for extracting mid-surfaces from segmented volumetric data. Our method produces smooth, uniformly triangulated meshes that accurately capture the structural features of interest. This repository provides the source code, which can be comiled into a plugin for [ParaView](https://www.paraview.org).

## Installation

The MidSurfer plugin requires the following software:

* The [CMake](https://cmake.org/) build system, and
* [ParaView source](https://www.paraview.org/download/?filter=Sources) for building the plugin.

In order to build the MidSurfer plugin, ParaView needs to be compiled from sources manually using these
[CMake instructions](https://gitlab.kitware.com/paraview/paraview/-/blob/master/Documentation/dev/build.md).

The MidSurfer plugin can be installed using the standard CMake procedure:

1. Create a build folder
2. Launch CMake with this repository as source folder
3. Configure the project with `-DParaView_DIR=<path to ParaView build folder>`.

## Using the MidSurfer Plugin

Once the plugin has been built, it can be loaded via `Tools->Manage Plugins...`, pressing `Load New..` and selecting the plugin binary (`MidSurfer.[so|dll]`, depending on the operating system). This adds an entry `MidSurfer` in the `Filters` menu which provides the following filters:

* `Extract Midsurface (fast)`: The basic algorithm, which can be used on segmentation masks with well defined structures (like `Fig7.vti`, `Fig13A.vti`, and `Fig13B.vti` in the data folder, corresponding to Fig.7, Fig13(A), and Fig.13(B) in the paper).
* `Extract Midsurface (TODO)`: The complete algorithm, which is slower but also works on segmentation masks with very thin structures (like `Fig11-IMM.vti`, and `Fig11-OMM.vti` in the data folder, corresponding to Fig.11 in the paper).

By toggling the advanced properties button (little gear icon on the top right in the Properties panel), parameters which are automatically set by the algorithm can be accessed for illustration and exploration of the algorithm (disclaimer: not all combinations have been tested or might even make sense).

In addition, the plugin provides the following filters, which are used in the MidSurfer algorithm itself. These filters can be used to explore individual steps of the algorithm, but might be useful as standalone filters as well:

* `Compute Connected Commponents in Binary Image`: TODO
* `Compute Eigen Vector Field`: TODO
* `Extract Center Line`: TODO
* `Generate Point Cloud From Segmentation Mask`: TODO
* `Signed Distance Field`: TODO
* `Zipper Triangulation`: TODO

## Data sets

Data sets used in the paper and summarized in Table 1 therein can be found in the [data](data/) folder.

## Method

![Method](images/method.png)

TODO: link to paper web page
