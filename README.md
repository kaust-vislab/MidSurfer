# MidSurfer
Paraview Plugin for extracting mid-surfaces from volumetric segmentation masks. 

![Teaser](images/teaser.png)

## Overview

TODO

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

* `Extract Midsurface (fast)`: The basic algorithm, which can be used on segmentation masks with well defined structures (like `Fig7.vti` and `Fig13.vti` in the data folder, corresponding to Fig.7 and Fig.13 in the paper).
* `Extract Midsurface (TODO)`: The complete algorithm, which is slower but also works on segmentation masks with very thin structures (TODO: add data sets to folder or provide download link)

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
