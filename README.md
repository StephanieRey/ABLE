# ABLE

This repository contains the current MATLAB implementation of 'ABLE' (an Activity-Based LEvel set segmentation algorithm), which detects the location of neurons from two-photon calcium imaging data. For technical details, see the [paper](https://doi.org/10.1101/190348). 

## Overview

In our framework, multiple coupled active contours evolve, guided by a model-based cost function, to identify cell boundaries. An active contour seeks to partition the local region into two subregions, a cell interior and exterior, in which all pixels have maximally 'similar' temporal activity. To see an example of an active contour evolving to segment a cell, see [evolution_video.avi](https://github.com/stephanierey/ABLE/evolution_video.avi). 

The algorithm is flexbile - we incorporate no prior assumptions regarding cell shape or temporal activity. This allows us to detect a [diverse array](https://github.com/stephanierey/ABLE/results_figures/large_and_small_ROIs.png) of ROIs, as illustrated in a demo on mouse in vivo imaging data. In this dataset, ABLE detects both cell bodies and cross-sections of dendrites. 

Algorithm runtime is virtually independent of video length - videos from long recording sessions can be processed if loaded in unsigned integer format. In the current implementation, the runtime is linear with the number of cells, approximately 15 minutes is required to segment 400 cells. In the future, we will exploit the natural parallelism in the algorithm to reduce this. 

## Demo

In this folder we provide a [demo](https://github.com/stephanierey/ABLE/demo.m) of ABLE along with some [documentation](https://github.com/stephanierey/ABLE/ABLE_documentation.pdf). For the demo we use one set of videos from the following publicly available dataset:  

Simon Peron, Jeremy Freeman, Vijay Iyer, Caiying Guo, Karel Svoboda (2015);  Volumetric calcium imaging data recorded during performance of a single whisker object localization task, sampling activity in the majority of the relevant superficial barrel cortex neurons (75 %, 12,000 neurons per mouse). CRCNS.org. [doi:10.6080/K0TB14TN](http://dx.doi.org/10.6080/K0TB14TN).

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* To regularize the Level Set functions we use Li et al.'s [DRLSE](http://dx.doi.org/10.1109/TIP.2010.2069690).




