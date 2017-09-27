
# ABLE

%%%% NOTICE: demo is being updated, will be online by 01.10.2017 %%%%%


This repository contains the current MATLAB implementation of 'ABLE' (an Activity-Based LEvel set segmentation algorithm), which detects the location of neurons from two-photon calcium imaging data. For technical details, see the paper at https://doi.org/10.1101/190348. 

In our framework, multiple coupled active contours evolve, guided by a model-based cost function, to identify cell boundaries. An active contour seeks to partition the local region into two subregions, a cell interior and exterior, in which all pixels have maximally 'similar' temporal activity. The algorithm is flexbile - we incorporate no prior assumptions regarding cell shape or temporal activity. This allows us to detect a diverse array of ROIs, as illustrated in the demo (see demo.m). To initialise the algorithm, we identify local peaks in the correlation and/or mean image. Algorithm runtime is virtually independent of video length - videos from long recording sessions can be processed if loaded in unsigned integer format. In the current implementation, the runtime is linear with the number of cells, approximately 15 minutes is required to segment 400 cells. In the future, we will exploit the natural parallelism in the algorithm to reduce this. 

See the file 'ABLE documentation.pdf' for guidance on how to use the algorithm. See 'evolution_video.avi' for an example of the process of segmenting one cell. 
