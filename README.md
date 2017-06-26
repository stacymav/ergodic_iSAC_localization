# ergodic_iSAC_localization

C++ code implementing ergodic iSAC for target localization: an ergodic control algorithm that uses complex agent dynamics to track and localize random number of targets by exploring an expected information density in real time. The current implementation uses 12-dimensional quadrotor dynamics and an Extended Kalman Filter to track two moving targets with bearing-only measurements. See videos of the results here: https://vimeo.com/stacymav

# Dependencies

The code requires the Boost and Eigen libraries.

# To compile and run

--- Update Makefile.txt with local Boost and Eigen paths
--- "make"
--- Run Quad_Euler.exe
--- Plot resulting trajectories in Matlab, using ./data/plots_matlab.m

# Customization

All possible changes, e.g. agent dynamics, number and motion of targets, estimation filter etc., can be made by updating the files included in the "user" folder.
