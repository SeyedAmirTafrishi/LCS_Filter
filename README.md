# Line-Circle-Square (LCS): A Multilayered Geometric Filter for Edge-Based Detection
## Introduction

Line-Circle-Square (LCS) filter applies detection, tracking and learning to each defined expert (Line, Circle and Square experts) to obtain more information for judging scene without over-calculation. The interactive learning feed between each expert creates a minimal error that fights against overwhelming landmark signs in crowded scenes without mapping. Our experts basically are dependent on trust factors' covariance with geometric definitions to ignore, emerge and compare detected landmarks.

Compared with other related filters in the literature, the proposed LCS filter (and the LC filter) has the following unique advantages: (1) it reduces computation demand; (2) it has the ability to minimize the problem of overconfidence during detection; (3) real-time process for detecting abnormal behaviors at outside world such as partial detection of incoming objects toward the camera/moving vehicle; (4) primary detection with geometrical computation which creates different level of information, i.e., low (edges) to high (layers) for mapping and localization; (5) the multi-layer nature makes it suitable for real-time processing with potential to be executed in parallel.

![](results/LCS_demo.gif)

A video of this work can be found at: https://www.youtube.com/watch?v=hh6B4aF7jDI


## Autonomous Car Example

The code of the autonomous car example in the evaluation can be found in the `autonomous_car_example` folder.

The dataset we used is Oxford RoboCar dataset and to be more specific, we used the following two datasets:

- Day: https://robotcar-dataset.robots.ox.ac.uk/datasets/2015-08-17-13-30-19/
- Night: https://robotcar-dataset.robots.ox.ac.uk/datasets/2014-12-17-18-18-43/

You can find the demonstration videos of the autonomous car example in the following link:

https://drive.google.com/drive/folders/1iwPiapS100yQHR1NgSOOHTOxFetvqOsA

The two videos are:

- RAS_NG_D.mp4: the night scenario
- RAS_DH_D.mp4: the day scenario

## Usage

To use the code, run the script `main_offline.m` from the MATLAB command line. This will run the filter with some example pictures in offline mode. 

LCS filter also supports to run in real-time. To use the real-time version of the filter, use the script `main_rt.m` (debugging stage).

## Dependencies

There is no dependencies for the offline code. However, in order to run the real-time version, the following packages are required:

- MATLAB `mobiledev` package
- MATLAB `webcam` package
- Requires a web camera and an IMU sensor. 

## Credits

This project depends on the *fast9* edge detection algorithm, which is developed by **Edward Rosten**, **Reid Porter** and **Tom Drummond**.

## Contributors

- **Seyed Amir Tafrishi**, Kyushu University, Japan
- **Xiaotian Dai**, Department of Computer Science, University of York, UK
- **Vahid E. Kandjani**, University College of Nabi Akram, Iran

##  Publications

- Tafrishi, S.A., Xiaotian, D., Kandjani, V.E., 2021. Line-Circle-Square (LCS): A Multilayered Geometric Filter for Edge-Based Detection. Robotics and Autonomous Systems (RAS).
- Tafrishi, S.A. and Kandjani, V.E., 2017, October. [Line-Circle: A Geometric Filter for Single Camera Edge-Based Object Detection](https://ieeexplore.ieee.org/abstract/document/8466193). In *2017 5th RSI International Conference on Robotics and Mechatronics (ICRoM)* (pp. 588-594). 
