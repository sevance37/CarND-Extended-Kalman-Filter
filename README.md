# CarND-Extended-Kalman-Filter
Estimate the location of a moving object by using an extended Kalman filter on noisy lidar and radar measurements.  The program is written in C++.  This Project is from Udacity's Self-Driving Car Engineer Nanodegree Program.

## Basic Set-up
1. Clone this repo.
2. Make sure [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) is installed.  Two install scripts are included for MAC and Linux.  For Windows use [Windows 10 Bash on Ubuntu](https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/) and following the Linux instructions.
3. Make a build directory: `mkdir build && cd build`
4. Compile: `cmake .. && make`
5. Run the programs: Run `./ExtendedKF `. Open [Term 2 Simulator](https://github.com/udacity/self-driving-car-sim/releases) and run the corresponding page.

## Program Input and output.
INPUT: values provided by the simulator to the c++ program  
["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)  

OUTPUT: values provided by the c++ program to the simulator  
["estimate_x"] <= kalman filter estimated position x  
["estimate_y"] <= kalman filter estimated position y  
["rmse_x"]  
["rmse_y"]  
["rmse_vx"]  
["rmse_vy"]  

## Other Important Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
