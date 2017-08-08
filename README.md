# Multizone Soundfield Reproduction (MSR) Toolbox

[![GitHub release](https://img.shields.io/github/release/JacobD10/MSR.svg?style=flat-square)](https://github.com/JacobD10/MSR/releases)
[![GitHub commits](https://img.shields.io/github/commits-since/JacobD10/MSR/1.0.0.svg?style=flat-square)](https://github.com/JacobD10/MSR/commits/master)
[![Github All Releases](https://img.shields.io/github/downloads/JacobD10/MSR/total.svg?style=flat-square)](https://github.com/JacobD10/MSR)
[![Github file size](https://img.shields.io/github/size/JacobD10/MSR/release/release.zip.svg?style=flat-square)](https://github.com/JacobD10/MSR/blob/master/release/release.zip)
[![license](https://img.shields.io/github/license/JacobD10/MSR.svg?style=flat-square)](https://github.com/JacobD10/MSR/blob/master/LICENSE)
[![Twitter URL](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2FJacobD10%2FMSR&via=_JacobDonley&text=Check%20out%20the%20Multizone%20Soundfield%20Reproduction%20toolbox%20for%20%23MATLAB%21&hashtags=software%20%23code%20%23audio)

## [Documentation Available Online](https://www.soundzones.com/jdonley/MSR/doc/html/)

## Installation
Installing the MSR Toolbox will by default override your MATLAB [startup](https://au.mathworks.com/help/matlab/ref/startup.html) routine. If this is **not** functionality you would like then you may at anytime rename `startup.m` to `initialise.m`.

1) Download the latest release of the MSR Toolbox and extract it to your desired location.
2) Open MATLAB and navigate to the MSR Toolbox directory.
3) Type `edit('+All_Settings' filesep 'Global_Settings.m'])` and press enter.
4) Read the comments in the `Global_Settings` file and set the parameters as desired.
5) Restart MATLAB.
6) The startup procedure will display a list of dependencies which are not installed currently and that are required for the toolbox to run. (If you renamed `startup.m`, run `initialise.m`)
	* Install the dependencies to a location on the MATLAB path and restart MATLAB.
7) Congratulations, the installation should be complete!

## Example
When you start MATLAB (or run `startup.m` or `initialise.m`) and editor will open with an MSR system evaluation script. This is a helpful script for running a complete MSR system from start to finish. There are some example system settings files in `+Current_Systems>+Archived`. When navigating among the files and folders click the `+` next to the folder to expand the contents rather than double clicking (which changes the current directory). 

1) Edit `+Current_Systems>loadCurrentSRsystem` and change the variable `Current_System` to a system settings filename of your choosing (without file extension) which is located under the filepath you set in the `Global_Settings` file (by default this is `+Current_Systems`).
2) Edit the system settings file which you specified in step 1 to suit your needs.
3) Run the entire MSR system evaluation script `+Current_Systems>Evaluate_Current_System` or run section by section in order.
4) Done. Each step should complete without errors and all data should be saved under the drives specified in the system settings file. 

There are also small examples on individual MSR classes such as `loudspeaker_setup.m` and `multizone_soundfield_OBE.m`. For in-depth descriptions of classes, functions and scripts please [read the documentation](https://www.soundzones.com/jdonley/MSR/doc/html/).


## Authors
### [*Jacob Donley*](https://www.soundzones.com/jdonley)&nbsp;&nbsp; :email: [Email](mailto:jrd089@uowmail.edu.au)&nbsp;&nbsp; [![Twitter Follow](https://img.shields.io/twitter/follow/_JacobDonley.svg?style=social&label=Follow)](https://twitter.com/intent/user?screen_name=_JacobDonley)

