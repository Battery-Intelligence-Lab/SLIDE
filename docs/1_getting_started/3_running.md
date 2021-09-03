---
layout: default
title: Running
nav_order: 3
---


# Running Slide

After installing SLIDE and compiling the code, running it is as easy as executing the executable output in ```Release``` (by default) folder. 

To have a nice working environment, we suggest using [Visual Studio Code (VScode)](https://code.visualstudio.com/) for all platforms. 

- After downloading VScode, you may install required extensions using ```C/C++ Extension Pack``` extension. 
- After restarting VScode, you can open the SLIDE folder in VScode and let it configure. 
- Then you can scan and select the installed kits (i.e., compilers) run by clicking launch button on the blue bar at the bottom. 

## Using Eclipse (Windows users)
You first have to import the project in Eclipse.
- In Eclipse, click ```file```, ```import```
- Under ```General```, select ```Existing projects into Workspace``` (and click on ```next```)
- Click on ```browse``` (next to ```Select root directory```), and brose to the workspace (or the folder where you had downloaded and extracted this project). Select the folder of this project (i.e. the unzipped folder) and click ```ok```.
- Click on ‘```Finish```, and the project should appear in the ```project explorer window``` (the left of the Eclipse screen).
Sometimes, Eclipse does not recognise the downloaded folder as an existing project. In that case, make a new project in your workspace
- Right click in the workspace area of Eclipse
- New -> C++ project
- With all the standard settings, except that you should use MinGW GCC as toolchain
- Finish
- Right click on the new project folder (in eclipse)
- Import -> file system -> browse to the unzipped folder and import all files
In case the MinGW GCC toolchain is not recognized or not present in the eclipse option, follow these steps:
- Right click on the project
- Properties
- C/C++ Build
- Environment
   - The value of the variable MINGW_HOME should be the folder where you installed the MinGW compiler, e.g. ```C:\TDM-GCC-64```
- Tool Chain Editor
   - From the dropdown menu in ```Current toolchain```, select ```MinGW GCC```
Then you have to ensure the correct settings are enabled. Note that these settings are to get the MinGW GNU g++ compiler, if you want to use a different compiler you have to make your own settings). Right click on the project folder (in Eclipse) and select ```properties``` all the way at the bottom. In the pop-up window on the left, extend ```C/C++ Build``` (click on the white triangle) and in the drop-down menu ```Settings```. On the tab ```tool settings``` check the following:
- GCC C++ Compiler: in the box called ```Command``` write the following text: ```g++ -std=c++17``` (this calls the g++ compiler for C++ and uses C++ standard 2017)
- GCC C++ Compiler -> Miscellaneous: in the box called ```other flags``` type the following text (without the quotation marks): ``` -c -fmessage-length=0 -std=c++17 ```
- MinGW C++ linker: in the box called ```Command``` type following text: ```g++```
Now the code should be ready to use.
- You can expand folders, the C++ code is in the subfolder ```src```. By double clicking on ```main.cpp``` you will open the code of the main function, where you can select what you want to simulate by uncommenting the line which calls the function.
- Then you have to build the code by clicking on the hammer icon in Eclipse (or press ```CTRL + b```). The first time you do this, it might take a while (over a minute). But because builds are incremental, it should go much faster from then onward (a few seconds). If you get an error message that the command ```g++``` or ```MinGW``` could not be found, it means your computer doesn’t have a C++ compiler. In that case, follow the steps from the previos section [installation](../1_getting_started/2_installation.html).
- Then you can run the code by clicking on the run-button (the green circle with a white triangle inside it). Or you can run it by expanding the subfolder ‘Binaries’, right-clicking on the .exe-file in there, selecting ```run as``` and ```1 Local C/C++ Application```
- Then the simulation should start.

## Quick configuration

- As said, a C++ always starts executing the ```main``` function in ```main.cpp```. So if all function calls are still commented, nothings is going to happen. If you have uncommented all function calls, many simulations will be done after each other. You uncomment by removing the double backslash (```\\```) at the start of a line. You comment something in by adding a double backslash in front of it. E.g. to simulate a few CCCV cycles, uncomment the line starting with ```//CCCV(M, pref, deg, cellType, verbose)``` in the code block ```cycling function calls```. Save the code after you have done so. An elaborate explanation about this is given in the next chapters.
- The results are written in ```*.csv``` files, often in subfolders of the project folders. The name of the folder is an identifier indicating which degradation models were used for the simulation.
- Various MATLAB functions are written to read the ```*.csv``` files and plot the outcomes. E.g. to read the results from the simulation of the CCCV cycles, execute the script ```readCCCV.m``` in MATLAB.
Every time you do a simulation, one or more subfolders are created and the results of the simulation are written in those folders. The name of the subfolder consists of three parts (in this order)
- the prefix: in main you have to define a string called prefix. The name of all subfolders will start with the value of this string, followed by an underscore
- the degradation identification: a series of numbers will indicate which degradation models were used during the simulation (the string is generated by the function print_DEG_ID in degradation.cpp or by the MATLAB script printDEGID.m). identifiers of the same mechanism are separated by a hyphen ```-```, while identifiers of different mechanisms are separated by an underscore ```_```. E.g. ```1-0_2-3-1_1-4_1``` means we use
   - 1-0 SEI model 1, no porosity changes due to SEI
   - 2-3-1 CS models 2 and 3, decrease the diffusion constant according to model 1
   - 1-4 LAM models 1 and 4
   - 1 Li plating model 1
- the identification string specified in the function you are running (e.g. CCCV if you simulate the CCCV cycles).
