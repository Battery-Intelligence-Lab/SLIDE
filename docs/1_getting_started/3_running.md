---
layout: default
title: Running
nav_order: 3
---


# Running Slide


## Downloading and running the code (Windows users)
Download the zip file of Slide, and extract it. Move the project to the folder where you had put the Eclipse working space.
You then have to import the project in Eclipse.
- In Eclipse, click ‘file’, ‘import’
- Under ‘General’, select ‘Existing projects into Workspace’ (and click on ‘next’)
- Click on ‘browse’ (next to ‘Select root directory’), and brose to the workspace (or the folder where you had downloaded and extracted this project). Select the folder of this project (i.e. the unzipped folder) and click ‘ok’.
- Click on ‘Finish’, and the project should appear in the ‘project explorer window’ (the left of the Eclipse screen).
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
   - The value of the variable MINGW_HOME should be the folder where you installed the MinGW compiler, e.g. ’ C:\TDM-GCC-64
- Tool Chain Editor
   - From the dropdown menu in ‘Current toolchain’, select MinGW GCC
Then you have to ensure the correct settings are enabled. Note that these settings are to get the MinGW GNU g++ compiler, if you want to use a different compiler you have to make your own settings). Right click on the project folder (in eclipse) and select ‘properties’ all the way at the bottom. In the pop-up window on the left, extend ‘C/C++ Build’ (click on the white triangle) and in the drop-down menu ‘Settings’. On the tab ‘tool settings’ check the following:
- GCC C++ Compiler: in the box called ‘Command’ write the following text (without the quotation marks): ‘ g++ -std=c++11 ‘ (this calls the g++ compiler for C++ and uses C++ version 11)
- GCC C++ Compiler -> Miscellaneous: in the box called ‘other flags’ type the following text (without the quotation marks): ‘ -c -fmessage-length=0 -std=c++11 ’
- GCC C Compiler: in the box called ‘Command’ write the following text (without the quotation marks): ‘ gcc -std=c++11 ‘ (this calls the gcc compiler for C and uses C version 11)
- GCC C Compiler -> Miscellaneous: in the box called ‘other flags’ type the following text (without the quotation marks): ‘ -c -fmessage-length=0 -std=c++11 ’
- MinGW C++ linker: in the box called ‘Command’ type following text (without the quotation marks): ‘ g++ ‘
- MinGW C++ Linker -> Miscellaneous: in the box called ‘Linker flags’ type the following text (without the quotation marks): ‘ -Wl,--stack,8000000000 ‘(this option increases the amount of RAM that the simulation can use. If your computer doesn’t have much RAM you might have to reduce this number; if the simulation runs out of RAM, it will stop working but no error message will be printed, i.e. you have no clue why it stopped).
Now the code should be ready to use
- You can expand folders, the c++ code is in the subfolder ‘src’. By double clicking on ‘Main.cpp’ you will open the code of the main-function, where you can select what you want to simulate by uncommenting the line which calls the function. You uncomment by removing the double backslash (‘\\’) at the start of a line. You comment something in by adding a double backslash in front of it. E.g. to simulate a few CCCV cycles, uncomment the line starting with ‘//CCCV(M, pref, deg, cellType, verbose) in the code block ‘cycling function calls’. Save the code after you have done so. An elaborate explanation about this is given in the next chapters.
- Then you have to build the code by clicking on the hammer-icon in Eclipse (or press ‘CTRL’ + ‘b’). The first time you do this, it might take a while (over a minute). But because builds are incremental, it should go much faster from then onward (a few seconds). If you get an error message that the command ‘g++’ or ‘mingw’ could not be found, it means your computer doesn’t have a c++ compiler. In that case, follow the steps from the section ‘installing a C++ compiler’ above.
- Then you can run the code by clicking on the run-button (the green circle with a white triangle inside it). Or you can run it by expanding the subfolder ‘Binaries’, right-clicking on the .exe-file in there, selecting ‘run as’ and ‘1 Local C/C++ Application’
   - As said, a C++ always starts executing the ‘main’ function in Main.cpp. So if all function calls are still commented, nothings is going to happen. If you have uncommented all function calls, many simulations will be done after each other.
- Then the simulation should start.
- The results are written in csv files, often in subfolders of the project folders. The name of the folder is an identifier indicating which degradation models were used for the simulation.
- Various Matlab functions are written to read the csv files and plot the outcomes. E.g. to read the results from the simulation of the CCCV cycles, execute the script readCCCV.m in Matlab.
Every time you do a simulation, one or more subfolders are created and the results of the simulation are written in those folders. The name of the subfolder consists of three parts (in this order)
- the prefix: in main you have to define a string called prefix. The name of all subfolders will start with the value of this string, followed by an underscore
- the degradation identification: a series of numbers will indicate which degradation models were used during the simulation (the string is generated by the function print_DEG_ID in Degradation.cpp or by the Matlab script printDEGID.m). identifiers of the same mechanism are separated by a hyphen (-), while identifiers of different mechanisms are separated by an underscore (‘_’). E.g. 1-0_2-3-1_1-4_1 means we use
   - 1-0 SEI model 1, no porosity changes due to SEI
   - 2-3-1 CS models 2 and 3, decrease the diffusion constant according to model 1
   - 1-4 LAM models 1 and 4
   - 1 Li lating model 1
- the identification string specified in the function you are running (e.g. CCCV if you simulate the CCCV cycles).
The code forbids overwriting of previous results. If you do a second simulation, you have two options: either you remove the folders with the old data, or you change the value of ‘prefix’ (e.g. from ‘0’ to ‘1’) such that the data will be written in folders with a different name.

## Downloading and running the code (Linux users)
Download the zip file of Slide, and extract it.
In the terminal, navigate to the project folder (the main folder, not src). The first time you build the code, you need to type three commands
```bash
cmake . # (configure using Cmake and CMakeLists.txt)
make    # (build)
./slide # (run)
```
You might get a seg fault because the program can not access enough memory to store all its data. If you don’t know how to solve this, you can try typing
ulimit -S -s 8000000000
However, this command defines the amount of memory ALL programs can use, not just Slide. Therefore, use this with care
Every time you do a simulation, one or more subfolders are created and the results of the simulation are written in those folders. The name of the subfolder consists of three parts (in this order)
- the prefix: in main you have to define a string called prefix. The name of all subfolders will start with the value of this string, followed by an underscore
- the degradation identification: a series of numbers will indicate which degradation models were used during the simulation (the string is generated by the function print_DEG_ID in Degradation.cpp or by the Matlab script printDEGID.m). identifiers of the same mechanism are separated by a hyphen (-), while identifiers of different mechanisms are separated by an underscore (‘_’). E.g. 1-0_2-3-1_1-4_1 means we use
   - 1-0 SEI model 1, no porosity changes due to SEI
   - 2-3-1 CS models 2 and 3, decrease the diffusion constant according to model 1
   - 1-4 LAM models 1 and 4
   - 1 Li lating model 1
- the identification string specified in the function you are running (e.g. CCCV if you simulate the CCCV cycles).
Various Matlab functions are written to read the csv files and plot the outcomes. E.g. to read the results from the simulation of the CCCV cycles, execute the script readCCCV.m in Matlab.
The code forbids overwriting of previous results. If you do a second simulation, you have two options: either you remove the folders with the old data, or you change the value of ‘prefix’ (e.g. from ‘0’ to ‘1’) such that the data will be written in folders with a different name.
To run later simulations, you only need the last two commands from the terminal

```bash
make
./slide
```