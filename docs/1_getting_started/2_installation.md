---
layout: default
title: Installation
nav_order: 2
---

# Installation

SLIDE is yet to offer any binaries or wrappers in any other languages; therefore, the only way to install SLIDE is to compile it from C++ source files. 

## Building from the source

SLIDE aims to be compatible with different compilers and platforms. You may easily install SLIDE using CMake (although it is not an absolute requirement). Therefore, you need a suitable C++ compiler (preferably [GCC](https://gcc.gnu.org/)) and [CMake](https://cmake.org/) to follow this installation guide.   


### Linux

Generally, both gcc and CMake are installed by default on Linux platforms. However, in some cases you may need to install them. For example, Ubuntu 18.04 comes with an older compiler that does not support some of the functionalities in this code directly. Therefore, you may want to install a newer version of GCC as shown [here](https://linuxize.com/post/how-to-install-gcc-compiler-on-ubuntu-18-04/). 

1. Then you may download the repository as [*.zip file](https://github.com/davidhowey/SLIDE/archive/refs/heads/master.zip) or clone it using following command: 
```bash
git clone https://github.com/davidhowey/SLIDE.git
```
2. After downloading source files, you need to create a build folder and go into the build folder.
```bash
cd SLIDE  # Go into the main directory of source files.
mkdir build # Create a build folder if it is not present.
cd build # Go into the build directory. 
```
3. Then you may create Makefiles by 
```bash
cmake -G "Unix Makefiles" .. 
```
4. Compile the files:
```bash
cmake --build . # Assuming that you are still in the build folder. 
```

```note
Then the executable will be ready to run at ```../Release/slide```. By default ```CMAKE_BUILD_TYPE``` is set to ```Release```. If you want, you may also use one of the other options as ```Debug```, ```Release```, ```RelWithDebInfo```, ```MinSizeRel```.To build using alternative build type you may explicitly define ```CMAKE_BUILD_TYPE``` variable. For example, for building with debug mode you may use the following command. For further information on using CMake, please refer to [CMake guide](https://cmake.org/cmake/help/git-stage/index.html)  


```cmake --build . -DCMAKE_BUILD_TYPE=Debug```
```


### Windows

On Windows platforms, you probably need to install CMake and a C++ compiler. 

1. Install the latest version of [CMake binary](https://cmake.org/download/#latest).
2. You may then install a compiler of your choice. However, if you are going to use GCC, we suggest installing [TDM-GCC](https://jmeubank.github.io/tdm-gcc/download/) (preferably **MinGW-w64 based edition**). Otherwise, you may also install [Visual Studio Community](https://visualstudio.microsoft.com/vs/community/) IDE which comes with its own compiler.
3. You may steps 1--4 on [Linux section](#linux) except in step 3, you should write following command: 
```bash
cmake -G "MinGW Makefiles" ..  # if you use MinGW GCC compiler.
```
or you can create ```*.sln``` file as well as build files via following command if you use Visual Studio Community. 
```bash
cmake -G "Visual Studio 16 2019" ..  # if you use Visual Studio's compiler.
```

```note
    If you are using Visual Studio Community, you may also open the folder in Visual Studio directly, without using CMake. 
    See [here](https://docs.microsoft.com/en-us/visualstudio/ide/develop-code-in-visual-studio-without-projects-or-solutions?view=vs-2019) for detailed explanation.
```


**Using Eclipse:**


SLIDE was originally developed using Eclipse IDE and it is also possible to use Eclipse for running SLIDE. After installing [TDM-GCC](https://jmeubank.github.io/tdm-gcc/download/) you may follow these steps: 

1. Install a Java Development Kit (JDK) or a Java Runtime Environment (JRE): Eclipse needs a JDK or JRE to install and run. If you haven’t installed a JDK or JRE yet, it can be downloaded from [Oracle](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). If you already have a JDK or JRE, you can skip this step. If you are unsure, skip this step, and the Eclipse installer will notify you that you first have to install a JDK or JRE. You might have to restart your computer after installing JDK or JRE before the Eclipse installer will recognize it. Please not that not all versions of Eclipse support all versions of JDK or JRE. The Eclipse installer will inform you which Java version to install.

2. Download a version of [Eclipse](https://www.eclipse.org/downloads/packages/). This project is made using Eclipse Neon, but later versions should work as well. You can either download only *Eclipse IDE for C/C++ Developers* or the full Eclipse-version (on the right).  Make sure to download the 32-bit version if your computer has a 32 bit version of Windows. (note: 64 bit is identified by ‘64’ or ‘x64’ while the 32 bit version is often identified by ‘86’ or ‘x86’. 
   The JDK or JRE version must have the same bit-version as the Eclipse version you are installing. Sometimes the Eclipse installer doesn’t recognise the JDK or JRE you have installed. You can try the following steps:
    - Try restarting the computer.
    - Double check that you have the same bit-version (i.e. 32-bit Java and 32-bit Eclipse or 64-bit Java and 64-bit Eclipse), 
    and your Java version is the one needed by Eclipse (not the latest version of Java available on the Oracle website).
    - If that doesn’t work, uninstall the Java you had installed in step 1 and try re-installing it, potentially try a different version (7/8 and JDK/JRE). 
    Usually Eclipse finds the Java installation after a few re-installations.

3. Execute the installer and if you downloaded the full Eclipse-version, you have to select the *Eclipse IDE for C/C++ Developers* when asked which Eclipse program to install.

<!--- // COMMENT out
To check afterwards whether a 32 or 64 bit version is installed, open Task Manager and find the Eclipse process (in the tab ‘Processes’). If the process is called ‘eclipse.exe’, it is the 64-bit. 
If it is called ‘eclipse.exe *32’ it is a 32-bit version. can be observed from Task manager (ctrl-shift-esc). In case it is the 32 bit version, behind ‘eclipse.exe’ is written ‘(32 bit)’

 When you run Eclipse, it will ask you to ‘select a directory as a workspace’. This directory will be the folder where all your C++ projects are stored, so choose an empty folder you can easily find back.
-->

### macOS

1. Install the latest version of [Xcode](https://developer.apple.com/support/xcode/).
2. Install command line tools, [Homebrew](https://brew.sh/) and CMake by executing following commands on the terminal: 
```bash
xcode-select --install
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install cmake
```
3. Follow steps 1--4 on [Linux section](#linux)
