---
layout: default
title: Installation
nav_order: 2
---

# Installation
This chapter explains installation. 

## Installing a C++ compiler (Windows users only)

<p style='text-align: justify;'>
Some computers and operating systems come with a c++ compiler. If you don’t have one yet, you won’t be able to ‘build the code’ 
(see below, you will get an error message in Eclipse that the command ‘g++’ or ‘minGW’ couldn’t be found). In that case, you have to install a C++ compiler. 
Any C++ compiler should work, but this code has been developed using g++ (from mingw). The steps below are to install the compiler from TDM, which can do both 32 and 64 bit.
</p>

- Download the [TMD-GCC](https://sourceforge.net/projects/tdm-gcc/) installer
- Run the installer, and select ‘Create: create a new TDM-GCC installation)
- Select the ‘MinGW-w64/TDM (32-bit and 64-bit)’ installation
- Install it in the standard directory. DO NOT MODIFY THE INSTALLATION DIRECTORY, else Eclipse might not find the compiler.
- Use any ‘mirror’ you like (this is the server from which the code will be downloaded)
- In ‘Choose Components’ you don’t have to change anything. Normally, the required packages are installed (the essential package is gcc -> c++ but you need the other packages to make the installation work)
- Click install. It might take a couple of minutes.


## Installing Eclipse (windows users only)

Slide is mostly written in C++ and requires a programming interface. Any c++ programming environment will work, below are the steps to install the environment from Eclipse.

1. Install a Java Development Kit (JDK) or a Java Runtime Environment (JRE): Eclipse needs a JDK or JRE to install and run. If you haven’t installed a JDK or JRE yet, 
it can be downloaded from [Oracle](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html). If you already have a JDK or JRE, you can skip this step. 
If you are unsure, skip this step, and the Eclipse installer will notify you that you first have to install a JDK or JRE. You might have to restart your computer after installing 
JDK or JRE before the Eclipse installer will recognize it. Please not that not all versions of Eclipse support all versions of JDK or JRE. The Eclipse installer 
will inform you which Java version to install. 

2. Download a version of [Eclipse](https://www.eclipse.org/downloads/packages/). This project is made using Eclipse Neon, but later versions should work as well. You can either download only ‘Eclipse IDE for C/C++ Developers’ or the full Eclipse-version (on the right).  Make sure to download the 32-bit version if your computer has a 32 bit version of Windows. (note: 64 bit is identified by ‘64’ or ‘x64’ while the 32 bit version is often identified by ‘86’ or ‘x86’. 
   The JDK or JRE version must have the same bit-version as the Eclipse version you are installing. Sometimes the Eclipse installer doesn’t recognise the JDK or JRE you have installed. 
   You can try the following steps

    - Try restarting the computer.
    - Double check that you have the same bit-version (i.e. 32-bit Java and 32-bit Eclipse or 64-bit Java and 64-bit Eclipse), 
    and your Java version is the one needed by Eclipse (not the latest version of Java available on the Oracle website).
    - If that doesn’t work, uninstall the Java you had installed in step 1 and try re-installing it, potentially try a different version (7/8 and JDK/JRE). 
    Usually Eclipse finds the Java installation after a few re-installations.

3. Execute the installer and if you downloaded the full Eclipse-version, you have to select the ‘Eclipse IDE for C/C++ Developers’ when asked which Eclipse program to install.

To check afterwards whether a 32 or 64 bit version is installed, open Task Manager and find the Eclipse process (in the tab ‘Processes’). If the process is called ‘eclipse.exe’, it is the 64-bit. 
If it is called ‘eclipse.exe *32’ it is a 32-bit version. can be observed from Task manager (ctrl-shift-esc). In case it is the 32 bit version, behind ‘eclipse.exe’ is written ‘(32 bit)’

When you run Eclipse, it will ask you to ‘select a directory as a workspace’. This directory will be the folder where all your C++ projects are stored, so choose an empty folder you can easily find back.

## Installing CMake (Linux users only)

Please download and install [CMake](https://cmake.org/download/) 
- Download the zipped file (.tar.gz) and unzip it
- In the terminal, navigate to the folder where you downloaded it using `cd` command
- Type the following commands in the terminal:

```bash
./configure
make
make install
```