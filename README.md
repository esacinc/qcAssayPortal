qcAssayPortal
================================================

A Python program to evaluate the data quality of the five experiments in CPTAC assay portal, the format of input files is *.sky.zip.
The experiments are: 
Experiment 1, response curve (or ResponseCurve)
Experiment 2, evaluation of repeatability (or ValidationSamples)
Experiment 3, selectivity (or Seletivity)
Experiment 4, stability (Stability)
Experiment 5, reproducible detection of endogenous analyte (Endogenous)
The current version supports Experiment 1 and Experiment 2.

The work flow of qcAssayPortal is shown in workflow.pdf. 

A report file will be generated which will capture the details of errors and warnings when running the R codes for the individual experiment.
The details of the errors and warnings are listed in issue_categories.pdf

Documentation
-------------

* [Installation](#installation)
  * [Install qcAssayPortal](#install-qcAssayPortal)
* [How to use it](#how-to-use-it)
  * [Command line](#command-line)
* [Changelog](#changelog)
* [Citation](#citation)

Installation
------------

### Install qcAssayPortal

qcAssayPortal is implemented as a Python program running on a Windows platform that requires pre-installation of Skyline https://skyline.ms/project/home/software/Skyline/begin.view and R https://cran.r-project.org/bin/windows/base/  <br /><br />
For Python (v2.7.\*) programming language, it requires Python-related libraries, including pandas(>= v0.21.0) and Jinja2(v2.9.6). <br /><br />
For R (v3.5.\* is recommended) programming language, it requires R-related libraries, including Cairo, evaluate, stringr, plyr, MASS ggplot2 and dplyr. In order to install required libraries, please run install.packages() in the R console. <br /><br />
For Skyline (v4.2 is recommended), it requires pre-installed Skyline command-line interface https://skyline.ms/_webdav/home/software/Skyline/@files/docs/Skyline%20Command-Line%20Interface-3_7.pdf. Or install Skyline administrator downloaded from https://skyline.ms/wiki/home/software/Skyline/page.view?name=install-administator-64. The command-line toolSkylineCmd.exe can be found in the installation directory. <br />

qcAssayPortal can be installed from the source code by pip:<br />
1) Download qcAssayPortal source code from URL and unzip the zipped file folder.<br />
2) If the package of wheel is not installed, run `pip install wheel` to install it.<br />
3) Change the directory to qcAssayPortal's directory and run `python setup.py bdist_wheel` to build a wheel file for the subsequent installation via pip.<br />
4) Run `pip install ./dist/qcAssayPortal-1.0-py2-none-any.whl` to install qcAssayPortal.<br />


How to use it
-------------

### Command line

    Usage: qcAssayPortal [-h] [<args>]

    Example 1:
    qcAssayPortal -ps "C:\Program Files\SkylineDailyRunner\SkylineDailyRunner.exe" -pr "C:\Program Files\R\R-3.3.1\bin\Rscript.exe" -i "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1?0309_MouseV2B1.5_refined_2018-07-03_14-59-18.sky.zip" -e exp1 -t "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp1\meta_data.tsv"
    
    Example 2:
    qcAssayPortal -ps "C:\Program Files\Skyline\SkylineCmd.exe" -pr "C:\Program Files\R\R-3.3.1\bin\Rscript.exe" -i "D:\Skyline_analysis\qcAssayPortal\data\UVicPC_Borchers-MousePlasma_Agilent6490_directMRM-Exp2" -e exp2

    optional arguments:
      -h, --help            show this help message and exit
      -ps <dir required>    the path to SkylineCmd.exe, SkylineRunner.exe, or SkylineDailyRunner.exe in the Windows OS
      -pr <dir required>    the path to Rscript.exe in the Windows OS
      -i input <required> [input <required> ...]
                            two options: 1. A directory where all the *.sky.zip files are located 2.*.sky.zip *.sky.zip ... (at least one input *.sky.zip)
      -e string <required>  the experiment type. Choose one from Options: exp1 , exp2
      -t peptide type file <required> the directory of the file whose first column is the *.sky.zip and second column is peptide type ("purified" or "crude"). When -e is exp1, it must be assigned. Otherwise, please leave it blank.
      -o <dir>              the directory of the outputs (default: current directory)
      --version             show program's version number and exit

Changelog
---------

Citation
--------
