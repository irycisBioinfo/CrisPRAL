# CrisPR Analysis. (CrisPRAL)
### A tool for CrisPR variant analysis of genome editing experiments and detection of genetical mosaicism in rare diseases.

The purpose of this tool is to allow for the identification, sorting and cuantification of alelic diversity given
crude .fastq data obtained through Next Generation Sequencing (NGS). It can both be used for the analysis of the 
alelic multiplicity that can arise from CrisPR editing protocols in animal models, as the characterization of somatic
mosaicism associated to the natural genetic therapy phenomenom.

We believe this project to be of interest in the attempts of translating CRISPR editing system from experimental 
to clinical practice.

## Getting Started.

*The main purpose of this app is to be accesible via web domain, but for the time being the app must be run locally in order
to function properly.*

This instructions will attempt to aid in the deployment of the app in your local system.

### Prerequisites

* Linux operating system
* interpreter:
  * R
  * perl
  * python
  
### Dependancies

* Cutadapt

In Linux sub-systems
```shell
	$ sudo apt install cutadapt
```
or
```shell
	$ sudo apt install python3-cutadapt
```
* Google Chrome & chromium-browser

```shell
	$ wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
	$ sudo dpkg -i google-chrome-stable_current_amd64.deb
```
  
### Installing

Clone this repository into your system:

```shell
	$ git clone https://github.com/irycisBioinfo/CrisPRAL.git    
```

or download .zip and unpack.

### Testing

Load app.R file through a R interpreter. 

* You can use either RStudio by opening the file and running the app through the "Run App" option at the top right corner of the script editor.

* Sourcing directly the file through RStudio's console.
```R
	source( "path/to/file/app.R" )
```
* Executing the file through your terminal window
```shell
	$ Rscript path/to/file/app.R
```
When executing through the terminal window you will have to manually copy the ip adress displayed and paste it in your browser of preference.

You will find yourself with the main page presenting the mosaic finder functionality. Enter the obligatory fields, Read 1, Read 2, and the reference. You can find demo files for these fields in */Trials* directory.

The rest are optional options but demo files are presented aswell for the Adapters filtering and the PCR Primers trimming.
