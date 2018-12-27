# WeMe
WeMe is a method to find a consensus subclonal structure from multiple subclonal structures (typically from different methods).  
For a description of the method, see Supplement section XX in the PCAWG Heterogeneity paper (available on BioArxiv at: https://www.biorxiv.org/content/early/2018/07/13/312041)
WeMe requires R and the following libraries: doMC,reshape2, ggplot2,and RColorBrewer.
WeMe has been tested on OS X Version 10.14 with the following package versions:
* R: version 3.5.0
* doMC: version 1.3.5
* reshape2: version 1.4.3
* ggplot2: version 3.1.0
* RColorBrewer: version 1.1-2

As input WeMe requires multiple reconstructions in the PCAWG-11 format (reference).   
WeMe does not require any installation.
The typical operation involves running WeMe in a folder containing all the subclonal structure files, with a directory per-method.  e.g:
```
$ ls -lrtd */
drwxr-xr-x 2 user group 1380352 Jun 11 19:39 method1/
drwxr-xr-x 2 user group 1945600 Jun 11 19:40 method2/
drwxr-xr-x 2 user group  327680 Jun 11 19:42 method3/
drwxr-xr-x 2 user group  327680 Jun 11 19:44 method4/
```
Inside each directory there should a subclonal reconstruction file with the following naming convention:

```
{samplename}_subclonal_structure.txt
```
The `weme_demo` folder includes an example of this layout.

# Running WeMe

```
source("weme.R")
```
Identify all subclonal structure files and the methods
```
sids = find_sids()
```
Generate consensus structure files for all samples found
```
genconsensus(sids,rounddown=FALSE)
```
The above call will generate a consensus subclonal structure file and a plot comparing the independent structures and the consensus structure.
Running the demo should take less than 1 minute on a regular computer.
Example output is provided in the root of the weme_demo folder.

In order to ensure consensus results are consistent with an independently determined purity levels, WeMe can optionally correct the consensus solution. 
Correct the top cluster to match input purity.  Can be run independently.  THe purity/ploidy file must contain columns samplename and purity.
```
correct_purity("purity.demo.txt")
```
