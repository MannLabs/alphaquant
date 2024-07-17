<!---
![Pip installation](https://github.com/MannLabs/alphaquant/workflows/Default%20installation%20and%20tests/badge.svg)
![GUI and PyPi releases](https://github.com/MannLabs/alphaquant/workflows/Publish%20on%20PyPi%20and%20release%20on%20GitHub/badge.svg)
[![Downloads](https://pepy.tech/badge/alphaquant)](https://pepy.tech/project/alphaquant)
[![Downloads](https://pepy.tech/badge/alphaquant/month)](https://pepy.tech/project/alphaquant)
[![Downloads](https://pepy.tech/badge/alphaquant/week)](https://pepy.tech/project/alphaquant)
-->


# AlphaQuant
AlphaQuant is an open-source Python package for sensitive proteomics quantification. You can process MS data analyzed by Spectronaut, DIANN, [AlphaPept](https://github.com/MannLabs/alphapept) or MaxQuant using a Graphical User Interface (GUI) or the python package. The current focus is on the comparison of two biological conditions (i.e. "making volcano plots"), with multi-condition functionality to be added soon.

It is part of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/). To enable all hyperlinks in this document, please view it at [GitHub](https://github.com/MannLabs/alphaquant).

* [**About**](#about)
* [**License**](#license)
* [**Installation**](#installation)
  * [**One-click GUI**](#one-click-gui)
  * [**Developer installer**](#developer)
* [**Usage**](#usage)
  * [**GUI**](#gui)
  * [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
* [**Troubleshooting**](#troubleshooting)
* [**Citations**](#citations)
* [**How to contribute**](#how-to-contribute)
* [**Changelog**](#changelog)

---
## About
The standard approach for proteomics quantification is the calculation of point estimates that reflect the abundance of a particular protein. This approach usually neglects a large part of the quantitative information that is available, including the type, quality and reliability of the underlying, quantified peptides. AlphaQuant introduces a collection of novel Bioinformatics algorithms for increased accuracy and sensitivity of proteomics quantification. It is built on the foundation of the [MS-EmpiRe](https://doi.org/10.1074/mcp.RA119.001509) algorithm. 
Alphaquant is an open-source Python package of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/).

---
## License

AlphaQuant was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/) and is freely available with an [Apache License](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---
## Installation

AlphaQuant can be installed and used on all major operating systems (Windows, macOS and Linux).
There are currently two different types of installation possible:

* [**One-click GUI installer:--under construction--**](#one-click-gui) Choose this installation if you only want the GUI and/or keep things as simple as possible. Note that this version is quite outdated and does not contain many of the features. Update is in the making.
<!---
* [**Pip installer:**](#pip) Choose this installation if you want to use AlphaQuant as a Python package in an existing python 3.9 environment (e.g. a Jupyter notebook). If needed, the GUI and CLI can be installed with pip as well.
-->
* [**Developer installer:**](#developer) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/) and Python. This installation allows access to all available features of AlphaQuant and even allows to modify its source code directly. Generally, the developer version of AlphaQuant outperforms the precompiled versions which makes this the installation of choice for high-throughput experiments.

### One-click GUI --under construction--

The GUI of AlphaQuant is a completely stand-alone tool that requires no knowledge of Python or CLI tools. **Note that this version is quite outdated and does not contain many of the features. Update is in the making.** Click on one of the links below to download the latest release for:

* [**Windows**](https://github.com/MannLabs/alphaquant/releases/latest/download/alphaquant_gui_installer_windows.exe)
* [**macOS**](https://github.com/MannLabs/alphaquant/releases/latest/download/alphaquant_gui_installer_macos.pkg)
* [**Linux**](https://github.com/MannLabs/alphaquant/releases/latest/download/alphaquant_gui_installer_linux.deb)

Older releases remain available on the [release page](https://github.com/MannLabs/alphaquant/releases), but no backwards compatibility is guaranteed.

<!---
### Pip

AlphaQuant can be installed in an existing python 3.9 environment with a single `bash` command. *This `bash` command can also be run directly from within a Jupyter notebook by prepending it with a `!`*:

```bash
pip install alphaquant
```

Installing AlphaQuant like this avoids conflicts when integrating it in other tools, as this does not enforce strict versioning of dependancies. However, if new versions of dependancies are released, they are not guaranteed to be fully compatible with AlphaQuant. While this should only occur in rare cases where dependencies are not backwards compatible, you can always force AlphaQuant to use dependancy versions which are known to be compatible with:

```bash
pip install "alphaquant[stable]"
```

NOTE: You might need to run `pip install pip==21.0` before installing AlphaQuant like this. Also note the double quotes `"`.

For those who are really adventurous, it is also possible to directly install any branch (e.g. `@development`) with any extras (e.g. `#egg=alphaquant[stable,development-stable]`) from GitHub with e.g.

```bash
pip install "git+https://github.com/MannLabs/alphaquant.git@development#egg=alphaquant[stable,development-stable]"
```
-->
### Developer

AlphaQuant can also be installed in editable (i.e. developer) mode with a few `bash` commands. This allows to fully customize the software and even modify the source code to your specific needs. When an editable Python package is installed, its source code is stored in a transparent location of your choice. While optional, it is advised to first (create and) navigate to e.g. a general software folder:

```bash
mkdir ~/folder/where/to/install/software
cd ~/folder/where/to/install/software
```

***The following commands assume you do not perform any additional `cd` commands anymore***.

Next, download the AlphaQuant repository from GitHub either directly or with a `git` command. This creates a new AlphaQuant subfolder in your current directory.

```bash
git clone https://github.com/MannLabs/alphaquant.git
```

For any Python package, it is highly recommended to use a separate [conda virtual environment](https://docs.conda.io/en/latest/), as otherwise *dependancy conflicts can occur with already existing packages*.

```bash
conda create --name alphaquant python=3.9 -y
conda activate alphaquant
```

Finally, install AlphaQuant:

```bash
pip install -e .
```
By using the editable flag `-e`, you can make modifications to the [alphaquant source code](alphaquant) and these modifications will be directly reflected when running AlphaQuant. We currently recommend the stable 

Some details: By default this installs loose dependancies (no explicit versioning). It is also possible to install additional [development dependencies](requirements/requirements_development.txt), which allows to make use of more features (the call is then a bit more complex and could be e.g. `pip install -e "./alphaquant[stable,development-stable]"`).


---
## Usage

There are two ways to use AlphaQuant:

* [**GUI** --under construction--](#gui)
<!---* [**CLI**](#cli)-->
* [**Python and Jupyter Notebooks**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of AlphaQuant, it is often quite slow because some functions might still need compilation on your local operating system and architecture. Subsequent use should be a lot faster.

### GUI --under construction--

The GUI is currently accessible through the one-click GUI installer. Currently it only does the pairwise analysis. Further functionalities can be accessed through python and jupyter notebooks.
<!-- If the GUI was not installed through a one-click GUI installer, it can be activate with the following `bash` command:

```bash
alphaquant gui
```

Note that this needs to be prepended with a `!` when you want to run this from within a Jupyter notebook. When the command is run directly from the command-line, make sure you use the right environment (activate it with e.g. `conda activate alphaquant` or set an alias to the binary executable (can be obtained with `where alphaquant` or `which alphaquant`)).     - [ ] Fix GUI-->

<!---
### CLI

The CLI can be run with the following command (after activating the `conda` environment with `conda activate alphaquant` or if an alias was set to the AlphaQuant executable):

```bash
alphaquant -h
```

It is possible to get help about each function and their (required) parameters by using the `-h` flag.
-->

### Python and Jupyter notebooks

We have compiled a set of Jupyter notebooks together with some example data in the [example_nbs folder](example_nbs). There, you can use very simple calls in order to:
 * perform very sensitive differential expression analysis on a single condition, analyze and visualize proteoforms [here](example_nbs/differential_expression.ipynb)
 * analyze multiple condition together and inspect proteoform profiles [here](example_nbs/multi_condition_analysis.ipynb)
 * perform phosphosite and ptm mapping with subsequent differential expression analysis, as well as proteome normalization of phospho sites [here](example_nbs/differential_expression_PTM.ipynb)
 * combine the AlphaQuant proteoform analysis with deep learning on sequences in order to infer regulated phospho peptides from un-enriched standard proteome data [here](example_nbs/phospho_inference_analysis.ipynb)
 * visualize the tree structure of differential expression analysis [here](example_nbs/visualizing_tree_structure.ipynb)


## Preparing input files

**note: AlphaQuant is currently under development, mostly using DIA-NN and Spectronaut outputs. Other search engines and the generic format might cause problems.**
### Spectronaut

AlphaQuant takes a Spectronaut .tsv table as input. When exporting from Spectronaut, the correct columns need to be selected. These can be obtained by downloading one of the export schemes available below. We provide one export scheme for sprecursor quantification and one export scheme for fragment ion quantification. Fragment ion quantification shows slightly more accuracy, but the files are around 10 times larger.

An export scheme can then simply be loaded into Spectronaut as follows:

Go to the "Report" perspective in Spectronaut, click "Import Schema" and provide the file.

The data needs to be exported in the normal long format as .tsv file. Please double check that the schema is actually selected, sometimes Spectronaut (or at least older versions) lags when you select the schema. You should see that the preview changes when you click on it.



<a href="https://github.com/MannLabs/AlphaQuant/raw/master/alphaquant/config/spectronaut_tableconfig_precursor.rs" download>Download Spectronaut export scheme for precursor quantification</a>

<a href="https://github.com/MannLabs/AlphaQuant/raw/master/alphaquant/config/spectronaut_tableconfig_fragion.rs" download>Download Spectronaut export scheme for fragment ion quantification</a>

<a href="https://github.com/MannLabs/AlphaQuant/raw/master/alphaquant/config/spectronaut_tableconfig_ptm_fragion.rs" download>Download Spectronaut export scheme for fragment ion quantification WITH PTM </a>

The **samplemap.tsv** file must map to the **R.Label** column.


### DIA-NN
Provide the path to the DIANN "report.tsv" output table.  
The **samplemap.tsv** file must map the the **File.Name** column.

### MaxQuant
Provide the path to the MaxQuant "peptides.txt" output table or the MaxQuant evidence.txt output table.  
For "peptides.txt", the **samplemap.tsv** file must map the names of the columns starting with "Intensity ", but **without** the "Intensity ". For example "Intensity sample1.raw" "Intensity sample2.raw"-> "sample1.raw" "sample2.raw".  
For "evidence.txt, the **samplemap.tsv** file must map the **Experiment** column.

### FragPipe
Provide the path to the "combined_ion.tsv" output table.  
For "peptides.txt", the **samplemap.tsv** file must map the names of the columns ending with " Intensity", but **without** the " Intensity". For example "sample1 Intensity" "sample2 Intensity"-> "sample1" "sample2".



## Output tables

### results.tsv

* *condition_pair:* the names of the two conditions that are compared against each other (condition1 _VS_ condition2).  The log2 fold change is calculated as condition1 - condition2
* *p_value:* the uncorrected(!) p-value of the differential expression analysis. It tests the null hypothesis: 'no change between condition1 and condition2'. Lower values mean higher significance.
* *fdr:* the multiple testing corrected p-value with the Benjamini-Hochberg method
* *log2fc:* the estimated log 2 fold change.
* *number_of_ions:* number of raw datapoints used for protein intensity estimation.
* *quality_score:* a quantitative score indicating the quality of quantification. Higher scores mean higher quality.
* *summed_intensity:* the summed (non-log) intensities of all base ions


### proteoforms.tsv

* *protein:* protein or gene name
* *proteoform_id:* the protein name with a number at the end, indicating the nth proteoform. For EGFR_0 would be the reference proteoform of EGFR, EGFR_1 would indicated a second group of EGFR peptides that behave differently to EGFR_0. Many proteins will only have one reference proteoform.
* *cluster:* proteoform number
* *is_reference:*	TRUE if proteoform is the reference proteoform
* *peptides:* sequences of all peptides that map
* *quality_score:* alphaquant quality score between 0 and 1 (higher is better)
* *log2fc:* the estimated log2 fold change
* *fraction_of_peptides:* the fraction of peptides within the whole protein that belongs to the proteoform
* *fcdiff:* fold change difference relative to the reference proteoform

---
## Troubleshooting

In case of issues, check out the following:

* [Issues](https://github.com/MannLabs/alphaquant/issues): Try a few different search terms to find out if a similar problem has been encountered before
* [Discussions](https://github.com/MannLabs/alphaquant/discussions): Check if your problem or feature requests has been discussed before.

---
## Citations

Manuscript in preparation.

---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/alphaquant/stargazers) to boost our visibility! All direct contributions are also welcome. Feel free to post a new [issue](https://github.com/MannLabs/alphaquant/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/alphaquant/pulls) with a new branch. For an even more interactive participation, check out the [discussions](https://github.com/MannLabs/alphaquant/discussions) and the [the Contributors License Agreement](misc/CLA.md).

---
## Changelog

See the [HISTORY.md](HISTORY.md) for a full overview of the changes made in each version.
