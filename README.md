![GitHub Release](https://img.shields.io/github/v/release/mannlabs/alphaquant?logoColor=green&color=brightgreen)
![Versions](https://img.shields.io/badge/python-3.10_%7C_3.11_%7C_3.12-brightgreen)
![License](https://img.shields.io/badge/License-Apache-brightgreen)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mannlabs/alphaquant/e2e_tests_quick_multiple_platforms.yml?branch=main&label=E2E%20Tests)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mannlabs/alphaquant/install_and_unit_tests.yml?branch=main&label=Unit%20Tests)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mannlabs/alphaquant/publish_on_pypi.yml?branch=main&label=Deploy%20PyPi)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MannLabs/alphaquant/main?urlpath=%2Fdoc%2Ftree%2Fexample_nbs%2Fdifferential_expression.ipynb)

<img src="release/images/alphaquant_gui.jpg" alt="preview" width="800"/>

==> [Run it on a Jupyter Notebook right now in your browser!](https://mybinder.org/v2/gh/MannLabs/alphaquant/main?urlpath=%2Fdoc%2Ftree%2Fexample_nbs%2Fdifferential_expression.ipynb) No login or installation required.
<== 

# AlphaQuant
AlphaQuant is an innovative open-source Python package for proteomics data analysis. It implements tree-based quantification - a hierarchical approach to organize and analyze quantitative data across multiple levels - from fragments and MS1 isotopes through charge states, modifications, peptides, and genes.

It is part of the AlphaPept ecosystem from the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/).


## Who Should Use AlphaQuant?

AlphaQuant is designed for proteomics researchers analyzing DDA or DIA experiments with multiple conditions (e.g., control vs. treatment, time-series, or multi-condition studies). If your goal is to compare and interpret quantitative proteomics data systematically, AlphaQuant provides:

- **All-in-one Statistical Analysis**: AlphaQuant delivers comprehensive statistical analysis of your differential experiments, performing all critical steps from normalization to multiple testing correction in one go, with results visualized through volcano plots and other informative displays.
- **Sensitive Detection of Changes**: AlphaQuant excels at capturing subtle patterns and handling missing values to ensure important biological signals are not overlooked. This is achieved by using Fragment and MS1-level analysis as well as intensity-dependent counting statistics.
- **Proteoform Analysis**: AlphaQuant automatically performs clustering of peptides with similar quantitative behavior to infer regulated proteoforms.
- **Support for Major Search Engines**: Direct support for all major search engines in DDA and DIA workflows (DIA-NN, Spectronaut, AlphaDIA, MaxQuant, FragPipe, AlphaPept) - just use their standard output files



## Table of Contents


* [**Installation**](#installation)
  * [**One-click GUI**](#one-click-gui-installation)
  * [**Pip**](#pip)
  * [**Developer installation**](#developer-installation)
  * [**Docker**](#docker)
* [**Usage**](#usage)
  * [**GUI**](#gui)
  * [**Python and jupyter notebooks**](#python-and-jupyter-notebooks)
* [**Troubleshooting**](#troubleshooting)
* [**Citations**](#citations)
* [**How to contribute**](#how-to-contribute)
* [**License**](#license)
* [**Changelog**](#changelog)


---
## Installation

AlphaQuant can be installed and used on all major operating systems (Windows, macOS and Linux).
There are currently two different types of installation possible:

* [**One-click GUI installer**](#one-click-gui-installation) Choose this installation if you only want the GUI and/or keep things as simple as possible.

* [**Pip installer:**](#pip) Choose this installation if you want to use AlphaQuant as a Python package in an existing python 3.11 environment (e.g. a Jupyter notebook). If needed, the GUI can be installed with pip as well.

* [**Developer installation:**](#developer-installation) Choose this installation if you are familiar with CLI tools, [conda](https://docs.conda.io/en/latest/) and Python. This installation allows access to all available features of AlphaQuant and even allows to modify its source code directly. Generally, the developer version of AlphaQuant outperforms the precompiled versions which makes this the installation of choice for high-throughput experiments.

* [**Docker**](#docker) Choose this installation if you want to use AlphaQuant without any installation to your system.
 
### One-click GUI installation
Currently available for **MacOS**, **Windows**.
You can download the latest release of alphaquant [here](https://github.com/Mannlabs/alphaquant/releases/latest).

* **Windows:** Download the latest `alphaquant-X.Y.Z-windows-amd64.exe` build and double click it to install. If you receive a warning during installation click *Run anyway*.
* **MacOS:** Download the latest build suitable for your chip architecture
(can be looked up by clicking on the Apple Symbol > *About this Mac* > *Chip* ("M1", "M2", "M3" -> `arm64`, "Intel" -> `x64`),
`alphaquant-X.Y.Z-macos-darwin-arm64.pkg` or `alphaquant-X.Y.Z-macos-darwin-x64.pkg`. Open the parent folder of the downloaded file in Finder,
right-click and select *open*. If you receive a warning during installation click *Open*. If you want to use `.raw` files on Thermo instruments alphaRaw is required, which depends on Mono. A detailed guide to installing alphaRaw with mono can be found [here](https://github.com/MannLabs/alpharaw#installation).
* **Linux:** Installers are provided, but undergo only limited testing: `alphaquant-X.Y.Z-linux-x64.deb` build and install it via `dpkg -i alphaquant-X.Y.Z-linux-x64.deb`. In case of issues, follow the steps for the
[developer installation](docs/installation.md#developer-installation) in order to use the GUI.



### Pip

AlphaQuant can be installed in an existing python 3.11 environment with

```bash
pip install alphaquant
```

Installing AlphaQuant like this avoids conflicts when integrating it in other tools, as this does not enforce strict versioning of dependancies. However, if new versions of dependancies are released, they are not guaranteed to be fully compatible with AlphaQuant. While this should only occur in rare cases where dependencies are not backwards compatible, you can always force AlphaQuant to use dependancy versions which are known to be compatible with:

```bash
pip install "alphaquant[stable]"
```

if you want to add the GUI to your environment, you can install it with the following command:

```bash
pip install "alphaquant[stable,gui-stable]"
```

For those who are really adventurous, it is also possible to directly install any branch (e.g. `@development`) with any extras (e.g. `#egg=alphaquant[stable,development-stable]`) from GitHub with e.g.

```bash
pip install "git+https://github.com/MannLabs/alphaquant.git@development#egg=alphaquant[stable,development-stable]"
```

### Developer installation

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
conda create --name alphaquant python=3.11 -y
conda activate alphaquant
```

Finally, install AlphaQuant:

```bash
pip install -e .
```
By using the editable flag `-e`, you can make modifications to the [alphaquant source code](alphaquant) and these modifications will be directly reflected when running AlphaQuant. We currently recommend the stable

Some details: By default this installs loose dependancies (no explicit versioning). It is also possible to install additional [development dependencies](requirements/requirements_development.txt), which allows to make use of more features (the call is then a bit more complex and could be e.g. `pip install -e "./alphaquant[stable,development-stable]"`).

### Docker
The containerized version can be used to run AlphaQuant without any installation to your system.

#### 1. Setting up Docker
Install the latest version of docker (https://docs.docker.com/engine/install/).

#### 2. Prepare folder structure
Set up your data to match the expected folder structure:Create a folder and store its name in a variable, 
e.g. `DATA_FOLDER=/home/username/data; mkdir -p $DATA_FOLDER`

#### 3. Start the container
```bash
docker run -v $DATA_FOLDER:/app/data -p 41215:41215 mannlabs/alphaquant:latest
```
After initial download of the container, alphaquant will start running immediately,
and can be accessed under [localhost:41215](localhost:41215). 

Note: in the app, the local `$DATA_FOLDER` needs to be referred to as "`/app/data`".

#### Alternatively: Build the image yourself
If you want to build the image yourself, you can do so by
```bash
docker build -t alphaquant .
```
and run it with
```bash
docker run -p 41215:41215 -v $DATA_FOLDER:/app/data -t alphaquant
```

---
## Usage

There are two ways to use AlphaQuant:

* [**GUI**](#gui)
<!---* [**CLI**](#cli)-->
* [**Python and Jupyter Notebooks**](#python-and-jupyter-notebooks)

NOTE: The first time you use a fresh installation of AlphaQuant, it is often quite slow because some functions might still need compilation on your local operating system and architecture. Subsequent use should be a lot faster.

### GUI usage

If the GUI was not installed through a one-click GUI installer, it can be activate with the following `bash` command:

```bash
alphaquant gui
```


### Python and Jupyter notebooks

Quickstart:

```python
import alphaquant.run_pipeline as aq_pipeline

aq_pipeline.run_pipeline(input_file=INPUT_FILE, samplemap_file=SAMPLEMAP_FILE, results_dir=RESULTS_DIRECTORY)
```

For more detailed examples and advanced use cases, we provide several Jupyter notebooks with example data in the [example_nbs folder](example_nbs): There, you can use very simple calls in order to:
 * perform very sensitive differential expression analysis on a single condition, analyze and visualize proteoforms [here](example_nbs/differential_expression.ipynb)
 * analyze multiple condition together and inspect proteoform profiles [here](example_nbs/multi_condition_analysis.ipynb)
 * perform phosphosite and ptm mapping with subsequent differential expression analysis, as well as proteome normalization of phospho sites [here](example_nbs/differential_expression_PTM.ipynb)



## Preparing input files

### The samplemap.tsv file
The samplemap.tsv is a **tab-separated** file that is always required. In the GUI, you can create it during the setup process. The samplemap.tsv maps the experiment names (i.e. the individual runs) to the condition names (e.g. "control" and "treatment"). The column names are **sample** and **condition**. A typical example of a samplemap.tsv file is:

```
sample	condition
run1	control
run2	control
run3	control
run4	treatment
run5	treatment
run6	treatment
```

### DIA-NN
Provide the path to the DIANN "report.tsv" output table.
The **samplemap.tsv** file must map the the **Run** column.

### AlphaDIA
Provide the path to "precursors.tsv", or "fragment_precursorfiltered.matrix.parquet"
The **samplemap.tsv** file must map to the **run** column.

### Spectronaut
AlphaQuant takes a Spectronaut .tsv table as input. When exporting from Spectronaut, the correct columns need to be selected. These can be obtained by downloading one of the export schemes available below. We provide one export scheme for sprecursor quantification and one export scheme for fragment ion quantification. Fragment ion quantification shows slightly more accuracy, but the files are around 10 times larger.

An export scheme can then simply be loaded into Spectronaut as follows:

Go to the "Report" perspective in Spectronaut, click "Import Schema" and provide the file.

The data needs to be exported in the normal long format as .tsv file. Please double check that the schema is actually selected, sometimes Spectronaut (or at least older versions) lags when you select the schema. You should see that the preview changes when you click on it.



<a href="https://github.com/MannLabs/AlphaQuant/raw/main/alphaquant/config/spectronaut_tableconfig_precursor.rs" download="spectronaut_tableconfig_precursor.rs">Download Spectronaut export scheme for precursor quantification</a>

<a href="https://github.com/MannLabs/AlphaQuant/raw/main/alphaquant/config/spectronaut_tableconfig_fragion.rs" download="spectronaut_tableconfig_fragion.rs">Download Spectronaut export scheme for fragment ion quantification</a>

<a href="https://github.com/MannLabs/AlphaQuant/raw/main/alphaquant/config/spectronaut_tableconfig_ptm_fragion.rs" download="spectronaut_tableconfig_ptm_fragion.rs">Download Spectronaut export scheme for fragment ion quantification WITH PTM</a>

The **samplemap.tsv** file must map to the **R.Label** column.

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
* *fraction_of_peptides:* the fraction of peptides within the whole protein that belongs to the proteoform_id
* *fcdiff:* fold change difference relative to the reference proteoform

---
## Troubleshooting

In case of issues, check out the following:

* [Issues](https://github.com/MannLabs/alphaquant/issues): Try a few different search terms to find out if a similar problem has been encountered before
* [Discussions](https://github.com/MannLabs/alphaquant/discussions): Check if your problem or feature requests has been discussed before.

---
## Citations

A manuscript has been submitted to bioRxiv:
> **Tree-based quantification infers proteoform regulation in bottom-up proteomics data**
> Constantin Ammar, Marvin Thielert, Caroline A M Weiss, Edwin H Rodriguez, Maximilian T Strauss, Florian A Rosenberger, Wen-Feng Zeng, Matthias Mann
> bioRxiv 2025.03.06.641844; doi: https://doi.org/10.1101/2025.03.06.641844 

---
## How to contribute

If you like this software, you can give us a [star](https://github.com/MannLabs/alphaquant/stargazers) to boost our visibility! All direct contributions are also welcome. Feel free to post a new [issue](https://github.com/MannLabs/alphaquant/issues) or clone the repository and create a [pull request](https://github.com/MannLabs/alphaquant/pulls) with a new branch. For an even more interactive participation, check out the [discussions](https://github.com/MannLabs/alphaquant/discussions) and the [the Contributors License Agreement](misc/CLA.md).

---

## License

AlphaQuant was developed by the [Mann Labs at the Max Planck Institute of Biochemistry](https://www.biochem.mpg.de/mann) and the [University of Copenhagen](https://www.cpr.ku.dk/research/proteomics/mann/) and is freely available with an [Apache License](LICENSE.txt). External Python packages (available in the [requirements](requirements) folder) have their own licenses, which can be consulted on their respective websites.

---
## Changelog

See the [HISTORY.md](HISTORY.md) for a full overview of the changes made in each version.
