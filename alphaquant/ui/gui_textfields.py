from turtle import width
import panel.pane
import panel as pn
import os
import pathlib


class Paths():
    CONFIGS_PATH = os.path.join(pathlib.Path(__file__).parent.parent.absolute(), "config")
    spectronaut_fragion_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_fragion.rs")
    spectronaut_precursor_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_precursor.rs")
    spectronaut_ptm_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_ptm_fragion.rs")


class ButtonConfiguration():
    width = 530


class DownloadSchemes():
    spectronaut_fragion = pn.widgets.FileDownload(
        file=Paths.spectronaut_fragion_path,
        filename="spectronaut_tableconfig_fragion.rs",
        label="spectronaut_tableconfig_fragion.rs",
        button_type="light",
    )

    spectronaut_precursor = pn.widgets.FileDownload(
        file=Paths.spectronaut_precursor_path,
        filename="spectronaut_tableconfig_precursor.rs",
        label="spectronaut_tableconfig_precursor.rs",
        button_type="light",
    )

    spectronaut_ptm = pn.widgets.FileDownload(
        file=Paths.spectronaut_ptm_path,
        filename="spectronaut_tableconfig_ptm_fragion.rs",
        label="spectronaut_tableconfig_ptm_fragion.rs",
        button_type="light",
    )


class Descriptions():
    who_should_use = panel.pane.Markdown("""
#### **Who Should Use AlphaQuant?**

AlphaQuant is designed for proteomics researchers analyzing DDA or DIA experiments with multiple conditions (e.g., control vs. treatment, time-series, or multi-condition studies).

Key features:
- **All-in-one Statistical Analysis**: Comprehensive statistical analysis from normalization to multiple testing correction, with volcano plots and other visualizations
- **Sensitive Detection of Changes**: Captures subtle patterns using Fragment and MS1-level analysis with intensity-dependent counting statistics
- **Proteoform Analysis**: Automatic clustering of peptides with similar quantitative behavior to infer regulated proteoforms
- **Support for Major Search Engines**: Direct support for DIA-NN, Spectronaut, AlphaDIA, MaxQuant, FragPipe, and AlphaPept
""",
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    run_pipeline_instruction = panel.pane.Markdown("""
#### **Run Pipeline**

Follow these steps to analyze your data:
1. Upload your proteomics data file
2. Set output folder
3. Map samples: Either use GUI table or create samplemap.tsv
4. Select mode: "Pairwise Comparison" or "Median Condition Analysis"
5. Choose conditions to compare (for pairwise mode)
6. Click "RUN PIPELINE" to execute

For detailed instructions, use the help icons (?) next to each control.
""",
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    basic_plots_instruction = panel.pane.Markdown("""
#### **Basic Plots**

1. Select results directory
2. Upload sample mapping file
3. Choose visualization options

""",
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    proteoform_plots_instruction = panel.pane.Markdown("""
#### **Proteoform Plots**

1. Select results directory
2. Upload sample mapping file
3. Choose organism
4. Select condition pair
5. Pick protein to visualize

""",
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    intro_text = panel.pane.Markdown("""
####

<a href="https://datashare.biochem.mpg.de/s/j98TnBzPJWLrtnC/download" target="_blank">Download example data</a>
""",
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    single_comparison_instruction = panel.pane.Markdown("""
        Here you can visualize comparisons of two conditions as a volcano plot. You can click or search proteins of interest and detail plots of the quantification will be shown.
        The displayed data is stored as text files in the output folder you specified.
        """,
        width=830,
        align='start',
        margin=(0, 80, 0, 10))

    table_instructions = pn.Column(
        pn.pane.Markdown("""

**AlphaDIA:**
Provide the path to "precursors.tsv", or "fragment_precursorfiltered.matrix.parquet". The samplemap.tsv file must map to the run column of the precursors.tsv.

**DIA-NN:**
Provide the path to the DIANN "report.tsv" output table. The samplemap.tsv file must map the Run column.

**AlphaPept:**
Provide the path to the AlphaPept results_peptides.csv output table.

**MaxQuant:**
Provide the path to the MaxQuant "peptides.txt" or "evidence.txt" output table.
- For "peptides.txt": Map samplemap.tsv to column names starting with "Intensity " (without "Intensity ")
  Example: "Intensity sample1.raw" → "sample1.raw"
- For "evidence.txt": Map samplemap.tsv to the Experiment column

**FragPipe:**
Provide the path to the "combined_ion.tsv" output table. Map samplemap.tsv to column names ending with " Intensity" (without " Intensity")
Example: "sample1 Intensity" → "sample1"

**Spectronaut:**
AlphaQuant takes a Spectronaut .tsv table as input. Export schemes are available below for:
- Precursor quantification
- Fragment ion quantification (more accurate but ~10x larger files)
- Fragment ion quantification with PTM

To use an export scheme:
1. Go to the "Report" perspective in Spectronaut
2. Click "Import Schema" and provide the file
3. Export data in long format as .tsv file
Note: Verify schema selection as Spectronaut may lag when selecting

The samplemap.tsv file must map to the R.Label column.
"""),
        pn.Row(
            DownloadSchemes.spectronaut_fragion,
            pn.pane.Markdown("Most detailed report, good for analyses where you need high statistical power (e.g. small fold changes, or few peptides)")
        ),
        pn.Row(
            DownloadSchemes.spectronaut_precursor,
            pn.pane.Markdown("About 10x less data heavy, good for analyses with clear regulation happening")
        ),
        pn.Row(
            DownloadSchemes.spectronaut_ptm,
            pn.pane.Markdown("For PTM analyses")
        ),
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 20)
    )

    # Add tooltips/help text for each control
    tooltips = {
        'file_input': """Provide the filepath to your proteomic datasets analyzed by DIA-NN, Spectronaut, AlphaDIA, AlphaPept, MaxQuant or FragPipe (details on which tables are needed are given in the "Instructions" tab)""",

        'sample_mapping': """Map the experiment names (i.e. the names of the MS runs, such as sample1_control_23_2025.raw) to the condition names (e.g. "control", "treatment").

1. You do not have a samplemap.tsv yet:
   - Provide the filepath to your proteomics dataset
   - Click "Generate Samplemap Template"
   - Experiment names will be displayed in an interactive table
   - Fill in the condition name for each sample
   - Alternatively, a file named samplemap_template.tsv is saved in your results directory. Fill it out and proceed with step 2.

2. You already have a samplemap.tsv
   - Set the switch to "Upload sample to condition file"
   - Upload samplemap.tsv file
""",

        'analysis_mode': """Choose between:
Decide the analysis mode. In most cases, this will be the "Pairwise Comparison" mode (e.g. treatment1 vs. control, treatment2 vs. control).
There is also a more global analysis, "Median Condition Analysis", where each condition will be compared against the median of all conditions.
This allows direct comparability of each condition.""",

        'filtering_options': """Available filtering modes:
- OR mode: ≥2 values in either condition (default)
- AND mode: ≥2 values in both conditions
- Custom: Specify values per condition

Note: Missing values are handled by AlphaQuant counting statistics.
""",

        'ptm_settings': """For Spectronaut PTM analysis:
1. Modification Type: Specify exactly as it appears in Spectronaut modified sequence
   Example: '[Phospho(STY)]' for phosphorylation
2. Organism: Select proteome for site mapping (human/mouse available at the moment)

Note: Requires Spectronaut table with correct PTM columns (see table instructions)."""
    }


class Cards():
    width = 530

    who_should_use = pn.Card(
        Descriptions.who_should_use,
        header='Who Should Use AlphaQuant?',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    table_instructions = pn.Card(
        Descriptions.table_instructions,
        header='Table instructions for different search engines',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    run_pipeline = pn.Card(
        Descriptions.run_pipeline_instruction,
        header='Run Pipeline',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    basic_plots = pn.Card(
        Descriptions.basic_plots_instruction,
        header='Basic Plots',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    proteoform_plots = pn.Card(
        Descriptions.proteoform_plots_instruction,
        header='Proteoform Plots',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )



