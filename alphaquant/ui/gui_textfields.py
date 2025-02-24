from turtle import width
import panel.pane
import panel as pn
import os
import pathlib


class Paths():
    CONFIGS_PATH = os.path.join(pathlib.Path(__file__).parent.parent.absolute(), "..","config")
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
    run_pipeline_instruction = panel.pane.Markdown("""
#### **Run Pipeline**

Follow these steps to analyze your data:
1. Upload your proteomics data file
2. Set output folder
3. Map samples to conditions
4. Choose analysis mode
5. Select condition pairs (for pairwise comparison)
6. Click RUN_PIPELINE

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

Use the help icons (?) for detailed instructions.
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

Use the help icons (?) for detailed instructions.
""",
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 10))

    intro_text = panel.pane.Markdown("""
####

[Download example data](https://datashare.biochem.mpg.de/s/m1qR1hbz7lOyIzn/download)
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
**DIA-NN:**
Provide the path to the DIANN report.tsv output table.

**AlphaPept:**
Provide the path to the AlphaPept results_peptides.csv output table.

**MaxQuant:**
Provide the path to the MaxQuant peptides.txt output table.

**Spectronaut:**
Spectronaut exports tables based on user specification. You can download predefined configs below.
Go to the "Report" perspective in Spectronaut, click "Import Schema" and provide the file.
The data needs to be exported in long format as .tsv or .csv file.
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
        'file_input': """Supported file formats:
- DIA-NN: report.tsv
- AlphaPept: results_peptides.csv
- MaxQuant: peptides.txt
- Spectronaut: custom export (see table config downloads)""",

        'sample_mapping': """Two options available:
1. GUI Mapping: Fill in conditions for each sample in the interactive table
2. Manual Upload: Prepare a samplemap.tsv with 'sample' and 'condition' columns""",

        'analysis_mode': """Choose between:
- Pairwise Comparison: Compare specific condition pairs
- Median Condition Analysis: Compare each condition against the median of all conditions""",

        'filtering_options': """Available filtering modes:
- OR mode: ≥2 values in either condition (default)
- AND mode: ≥2 values in both conditions
- Custom: Specify values per condition

Note: Missing values handled by AlphaQuant counting statistics.""",

        'ptm_settings': """For Spectronaut PTM analysis:
1. Specify modification type (e.g., '[Phospho(STY)]')
2. Select organism for proteome mapping
Currently supports human and mouse."""
    }


class Cards():
    width = 530


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



