from turtle import width
import panel.pane
import panel as pn
import os
import pathlib


class Paths():
    CONFIGS_PATH = os.path.join(pathlib.Path(__file__).parent.absolute(), "configs")
    spectronaut_fragion_path = os.path.join(CONFIGS_PATH, "spectronaut_tableconfig_fragion.rs")

class ButtonConfiguration():
    width = 530


class Descriptions():


    project_instruction = panel.pane.Markdown("""
#### How to use AlphaQuant:

[Download example data](https://github.com/MannLabs/alphaquant/raw/main/example_nbs/data.zip)

You can also find example input data here "https://github.com/MannLabs/alphaquant/tree/main/example_nbs/data".

* **Run Pipeline**
    The run pipeline tab allows you to run differential expression analysis on your data and will write out results tables for you. To run the pipeline, you need to:
    1. Provide the filepath to your proteomic datasets analyzed by DIA-NN, Spectronaut, AlphaDIA, AlphaPept, MaxQuant or FragPipe (detailed instructions on which tables are needed are given below).
    2. Specify a folder where results are saved
    3. Map the experiment names (i.e. the names of the MS runs, such as sample1_control_23_2025.raw) to the condition names (e.g. "control" and "treatment"). You have two options here:
        1) Do the sample mapping here in the GUI. Once you provide the filepath to your proteomics dataset, the experiment names will be displayed in a small interactive table and you can fill in the condition name for each sample.
        2) You can prepare a samplemap.tsv yourself, e.g. with Excel or any text editor. The column names are *sample* and *condition* and they are separated by tabs. The names of the MS runs are extracted from your input file and they differ for the different search engines. You can check which column contains the experiments in the table instructions below (e.g. in the DIA-NN table, the column is called 'Run').
    4. Now, you need to decide the analysis mode. In most cases, this will be the "Pairwise Comparison" mode (e.g. treatment1 vs. control, treatment2 vs. control). There is also a more global analysis, "Median Condition Analysis", where each condition will be compared against the median of all conditions. This allows direct comparability of each condition.
    5. If you chose "Pairwise Comparison", you need to specify which conditions you want to compare. To select the condition pairs you want, highlight them in the box on the left and use the arrow to move them to the right.
    6. Now you can click the _RUN_PIPELINE_ button and the analysis will be carried out. You can track the progress via the GUI.

    **Basic Settings for Run Pipeline:**
    S1. Filtering options. It often happens that one protein is detected in several experiments in one condition and in very few or none in the other condition. We need to specify how to handle these cases. Per default, all proteins will be included that have at least 2 datapoints in one of the two conditions (`min. valid values in condition1 OR condition2`). If the datapoints are missing in the other condition, the AlphaQuant counting statistics module is used. There is also the option to filter only for proteins that have at least two datapoints in both conditons (`min. valid values in condition1 AND condition2`), or you can individually specify the number of datapoints in each condition (carefull if you have multiple condition pairs!)

    **PTM Settings for Run Pipeline:**
    1. Modification Type. If you use Spectronaut tables, you can run a PTM analysis including site mapping. You need to export the Spectronaut table with the correct columns (see instructions below). Then you need to specify the type of modification you want to perform site mapping on, in the way it appears in the Spectronaut modified sequence. In the case of phospho, for example, this is `[Phospho(STY)]`
    2. Organism. AlphaQuant needs to know which proteome .fastas to use in order to perform site mapping. Currently you can chose between human and mouse. If you need broader support, please reach out.

* **Basic Plots**
    The basic plots will show you a volcano plot for every condition pair. Additionally you can visualize the peptides underlying each protein individually. For this, you need to:
    1. Provide the filepath to the results directory
    2. Upload the sample mapping file. If you have created the samplemap from the GUI, you can find the samplemap.tsv in the results directory.

* **Proteoform Plots**
    In this tab you can re-create peptide resolved plots mapped to the protein sequence, as described in the AlphaQuant paper. For this, you need to
    1. Provide the filepath to the results directory
    2. Provide the fielpath to the sample mapping file. If you have created the samplemap from the GUI, you can find the samplemap.tsv in the results directory.
    3. Specify the organism. AlphaQuant needs to know which proteome .fastas to use in order to create the proteoform plots. Currently you can chose between human, mouse and yeast. If you need broader support, please reach out.
    4. Select the condition pair, where you want to investigate proteoform candidates
    5. Select which protein to visualize. A table with all proteoform candidates will appear and you can select by either clicking on the row in the table, or by typing the protein name in the field below.
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

    alphapept = pn.pane.Markdown(
        """
        Provide the path to the AlphaPept results_peptides.csv output table.

        """,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 20)
    )

    spectronaut = pn.pane.Markdown(
        """
        To get the most out of the Spectronaut data, AlphaQuant utilizes more than 30 different columns.
        These can be obtained by downloading the export scheme "spectronaut_tableconfig_fragion.rs",
        which can then simply be loaded into Spectronaut as follows:

        Go to the "Report" perspective in Spectronaut, click "Import Schema" and provide the file.

        The data needs to be exported in the **normal long** format as .tsv or .csv file.

        """,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 0, 20)
    )

    diann = pn.pane.Markdown(
            """
            Provide the path to the DIANN report.tsv output table.
            """,
            width=ButtonConfiguration.width,
            align='start',
            margin=(0, 80, 0, 20)
        )

    maxquant = pn.pane.Markdown(
            """
            Provide the path to the MaxQuant peptides.txt output table.
            """,
            width=ButtonConfiguration.width,
            align='start',
            margin=(0, 80, 0, 20)
        )




class DownloadSchemes():

    spectronaut = pn.widgets.FileDownload(
    file=Paths.spectronaut_fragion_path,
    filename="spectronaut_tableconfig_fragion.rs",
    button_type='default',
    auto=True,
    css_classes=['button_options'],
)


class Cards():
    width = 530

    alphapept = pn.Card(
        Descriptions.alphapept,
        header='AlphaPept instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )


    spectronaut = pn.Card(
        Descriptions.spectronaut,
        DownloadSchemes.spectronaut,
        header='Spectronaut instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(0, 80, 5, 10),
        css_classes=['spectronaut_instr']
    )
    diann = pn.Card(
        Descriptions.diann,
        header='DIANN instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )

    maxquant = pn.Card(
        Descriptions.maxquant,
        header='MaxQuant instructions',
        collapsed=True,
        width=ButtonConfiguration.width,
        align='start',
        margin=(20, 0, 20, 0),
        css_classes=['spectronaut_instr']
    )



