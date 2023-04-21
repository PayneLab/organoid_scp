import os
import pandas as pd

def load_protein_table(
    source: str,
) -> pd.DataFrame:
    """
    Load the protein abundances from a protein table with the original
    sample names.

    Parameters
    ----------
    source : {"pd", "mm"}
        Which tool you want the protein abundance data from. "pd" for
        Proteome Discoverer, "mm" for MetaMorpheus.

    Returns
    -------
    pd.DataFrame
        The formatted protein abundance table.
    """

    if source == "mm":

        df = _load_original_table("mm", "AllQuantifiedProteinGroups") # Load protein table
        ed = _load_original_table("mm", "ExperimentalDesign") # Load experimental design table

        df = df.set_index("Protein Accession") # Index by protein accession
        df = df[df.columns[df.columns.str.startswith("Intensity_")]] # Select just the protein abundance columns

        df.columns = df.columns.to_series().str.split("_", n=1, expand=True)[1].values # Cut off "Intensity_" from column names
        df.columns.name = "sample"
        df.index.name = "protein"

        df = df.\
        transpose().\
        reset_index(drop=False)

        # Cut the ".raw" off the end of the filenames, and generate a "sample" column for joining with protein table
        ed = ed.assign(
            FileName=ed["FileName"].str.rsplit(".", n=1, expand=True)[0],
            sample=ed["Condition"] + "_" + ed["Biorep"].astype(str),
        )

        # Mark the "bad" and "np" (no protein) samples
        ed = ed.assign(
            bad=ed["FileName"].str.endswith("bad"),
            np=ed["FileName"].str.endswith("np"),
        )

        # Select only the needed columns
        ed = ed[["sample", "bad", "np"]]

        # Join into protein table
        df = ed.merge(
            right=df,
            on="sample",
            how="outer",
        )

    elif source == "pd":
        # Has both scaled abundance and found/not found columns--extract out separately
        # Check Contaminant=True?
        # Also has GO annotations and Reactome and WikiPathways pathways
        # Are all of these Master proteins for protein groups? And if so, is this table equivalent to groups tables?
        df = oscutils._load_original_table("pd", "Proteins")

    return df

def _load_original_table(
    source: str,
    name: str,
) -> pd.DataFrame:
    """
    Load one of the data tables without any additional formatting.

    Parameters
    ----------
    source : {"pd", "mm"}
        Which tool your desired table was generated by. "pd" for Proteome 
        Discoverer, "mm" for MetaMorpheus.
    name : {
        "AllQuantifiedProteinGroups",
        "AnnotationProteinGroups",
        "ExperimentalDesign",
        "PathwayProteinGroups",
        "ProteinGroups",
        "Proteins",
        "StudyInformation",
    }
        The name of the table you want. "AllQuantifiedProteinGroups" is the
        only option for MetaMorpheus; the rest are from Proteome
        Discoverer.

    Returns
    -------
    pd.DataFrame
        The requested table.
    """
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(path_here, "data")

    if source == "mm":
        data_path = os.path.join(data_path, "metamorpheus", "Task1-SearchTask")
        if name in [
            "AllQuantifiedProteinGroups",
            "ExperimentalDesign",
        ]:
            file_path = os.path.join(data_path, f"{name}.tsv")
        else:
            raise ValueError(f"Invalid argument '{name}' for name parameter with '{source}' as source.")
    elif source == "pd":
        data_path = os.path.join(data_path, "proteome_discoverer")
        if name in [
            "AnnotationProteinGroups", # Coverage numbers and percentages for GO annotations
            "PathwayProteinGroups", # Coverage numbers and percentages for Reactome and WikiPathways pathways
            "ProteinGroups", # Coverage for protein groups
            "Proteins",
            "StudyInformation",
        ]:
            file_path = os.path.join(data_path, f"Payne_Organoids_{name}.txt")
        else:
            raise ValueError(f"Invalid argument '{name}' for name parameter with '{source}' as source.")
    else:
        raise ValueError(f"Invalid argument '{source}' for source parameter.")

    if source == "mm" and name == "AllQuantifiedProteinGroups":
        # All the lines in the tsv for this table have an extra tab
        # character at the end of every line except for the first one,
        # which causes the headers to be read in offset by one column.
        # This line addresses that.
        df = pd.read_csv(file_path, sep="\t", usecols=range(136))
    else:
        df = pd.read_csv(file_path, sep="\t")

    return df
