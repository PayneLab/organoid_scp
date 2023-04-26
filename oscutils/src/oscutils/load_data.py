import numpy as np
import os
import pandas as pd

def load_protein_table(
    source: str,
    quant_or_found: str,
    clean: bool = True,
) -> pd.DataFrame:
    """
    Load the protein abundances from a protein table with the original
    sample names.

    Parameters
    ----------
    source : {"pd", "mm"}
        Which tool you want the protein abundance data from. "pd" for
        Proteome Discoverer, "mm" for MetaMorpheus.
    quant_or_found : {"quant", "found"}
        Whether to return "proteins quantified" or just "proteins found"
        data.
    clean : bool, default True
        Whether to drop contaminated and no protein samples from the table.

    Returns
    -------
    pd.DataFrame
        The formatted protein abundance table.
    """
    if source not in ["mm", "pd"]:
        raise ValueError(f"Invalid argument '{source}' for source parameter.")
    if quant_or_found not in ["quant", "found"]:
        raise ValueError(f"Invalid argument '{quant_or_found}' for quant_or_found parameter.")

    if source == "mm":

        if quant_or_found == "found":
            raise ValueError("'found' data not available from MetaMorpheus.")

        df = _load_original_table("mm", "AllQuantifiedProteinGroups") # Load protein table
        ed = _load_original_table("mm", "ExperimentalDesign") # Load experimental design table

        df = df.set_index("Protein Accession") # Index by protein accession
        df = df[df.columns[df.columns.str.startswith("Intensity_")]] # Select just the protein abundance columns

        df.columns = df.columns.to_series().str.split("_", n=1, expand=True)[1].values # Cut off "Intensity_" from column names

        # Zero-pad the sample numbers
        cols_split = df.columns.to_series().str.rsplit("_", n=1, expand=True)
        df.columns = cols_split[0] + "_" + cols_split[1].map("{:0>2}".format)

        df.columns.name = "sample"
        df.index.name = "protein"

        df = df.\
        transpose().\
        reset_index(drop=False)

        # Cut the ".raw" off the end of the filenames, and generate a "sample" column for joining with protein table
        ed = ed.assign(
            FileName=ed["FileName"].str.rsplit(".", n=1, expand=True)[0],
            sample=ed["Condition"] + "_" + ed["Biorep"].map("{:0>2}".format),
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
        df = _load_original_table("pd", "Proteins") # Load protein table
        si = _load_original_table("pd", "StudyInformation") # Load study information table
        mm_ed = _load_original_table("mm", "ExperimentalDesign") # Load MetaMorpheus experimental design table so we can get same sample numbers

        df = df.set_index("Accession") # Index by protein accession

        if quant_or_found == "quant":
            df = df[df.columns[df.columns.str.startswith("Abundances Scaled F")]]
            df.columns = df.columns.to_series().str.extract(r"Abundances Scaled (F\d+) Sample")[0].values
        elif quant_or_found == "found":
            df = df[df.columns[df.columns.str.startswith("Found in Sample F")]]
            df.columns = df.columns.to_series().str.extract(r"Found in Sample (F\d+) Sample")[0].values
            
        df.columns.name = "File ID"
        df.index.name = "protein"

        df = df.\
        transpose().\
        reset_index(drop=False)

        si = si.assign(FileName=si["File Name"].str.rsplit(".", n=1, expand=True)[0].str.rsplit("\\", n=1, expand=True)[1])
        si = si[["File ID", "FileName"]]

        # Cut the ".raw" off the end of the filenames, and generate a "sample" column for joining with protein table
        mm_ed = mm_ed.assign(
            FileName=mm_ed["FileName"].str.rsplit(".", n=1, expand=True)[0],
            sample=mm_ed["Condition"] + "_" + mm_ed["Biorep"].map("{:0>2}".format),
        )

        # Mark the "bad" and "np" (no protein) samples
        mm_ed = mm_ed.assign(
            bad=mm_ed["FileName"].str.endswith("bad"),
            np=mm_ed["FileName"].str.endswith("np"),
        )

        # Select only the needed columns
        mm_ed = mm_ed[["sample", "FileName", "bad", "np"]]

        si = si.merge(
            right=mm_ed,
            on="FileName",
            how="outer",
        )

        si = si[["File ID", "sample", "bad", "np"]]

        df = si.\
        merge(
            right=df,
            on="File ID",
            how="outer",
        ).\
        drop(columns="File ID")

    if clean:
        df = df[~df["bad"] & ~df["np"] & ~df["sample"].str.startswith("qc_")]
        df = df.drop(columns=["bad", "np"])

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

def _get_proteins_found_count(
    source: str,
    quant_or_found: str,
    clean: bool,
) -> pd.DataFrame:
 
    if source == "both":

        df_mm = _get_proteins_found_count(source="mm", quant_or_found=quant_or_found, clean=clean).\
        assign(software="mm")

        df_pd = _get_proteins_found_count(source="pd", quant_or_found=quant_or_found, clean=clean).\
        assign(software="pd")

        df = pd.concat(
            objs=[df_mm, df_pd],
            axis="index",
            join="outer",
            ignore_index=True,
        )

    else:

        if clean:
            index_cols = ["sample"]
        else:
            index_cols = ["sample", "bad", "np"]

        df = load_protein_table(source=source, quant_or_found=quant_or_found, clean=clean).\
        set_index(index_cols).\
        notna().\
        sum(axis=1).\
        sort_values().\
        to_frame().\
        reset_index().\
        rename(columns={0: "quant_proteins_count"})

        sample_split = df["sample"].str.split("_", expand=True, regex=False)
        df = df.assign(
            sample_type=sample_split[0],
            sample_condition=sample_split[1].replace(to_replace=r"^0[12]$", value=np.nan, regex=True) # Replace ID numbers from blanks with NaN
        )

        if not clean:
            df = df.assign(status=np.where(df["bad"] & ~df["np"], "Contaminated", "Normal"))
            df = df.assign(status=np.where(df["np"] & ~df["bad"], "No protein", df["status"]))
            df = df.assign(status=np.where(df["np"] & df["bad"], "Contaminated, No protein", df["status"]))
            df = df.drop(columns=["bad", "np"])

    return df
