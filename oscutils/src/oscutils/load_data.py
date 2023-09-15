import numpy as np
import os
import pandas as pd

def load_protein_table(
    source: str,
    quant_or_found: str,
    clean: bool = True,
    load_old: bool = False,
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
        Whether to drop contaminated, no protein, QC, and blank samples from the table.
        This is only relevant if you're loading the old run, because those samples
        weren't included in the new runs.
    load_old: bool, default False
        Whether to load the data from the old run instead of the most recent one.
        This is only for backwards compatibility and verification as we switched to
        the data from the new runs; the old runs should not be used for any analysis.

    Returns
    -------
    pd.DataFrame
        The formatted protein abundance table.
    """

    # Validate arguments
    if source not in ["mm", "pd"]:
        raise ValueError(f"Invalid argument '{source}' for source parameter.")
    if quant_or_found not in ["quant", "found"]:
        raise ValueError(f"Invalid argument '{quant_or_found}' for quant_or_found parameter.")

    # Load and format tables
    if source == "mm":

        if quant_or_found == "found":
            raise ValueError("'found' data not available from MetaMorpheus.")

        df = _load_original_table("mm", "AllQuantifiedProteinGroups", load_old=load_old) # Load protein table
        ed = _load_original_table("mm", "ExperimentalDesign", load_old=load_old) # Load experimental design table.

        df = df.set_index("Protein Accession") # Index by protein accession
        df = df[df.columns[df.columns.str.startswith("Intensity_")]] # Select just the protein abundance columns

        df.columns = df.columns.to_series().str.split("_", n=1, expand=True)[1].values # Cut off "Intensity_" from column names

        # Standardize axis names
        df.columns.name = "sample"
        df.index.name = "protein"

        df = df.\
        transpose().\
        reset_index(drop=False)

        # Cut the ".raw" off the end of the filenames in the experimental design table, and generate a 
        # temporary "sample" column for joining with protein table. We'll then drop this "sample" column
        # and create another one based off of the old experimental design table for consistency between
        # loading versions. But we need this one to join in the case of the new tables.
        ed = ed.assign(
            FileName=ed["FileName"].str.rsplit(".", n=1, expand=True)[0],
            sample=ed["Condition"] + "_" + ed["Biorep"].astype(str),
        )[["FileName", "sample"]]

        df = ed.merge(
            right=df,
            on="sample",
            how="outer",
        ).drop(columns="sample") # Don't need it anymore

        # We're going to use the sample identifiers from the old experimental design table for consistency
        old_ed = _load_original_table("mm", "ExperimentalDesign", load_old=True)

        # Cut the ".raw" off the end of the filenames, and generate a "sample" column for identifying samples
        # across the whole project
        old_ed = old_ed.assign(
            FileName=old_ed["FileName"].str.rsplit(".", n=1, expand=True)[0],
            sample=old_ed["Condition"] + "_" + old_ed["Biorep"].map("{:0>2}".format),
        )

        # Mark the contaminated and no protein samples
        old_ed = old_ed.assign(
            contaminated=old_ed["FileName"].str.endswith("bad"),
            no_protein=old_ed["FileName"].str.endswith("np"),
        )

        # Select only the needed columns
        old_ed = old_ed[["sample", "FileName", "contaminated", "no_protein"]]

        # Join old experimental design into the protein table
        df = old_ed.merge(
            right=df,
            on="FileName",
            how="right",
        )

    elif source == "pd":

        df = _load_original_table("pd", "Proteins", load_old=load_old) # Load protein table
        inpf = _load_original_table("pd", "InputFiles", load_old=load_old) # Load input files table
        mm_ed = _load_original_table("mm", "ExperimentalDesign", load_old=True) # Load MetaMorpheus experimental design table so we can get same sample numbers. Note that we always load the old one because the new one doesn't have all of the sample numbers.

        df = df.set_index("Accession") # Index by protein accession

        # Select only rows with q <= 0.01 and Contaminant=False
        if load_old:
            qval_col = "Exp q-value Combined"
        else:
            qval_col = "Exp. q-value: Combined"

        df = df[~df["Contaminant"] & (df[qval_col] <= 0.01)]

        # Get either the protein quantification or the protein found columns
        if load_old:
            quant_startswith = "Abundances Scaled F"
            quant_extract = r"Abundances Scaled (F\d+) Sample"
            found_startswith = "Found in Sample F"
            found_extract = r"Found in Sample (F\d+) Sample"
        else:
            quant_startswith = "Abundance: F"
            quant_extract = r"Abundance: (F\d+): Sample"
            found_startswith = "Found in Sample: F"
            found_extract = r"Found in Sample: (F\d+): Sample"

        if quant_or_found == "quant":
            df = df[df.columns[df.columns.str.startswith(quant_startswith)]]
            df.columns = df.columns.to_series().str.extract(quant_extract)[0].values
        elif quant_or_found == "found":
            df = df[df.columns[df.columns.str.startswith(found_startswith)]]
            df.columns = df.columns.to_series().str.extract(found_extract)[0].values

        # Get axes names ready for joining
        df.columns.name = "File ID"
        df.index.name = "protein"

        df = df.\
        transpose().\
        reset_index(drop=False)

        # Get just the file names from the full file paths in the input files table
        inpf = inpf.assign(FileName=inpf["File Name"].str.rsplit(".", n=1, expand=True)[0].str.rsplit("\\", n=1, expand=True)[1])
        inpf = inpf[["File ID", "FileName"]]

        # Cut the ".raw" off the end of the filenames, and generate a "sample" column for joining with protein table
        mm_ed = mm_ed.assign(
            FileName=mm_ed["FileName"].str.rsplit(".", n=1, expand=True)[0],
            sample=mm_ed["Condition"] + "_" + mm_ed["Biorep"].map("{:0>2}".format),
        )

        # Mark the contaminated and no protein samples
        mm_ed = mm_ed.assign(
            contaminated=mm_ed["FileName"].str.endswith("bad"),
            no_protein=mm_ed["FileName"].str.endswith("np"),
        )

        # Select only the needed columns
        mm_ed = mm_ed[["sample", "FileName", "contaminated", "no_protein"]]

        # Merge the input files table with the experimental design table from MetaMorpheus so we can get the same sample numbers
        # We join based off of the original file name, which had the well plate number and coordinates of the sample
        inpf = inpf.merge(
            right=mm_ed,
            on="FileName",
            how="left",
        )

        # Get just the columns we want after the join
        inpf = inpf[["File ID", "FileName", "sample", "contaminated", "no_protein"]]

        # Join with the protein table
        df = inpf.\
        merge(
            right=df,
            on="File ID",
            how="outer",
        ).\
        drop(columns="File ID")


    # Split out the sample type, condition, and number from the metadata columns into separate columns
    filename_split = df["FileName"].str.replace("_2_", "_").str.split("_", expand=True, regex=False)
    sample_rsplit = df["sample"].str.rsplit("_", n=1, expand=True)

    # Since the blank and QC samples aren't "healthy" or "unhealthy", we're going to mark the
    # sample condition as NaN for them. This variable is used to do so in the assign statement below.
    conditions_to_replace = filename_split[1][~filename_split[1].isin(["Healthy", "Unhealthy"])].unique()

    # Create metadata columns
    df = df.assign(
        sample_type=filename_split[0].replace(to_replace="Psuedo-bulk", value="pbulk").str.lower(),
        sample_condition=filename_split[1].replace(to_replace=conditions_to_replace, value=np.nan).str.lower(),
        sample_num=sample_rsplit[1].astype(int),
    )

    # Move all the metadata columns to the beginning of the table
    df = df.\
    set_index(["sample", "sample_type", "sample_condition", "sample_num", "contaminated", "no_protein"]).\
    reset_index(drop=False)

    # Drop contaminated, no protein, blank, and QC samples if requested
    if clean:
        df = df[~df["contaminated"] & ~df["no_protein"] & ~df["sample"].str.startswith("qc_") & ~df["sample"].str.startswith("blank_")]
        df = df.\
        drop(columns=["contaminated", "no_protein"])

    # Sort to look nice
    df = df.\
    sort_values(by="sample").\
    reset_index(drop=True)

    return df

def _load_original_table(
    source: str,
    name: str,
    load_old: bool,
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
        "InputFiles",
    }
        The name of the table you want. "AllQuantifiedProteinGroups" is the
        only option for MetaMorpheus; the rest are from Proteome
        Discoverer.
    load_old: bool, default False
        Whether to load the data from the old run instead of the most recent one.
        This is only for backwards compatibility and verification as we switched to
        the data from the new runs; the old runs should not be used for any analysis.

    Returns
    -------
    pd.DataFrame
        The requested table.
    """
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(path_here, "data")

    # Validate and process arguments
    if source == "mm":

        if load_old:
            run_dir = "run20230411"
            usecols_rng = range(136)
        else:
            run_dir = "run20230508"
            usecols_rng = range(105)
        
        data_path = os.path.join(data_path, "metamorpheus", run_dir, "Task1-SearchTask")
        if name in [
            "AllQuantifiedProteinGroups",
            "ExperimentalDesign",
        ]:
            file_path = os.path.join(data_path, f"{name}.tsv")
        else:
            raise ValueError(f"Invalid argument '{name}' for name parameter with '{source}' as source.")
    elif source == "pd":

        if load_old:
            run_dir = "run20230410"
            filename = f"Payne_Organoids_{name}.txt"
        else:
            run_dir = "run20230809"
            filename = f"Payne_organoids2.0_{name}.txt"
        
        data_path = os.path.join(data_path, "proteome_discoverer", run_dir)
        if name in [
            "AnnotationProteinGroups", # Coverage numbers and percentages for GO annotations
            "PathwayProteinGroups", # Coverage numbers and percentages for Reactome and WikiPathways pathways
            "ProteinGroups", # Coverage for protein groups
            "Proteins",
            "InputFiles",
        ]:
            file_path = os.path.join(data_path, filename)
        else:
            raise ValueError(f"Invalid argument '{name}' for name parameter with '{source}' as source.")
    else:
        raise ValueError(f"Invalid argument '{source}' for source parameter.")

    if source == "mm" and name == "AllQuantifiedProteinGroups":
        # Whenever MetaMorpheus generates this file, all the lines except the
        # first one have an extra tab character at the end for some reason.
        # This causes the headers to be read in offset by one column. This line
        # addresses that.
        df = pd.read_csv(file_path, sep="\t", usecols=usecols_rng, low_memory=False)
    else:
        df = pd.read_csv(file_path, sep="\t", low_memory=False)

    return df

def _get_proteins_found_count(
    source: str,
    quant_or_found: str,
    clean: bool,
) -> pd.DataFrame:
    """
    Get the number of proteins found for each sample.

    Parameters
    ----------
    source : {"both", "pd", "mm"}
        Which tool you want to use protein abundance data from. "pd" for
        Proteome Discoverer, "mm" for MetaMorpheus, "both" for both.
    quant_or_found : {"quant", "found"}
        Whether to use "proteins quantified" or just "proteins found" data.
        "found" is only available from Proteome Discoverer, so do not use
        "found" with the "mm" or "both" source options.
    clean : bool
        Whether to drop contaminated and no protein samples from the table.
        Setting this parameter to False is only relevant if you're loading
        the old runs, because those samples weren't included in the new runs.
        However, the old runs should not be used for any analysis, so I didn't
        even write an ability in this function to use them. So just always
        leave this parameter as True, False won't do anything useful now.

    Returns
    -------
    pd.DataFrame
        The counts of proteins found for each sample.
    """

    # If data from both MetaMorpheus and Proteome Discoverer are requested, recursively call this
    # function to get the individual tables, add a column to specify which tool each count came
    # from, and concatenate the tables
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

        index_cols = ["sample", "sample_type", "sample_condition", "sample_num"]
        if not clean:
            index_cols += ["contaminated", "no_protein"]

        # Get the count of proteins found for each sample
        df = load_protein_table(source=source, quant_or_found=quant_or_found, clean=clean).\
        set_index(index_cols).\
        notna().\
        sum(axis=1).\
        sort_values().\
        to_frame().\
        reset_index().\
        rename(columns={0: "quant_proteins_count"})

        # Replace the "contaminated" and "no_protein" tables with a single column called "status" that has
        # a value of either "Normal", "Contaminated", "No Protein", or "Contaminated, No protein" for each
        # sample
        if not clean:
            df = df.assign(status=np.where(df["contaminated"] & ~df["no_protein"], "Contaminated", "Normal"))
            df = df.assign(status=np.where(df["no_protein"] & ~df["contaminated"], "No protein", df["status"]))
            df = df.assign(status=np.where(df["no_protein"] & df["contaminated"], "Contaminated, No protein", df["status"]))
            df = df.drop(columns=["contaminated", "no_protein"])

    return df
