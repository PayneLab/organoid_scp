import os
import pandas as pd

def load_table(
    source: str,
    name: str,
) -> pd.DataFrame:
    """
    Load one of the data tables.

    Parameters
    ----------
    source : {"pd", "mm"}
        Which tool your desired table was generated by. "pd" for Proteome 
        Discoverer, "mm" for MetaMorpheus.
    name : {
        "AllQuantifiedProteinGroups",
        "AnnotationProteinGroups",
        "PathwayProteinGroups",
        "ProteinGroups",
        "Proteins",
        "StudyInformation",
    }
        The name of the table you want. "AllQuantifiedProteinGroups" is the
        only option for MetaMorpheus; the rest are from Proteome
        Discoverer.
    """
    path_here = os.path.abspath(os.path.dirname(__file__))
    data_path = os.path.join(path_here, "data")

    if source == "mm":
        data_path = os.path.join(data_path, "metamorpheus", "Task1-SearchTask")
        if name == "AllQuantifiedProteinGroups":
            file_path = os.path.join(data_path, f"{name}.tsv")
        else:
            raise ValueError(f"Invalid argument '{name}' for name parameter with '{source}' as source.")
    elif source == "pd":
        data_path = os.path.join(data_path, "proteome_discoverer")
        if name in [
            "AnnotationProteinGroups",
            "PathwayProteinGroups",
            "ProteinGroups",
            "Proteins",
            "StudyInformation",
        ]:
            file_path = os.path.join(data_path, f"Payne_Organoids_{name}.txt")
        else:
            raise ValueError(f"Invalid argument '{name}' for name parameter with '{source}' as source.")
    else:
        raise ValueError(f"Invalid argument '{source}' for source parameter.")

    df = pd.read_csv(file_path, sep="\t")

    return df
