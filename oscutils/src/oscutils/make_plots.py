import altair as alt
import numpy as np

from .load_data import load_protein_table

def make_proteins_found_plot(
    source: str,
    clean: bool,
) -> alt.Chart:
    
    if clean:
        index_cols = "sample"
    else:
        index_cols = ["sample", "bad", "np"]
    
    df = load_protein_table(source=source, clean=clean).\
    set_index(index_cols).\
    notna().\
    sum(axis=1).\
    sort_values().\
    to_frame().\
    reset_index().\
    rename(columns={0: "quant_proteins_count"})

    if not clean:
        df = df.assign(status=np.where(df["bad"] & ~df["np"], "Contaminated", "Normal"))
        df = df.assign(status=np.where(df["np"] & ~df["bad"], "No protein", df["status"]))
        df = df.assign(status=np.where(df["np"] & df["bad"], "Contaminated, No protein", df["status"]))
        color = alt.Color(
            "status",
            scale=alt.Scale(
                domain=["Normal", "Contaminated", "No protein", "Contaminated, No protein"],
                range=["#4c78a8", "#e45756", "#f58518", "#b279a2"],
            ),
        )

    chart = alt.Chart(df).mark_bar().encode(
        x="sample",
        y="quant_proteins_count",
    )
    
    if not clean:
        chart = chart.encode(color=color)
    
    return chart