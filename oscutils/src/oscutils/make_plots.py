import altair as alt

from .load_data import (
    load_protein_table,
    _get_proteins_found_count,
)

def make_proteins_counts_plot(
    quant_or_found: str,
    clean: bool,
) -> alt.Chart:

    df = _get_proteins_found_count(source="both", quant_or_found=quant_or_found, clean=clean)

    if clean:
        df = df.assign(
            sample_type_condition=(df["sample_type"] + "_" + df["sample_condition"].fillna("")).str.replace(r"_$", "", regex=True),
        )

    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X(
            "software",
            axis=alt.Axis(
                title=" ",
                titleY=-320,
            ),
        ),
        y=alt.Y(
            "quant_proteins_count",
            axis=alt.Axis(
                title=None,
            ),
        ),
    )    

    if clean:
        chart = chart.encode(
            color="sample_type_condition",
        )
    else:
        chart = chart.encode(
            color=alt.Color(
                "status",
                scale=alt.Scale(
                    domain=["Normal", "Contaminated", "No protein", "Contaminated, No protein"],
                    range=["#4c78a8", "#e45756", "#f58518", "#b279a2"],
                ),
            ),
        )

    chart = chart.facet(
        facet=alt.Facet(
            "sample",
            header=alt.Header(
                labelOrient="bottom",
            ),
        ),
        columns=15,
    ).configure_facet(
        spacing=0
    ).resolve_axis(
        x="independent",
        y="independent",
    )

    return chart
