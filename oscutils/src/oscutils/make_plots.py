import altair as alt

from .load_data import (
    load_protein_table,
    _get_proteins_found_count,
)

def make_proteins_counts_plot(
    quant_or_found: str,
    clean: bool = True,
) -> alt.Chart:
    """
    Generate a plot showing the count of proteins found for each sample,
    including results for both the MetaMorpheus and the Proteome Discoverer
    data.

    Parameters
    ----------
    quant_or_found : {"quant", "found"}
        Whether to use "proteins quantified" or just "proteins found" data.
        Currently only works with "quant", as "proteins found" data isn't
        available from MetaMorpheus.
    clean : bool, default True
        Whether to drop contaminated and no protein samples from the table.
        If they are not dropped, bars are colored based on whether the
        sample was normal, contaminated, or no protein. If they are
        dropped, bars are colored based on the sample type (healthy or
        unhealthy HFL1 cell, healthy or unhealthy psuedo-bulk, healthy or
        unhealthy HFL1 boost, or blank).

        Setting this parameter to False is only relevant if you're loading
        the old runs, because those samples weren't included in the new runs.
        However, the old runs should not be used for any analysis, so I didn't
        even write an ability in this function to use them. So just always
        leave this parameter as True, False won't do anything useful now.

    Returns
    -------
    alt.Chart
        A faceted bar chart showing the protein counts from each software
        tool for each sample.
    """

    df = _get_proteins_found_count(source="both", quant_or_found=quant_or_found, clean=clean)

    if clean:

        # Join the "sample_type" and "sample_condition" columns into a single column by which we'll color code the bars
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
        # Color the bars based on the "sample_type_condition" column we made
        chart = chart.encode(
            color="sample_type_condition",
        )
    else:
        # Color code the bars based on whether the sample is normal, contaminated, or no protein
        chart = chart.encode(
            color=alt.Color(
                "status",
                scale=alt.Scale(
                    domain=["Normal", "Contaminated", "No protein", "Contaminated, No protein"],
                    range=["#4c78a8", "#e45756", "#f58518", "#b279a2"],
                ),
            ),
        )

    # Break the chart into individual facets for each sample
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
