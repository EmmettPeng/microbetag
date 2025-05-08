import os
import ast
import numpy as np
import pandas as pd

from .utils import safe_literal_eval
from .networks import get_edgelist

if __name__ == "__main__":
    # from variables import MAPPINGS, GENERA_NCBI_IDS, FAMILIES_NCBI_IDS
    from seed_complementarity import mtg_logger
    from utils import flatten_list
else:
    from .seed_complementarity import mtg_logger
    # from .variables import MAPPINGS, GENERA_NCBI_IDS, FAMILIES_NCBI_IDS
    from .utils import flatten_list, remove_nan_from_list


logger = mtg_logger(__file__)
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def map_seq_to_ncbi_tax_level_and_id(
    abd_tab, tax_col, seqId, tax_scheme, tax_delim, mappings, get_children=None
):
    """

    Returns 
        splt_tax: A pd.DataFrame with taxonomy levels, NCBI accession ids and
                  GTDB representative genomes (if available). A file callsed seq_map.tsv
                  with its content is saved in the main output folder.
    Example:
        splt_tax:
        >>> df.iloc[4,:]
        Unnamed: 0                              5
        domain                           Bacteria
        phylum                         Firmicutes
        class                          Clostridia
        order                     Oscillospirales
        family                   Oscillospiraceae
        genus                      Flavonifractor
        species            Flavonifractor plautii
        seqId                             ASV0012
        microbetag_id                     ASV0012
        extendedSpecies    Flavonifractor plautii
        refSpecies         Flavonifractor plautii
        species_ncbi_id                  411475.0
        gtdb_gen_repr         ['GCA_000239295.1']
        genus_ncbi_id                    946234.0
        family_ncbi_id                   216572.0
        ncbi_tax_id                      411475.0
        ncbi_tax_level                   mspecies
        Name: 4, dtype: object
    """
    # To-do [2025-04-23, Haris Zafeiropoulos]:
    # Check if there is an easy way to go on with get_children

    splt_tax = split_and_validate_taxonomy(abd_tab, seqId, tax_col, tax_delim)

    splt_tax = normalize_species_column(splt_tax)

    splt_tax, gtdb_on = add_species_ncbi_ids(splt_tax, tax_scheme, mappings)

    splt_tax = add_higher_level_ncbi_id(
        splt_tax, "genus", mappings.GENERA_NCBI_IDS, underscore=True
    )
    splt_tax = add_higher_level_ncbi_id(
        splt_tax, "family", mappings.FAMILIES_NCBI_IDS, underscore=True
    )

    get_ncbi_tax_level(splt_tax)

    splt_tax.loc[:, "gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(
        lambda x: str(x) if isinstance(x, list) else x
    )
    splt_tax = splt_tax.drop_duplicates()

    # Convert string representations of lists to actual lists
    splt_tax["gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )
    repr_genomes_present = remove_nan_from_list(flatten_list(splt_tax["gtdb_gen_repr"]))

    # # Make sure you only have one
    # splt_tax["extendedSpecies"] = splt_tax["extendedSpecies"].apply(
    #     lambda x: x[0] if isinstance(x, list) else x
    # )

    splt_tax["gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith("[") else x
    )

    return splt_tax, repr_genomes_present  # , None


def split_and_validate_taxonomy(abd_tab, seqId, tax_col, tax_delim):
    taxonomies = abd_tab[[seqId, tax_col, "microbetag_id"]]
    splt_tax = taxonomies[tax_col].str.split(tax_delim, expand=True)
    #
    if splt_tax[0].str.contains("Root").all():
        splt_tax = splt_tax.drop(columns=[0])
    #
    if splt_tax.shape[1] != 7:
        raise ValueError("The taxonomy scheme provided is not a 7-level one.")
    #
    splt_tax.columns = RANKS

    splt_tax[seqId] = taxonomies[seqId].values
    splt_tax["microbetag_id"] = taxonomies["microbetag_id"].values
    return splt_tax


def normalize_species_column(splt_tax):
    #
    species = splt_tax["species"].apply(process_underscore_taxonomy)
    genus = splt_tax["genus"].apply(process_underscore_taxonomy)
    #
    if splt_tax["species"].str.contains("s__").all():
        if splt_tax["species"].str.contains(r"__[\s_]").all():
            return splt_tax.assign(extendedSpecies=species)
        if any(g in s for g, s in zip(genus, species)):
            return splt_tax.assign(extendedSpecies=species)
        return splt_tax.assign(extendedSpecies=genus + " " + species)
    return splt_tax.assign(extendedSpecies=splt_tax["species"])


def add_species_ncbi_ids(splt_tax, tax_scheme, mappings):
    """
    Map taxonomy provided to microbetag taxonomy scheme to map NCBI Taxonomy id and level.

    """

    if tax_scheme     == "GTDB":
        taxon_map_file = "gtdbSpecies2ncbiId2accession.tsv"

    elif tax_scheme   == "Silva":
        taxon_map_file = "gtdb_silvaSpecies2ncbi2accession.tsv"

    elif tax_scheme   == "microbetag_prep":
        taxon_map_file = "gc_accession_16s_gtdb_ncbid.tsv"

    else:
        taxon_map_file = "species2ncbiId.tsv"

    #
    taxon_map_path  = os.path.join(mappings.TAXONOMY, taxon_map_file)
    taxon_to_ncbiId = pd.read_csv(taxon_map_path, sep="\t")

    #
    if tax_scheme in {"GTDB", "Silva", "microbetag_prep"}:
        taxon_to_ncbiId.columns = ["refSpecies", "species_ncbi_id", "gtdb_gen_repr"]
        splt_tax = pd.merge(
            splt_tax,
            taxon_to_ncbiId,
            left_on="extendedSpecies",
            right_on="refSpecies",
            how="left",
        )
        splt_tax["gtdb_gen_repr"] = splt_tax["gtdb_gen_repr"].apply(
            lambda x: [x] if pd.notna(x) else np.nan
        )
        return splt_tax, True

    else:

        # Fallback logic for fuzzy mapping if needed
        taxon_to_ncbiId.columns = ["refSpecies", "species_ncbi_id"]
        splt_tax = pd.merge(
            splt_tax,
            taxon_to_ncbiId,
            left_on="extendedSpecies",
            right_on="refSpecies",
            how="left",
        )
        return splt_tax, False


def add_higher_level_ncbi_id(splt_tax, level, mapping_file, underscore=True):
    """
    Will add column with NCBI Taxonomy Id of a taxonomic level higher than species/strain.
    """
    df = pd.read_csv(mapping_file, sep="\t", names=[level, "ncbi_tax_id"])
    if underscore:
        splt_tax[f"extended{level}"] = splt_tax[level].apply(
            process_underscore_taxonomy
        )
    else:
        splt_tax[f"extended{level}"] = splt_tax[level]

    splt_tax = pd.merge(
        splt_tax,
        df,
        left_on=f"extended{level}",
        right_on=level,
        how="left",
        suffixes=("", "_y"),
    ).rename(columns={"ncbi_tax_id": f"{level.lower()}_ncbi_id"})

    return splt_tax.drop(columns=[f"extended{level}", f"{level}_y"], errors="ignore")


def process_underscore_taxonomy(entry):
    """
    Convert "s__" like taxa to ncbi like ones
    """
    parts = entry.split("__")
    if len(parts) == 1:
        return entry
    if len(parts[1]) > 0:
        return parts[1].replace("_", " ")
    else:
        return np.nan


def get_ncbi_tax_level(splt_tax):
    """
    Assign the lowest available taxonomic level for each sequence in a
    splt_tax DataFrame. The function will first check if the sequence has
    a species level assignment, then a genus level assignment, then a family
    level assignment. If none of the above are available, the ncbi_tax_id
    column will be left as NaN.

    Parameters
    ----------
    splt_tax : pd.DataFrame

    Returns
    -------
    pd.DataFrame
        The same DataFrame with the additional column ncbi_tax_id
        containing the lowest available taxonomic level assignment for
        each sequence. If no assignment could be made, the value will be
        NaN.
    """

    # Apply the taxonomic assignment function
    # Assign taxonomic ID using the assign_tax_id_for_node_level function
    splt_tax["ncbi_tax_id"] = splt_tax.apply(assign_tax_id_for_node_level, axis=1)

    # Convert ncbi_tax_id to numeric, handling errors by coercing them to NaN
    splt_tax["ncbi_tax_id"] = (
        pd.to_numeric(splt_tax["ncbi_tax_id"], errors="coerce")
        .astype("Int64")
        .astype(str)
    )

    splt_tax["ncbi_tax_level"] = splt_tax.apply(assign_tax_level, axis=1)

    return splt_tax


def assign_tax_id_for_node_level(row):
    # A generator expression that returns the first available non-NaN value
    return next(
        (
            row[level]
            for level in ["species_ncbi_id", "genus_ncbi_id", "family_ncbi_id"]
            if not pd.isna(row[level])
        ),
        np.nan,
    )


def assign_tax_level(row):
    # Check if 'gtdb_gen_repr' is not NaN or empty, return 'mspecies' if valid
    if pd.notna(row["gtdb_gen_repr"]):
        return "mspecies"  # Or any other relevant level you wish to assign

    # Proceed to check for the lowest available taxonomic level
    taxonomic_levels = ["species", "genus", "family"]
    taxonomic_columns = ["species_ncbi_id", "genus_ncbi_id", "family_ncbi_id"]

    # Iterate over the taxonomic levels and columns, returning the first non-NaN taxonomic level
    for level, column in zip(taxonomic_levels, taxonomic_columns):
        if not pd.isna(row[column]):
            return level

    return "higher than family"


def otf_seqid_ncbi_gtdb__map(config):
    """
    Builds a dataframe with sequence ids of nodes A and B found associated in the 
    co-occurrence network followed by their corresponding NCBI Taxonomy ids and the 
    representative GTDB genomes.

    Returns:
        pairs_of_interest: {("",""). ("","")}
        relative_genomes: {ncbi_id: [gc, gc, gc], ..}
        mspecies_map_df: 
    """
    edgelist_df = get_edgelist(config.network)
    seq_map_df  = config.otf_seq_tax_df.dropna(subset=["species_ncbi_id"]).copy()

    merged = edgelist_df.merge(seq_map_df[['microbetag_id', 'species_ncbi_id', 'gtdb_gen_repr']],
                               left_on='node_A', right_on='microbetag_id', how='left') \
        .rename(columns={'species_ncbi_id': 'species_ncbi_id_A',
                         'gtdb_gen_repr': 'gtdb_gen_repr_A'}) \
        .drop(columns='microbetag_id')

    # Merge again to get info for nodeB
    merged = merged.merge(seq_map_df[['microbetag_id', 'species_ncbi_id', 'gtdb_gen_repr']],
                          left_on='node_B', right_on='microbetag_id', how='left') \
        .rename(columns={'species_ncbi_id': 'species_ncbi_id_B',
                         'gtdb_gen_repr': 'gtdb_gen_repr_B'}) \
        .drop(columns='microbetag_id')

    # Clean and apply transformations
    filtered = merged.dropna(subset=['gtdb_gen_repr_A', 'gtdb_gen_repr_B']).copy()

    filtered[['gtdb_gen_repr_A', 'gtdb_gen_repr_B']] = filtered[
        ['gtdb_gen_repr_A', 'gtdb_gen_repr_B']].applymap(safe_literal_eval)

    filtered[["species_ncbi_id_A", "species_ncbi_id_B"]] = filtered[
        ["species_ncbi_id_A", "species_ncbi_id_B"]].astype(int)

    # Explode and drop duplicates
    exploded        = filtered.explode('gtdb_gen_repr_A').explode('gtdb_gen_repr_B').reset_index(drop=True)
    mspecies_map_df = exploded.drop_duplicates(subset=['gtdb_gen_repr_A', 'gtdb_gen_repr_B'])

    # Build pairs of interest and relative genomes
    pairs_of_interest = {
        (str(row['species_ncbi_id_A']), str(row['species_ncbi_id_B']))
        for _, row in mspecies_map_df.iterrows()
    }
    pairs_of_interest.update({(b, a) for a, b in pairs_of_interest})
    relative_genomes = {
        str(row['species_ncbi_id_A']): {str(row['gtdb_gen_repr_A'])}
        for _, row in mspecies_map_df.iterrows()
    }
    relative_genomes.update(
        {str(row['species_ncbi_id_B']): {str(row['gtdb_gen_repr_B'])}
         for _, row in mspecies_map_df.iterrows()}
    )

    outfile = os.path.join(config.output_dir, "edge_map.tsv")
    mspecies_map_df.to_csv(outfile)

    return pairs_of_interest, relative_genomes, mspecies_map_df
