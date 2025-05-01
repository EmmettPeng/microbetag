# microbetag : a software suite to annotate microbial co-occurrence networks

# Copyright (c) 2025 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file


import os
import pickle

import datetime
import ndex2.cx2
import pyshorteners
import pandas as pd

from .networks import get_edgelist, read_cyjson

from .utils import (
    mtg_logger,
    load_phenotypic_traits,
    extend_faprotax,
    extend_complements,
)
from .seed_complementarity import (
    load_seed_complement_files,
    build_url_with_seed_complements,
)


logger = mtg_logger(__name__)
COMPLEMENTARITY_TYPE, COMPLEMENTING = "complementarity", "(completed by)"
COOCCURENCE_DEPLETION = "co-occurrence/co-exclusion"
COOCCURRENCE, COOCCURYING = "co-occurrence", "(cooccurs with)"
DEPLETION, DEPLETING = "co-exclusion", "(negatively correlated with)"


# -----------------------------
# NODES
# -----------------------------
def taxonomy_levels(node):
    """
    Given a node with a taxonomy attribute, split to 7 taxonomic levels and gives their values to the corresponding ones

    node (dict):

    """
    try:
        levels = node["v"]["microbetag::taxonomy"].split(";")
        if len(levels) == 7:
            ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
            levels = [lvl.strip() for lvl in levels]

            for rank, value in zip(ranks, levels):
                node["v"][f"taxonomy::{rank}"] = value

            # Find last non-empty level
            for i in reversed(range(7)):
                if levels[i] and levels[i].split("_")[-1].strip():
                    node["v"]["microbetag::ncbi-tax-level"] = ranks[i]
                    break
    except Exception as e:
        pass  # or log the exception for debugging


def init_nodes_and_edges(edgelist, seq_id_to_taxonomy):
    """
    Builds nodes and edges.

    Args:
        edgelist (pd.DataFrame): DataFrame containing 'node_A', 'node_B', and 'microbetag::weight'.
        seq_id_to_taxonomy (Dict): Mapping from sequence IDs to taxonomy.

    Returns:
        Tuple[List[Dict], List[str], List[Dict]]: Nodes, node names, and edges.
    """
    # Extract unique node names efficiently
    node_names = list(set(edgelist["node_A"]).union(edgelist["node_B"]))

    # Create nodes
    nodes = [
        {
            "id": i,
            "v": {
                "name": name,
                "microbetag::taxonomy": seq_id_to_taxonomy.get(name, "Unknown"),
            },
        }
        for i, name in enumerate(node_names)
    ]

    # Enrich nodes with taxonomy levels
    for node in nodes:
        taxonomy_levels(node)

    # Create edges
    edges = [
        {
            "id": i,
            "s": node_names.index(node_a),
            "t": node_names.index(node_b),
            "v": {
                "interaction type": COOCCURRENCE if weight > 0 else DEPLETION,
                "interaction": COOCCURENCE_DEPLETION,
                "shared name": f"{node_a} {COOCCURYING if weight > 0 else DEPLETING} {node_b}",
                "microbetag::weight": weight,
            },
        }
        for i, node_a, node_b, weight in edgelist.itertuples(index=True, name=None)
    ]

    return nodes, node_names, edges


def update_with_phen_traits(config, nodes, node_names):
    """
    Updates nodes with phenotypic traits.

    Args:
        config (dict): Configuration for loading phenotypic traits.
        nodes (List[Dict]): List of node dictionaries.
        node_names (List[str]): List of node names corresponding to node IDs.
    """
    bin_phen_traits, _ = load_phenotypic_traits(config)

    for bin_name, phen_attributes in bin_phen_traits.items():
        try:
            node_index = node_names.index(bin_name)
        except ValueError:
            print(f"Warning: {bin_name} not found in node names.")
            continue  # Skip if bin_name is not in node_names

        node = nodes[node_index]

        node["v"]["microbetag::ncbi-tax-level"] = "mspecies"

        for trait, values in phen_attributes.items():
            node["v"][f"phendb::{trait}"]      = values["presence"]
            node["v"][f"phendbScore::{trait}"] = values["confidence"]


def update_with_faprotax_traits(config, nodes, node_names):
    """
    Updates nodes with FAPROTAX traits.

    Args:
        config (dict): Configuration for loading FAPROTAX traits.
        nodes (List[Dict]): List of node dictionaries.
        node_names (List[str]): List of node names corresponding to node IDs.
    """
    bin_faprotax_traits, _ = extend_faprotax(config)

    for bin_name, faprotax_attributes in bin_faprotax_traits.items():
        try:
            node_index = node_names.index(bin_name)
        except ValueError:
            print(f"Warning: {bin_name} not found in node names.")
            continue  # Skip if bin_name is not in node_names

        node = nodes[node_index]

        for trait in faprotax_attributes:
            node["v"][f"faprotax::{trait}"] = True


def update_with_manta(config, nodes, node_names):
    """
    Updates nodes with MANTA network cluster, assignment, and position data.

    Args:
        config (dict): Configuration containing the path to MANTA network data.
        nodes (List[Dict]): List of node dictionaries.
        node_names (List[str]): List of node names corresponding to node IDs.

    Returns:
        List[Dict]: MANTA layout with node positions.
    """
    manta_net = read_cyjson(config.manta_net)

    # Update nodes with cluster and assignment data
    for bin_name, cluster in manta_net.nodes(data="cluster"):

        try:
            nodes[node_names.index(bin_name)]["manta::cluster"] = cluster
        except ValueError:
            print(f"Warning: {bin_name} not found in node names.")

    for bin_name, assignment in manta_net.nodes(data="assignment"):

        try:
            nodes[node_names.index(bin_name)]["manta::assignment"] = assignment
        except ValueError:
            print(f"Warning: {bin_name} not found in node names.")

    # Store MANTA layout (node positions)
    manta_layout = [
        {"node": bin_name, "x": position["x"], "y": position["y"]}
        for bin_name, position in manta_net.nodes(data="position")
    ]

    return manta_layout


# -----------------------------
# EDGES - PATHWAY COMPLEMENTARITY
# -----------------------------
def pathway_complement_edge(
    edge_id, beneficiary, donor, complement, node_names, update=False, cx_edges=None
):
    """
    Creates an edge representing pathway complementarity between a donor and a beneficiary.

    Conceptually, the donor is the source, since a compound would be secreted from it and drive to the beneficiary (target)

    Args:
        beneficiary (str): The recipient of the complement.
        donor (str): The source of the complement.
        complement (any): The complement value.
        node_names (List[str]): List of node names to determine node indices.
        interaction_type (str): Type of interaction (e.g., "complementarity").
        interacting (str): Descriptor of the interaction.

    Returns:
        Dict: Edge representation.
    """
    if update:
        edge = cx_edges[edge_id]
    else:
        interacting = "".join(["(", COMPLEMENTING, ")"])
        edge = {
            "id": edge_id,  # Unique ID
            "s": node_names.index(donor),  # Source (donor)
            "t": node_names.index(beneficiary),  # Target (beneficiary)
            "v": {
                "interaction type": COMPLEMENTARITY_TYPE,
                "shared name": f"{beneficiary} {COMPLEMENTING} {donor}",
            },
        }

    column = f"compl::{beneficiary}:{donor}"
    edge["v"].update({column: _hat_complement(complement)})
    return edge


def _hat_complement(complements):
    """ """
    hat_compl = []
    for compl in complements.values():
        hat_compl.append("^".join(compl))
    return hat_compl


def _get_edge_id(beneficiary, donor, node_names, edges, interaction_type):
    """
    Checks if an edge is already there of a specific source - target -interaction type.
    If yes, it returns the index of the edge in the edges list.

    """
    target, source = node_names.index(beneficiary), node_names.index(donor)
    edge = next(
        (
            e
            for e in edges
            if e["s"] == source
            and e["t"] == target
            and e["v"]["interaction type"] == interaction_type
        ),
        None,
    )
    return (edges.index(edge), True) if edge is not None else (len(edges) + 1, False)


def _update_or_append(lst, index, item, update):
    if update:
        lst[index] = item
    else:
        lst.append(item)


def pathway_complements(config, edgelist_df, node_names, cx_edges):

    complements_dict = extend_complements(
        complements_json=config.compl_file,
        descrps_path=config.module_descriptions,
        max_scratch_alt=config.max_scratch_alt,
        pathway_complements_dir=config.pathway_complements_dir,
        pathway_complement_percentage=config.pathway_complement_percentage,
    )

    nodes_in_compls_dict = set(complements_dict.keys())

    for record in edgelist_df.iterrows():

        _, link = record
        node_a, node_b, _ = link

        if node_a in nodes_in_compls_dict and node_b in nodes_in_compls_dict:

            for beneficiary, donor in [(node_a, node_b), (node_b, node_a)]:

                edge_id, update = _get_edge_id(
                    beneficiary, donor, node_names, cx_edges, COMPLEMENTARITY_TYPE
                )

                complement = complements_dict[beneficiary][donor]

                pe = pathway_complement_edge(
                    edge_id,
                    beneficiary,
                    donor,
                    complement,
                    node_names,
                    update,
                    cx_edges,
                )

                _update_or_append(cx_edges, edge_id, pe, update)


# -----------------------------
# EDGES - SEED COMPLEMENTARITY
# -----------------------------
def _verbose_seed_complement(complements, beneficiarys_nonseed, kmap, shortener):
    """
    Appends the seed complementarities between two taxa as attributes to their corresponding edge
    id_x:
    id_y:
    """
    maps_in = list(kmap[kmap["modelseed"].isin(complements)]["map"].unique())
    complements_map = kmap[kmap["modelseed"].isin(complements)]
    beneficiarys_nonseeds_map = kmap[kmap["modelseed"].isin(beneficiarys_nonseed)]

    complements_verbose = [
        [
            kmap.loc[kmap["map"] == kegg_map, "category"]
            .unique()
            .item(),  # KEGG metabolism category
            kmap.loc[kmap["map"] == kegg_map, "description"]
            .unique()
            .item(),  # category description
            ";".join(
                set(
                    complements_map.loc[complements_map["map"] == kegg_map, "modelseed"]
                )
            ),  # seed coomplements in modelseed
            ";".join(
                set(
                    complements_map.loc[
                        complements_map["map"] == kegg_map, "kegg_compound"
                    ]
                )
            ),  # seed complements in kegg
            build_url_with_seed_complements(  # URL
                list(
                    complements_map.loc[
                        complements_map["map"] == kegg_map, "kegg_compound"
                    ]
                ),
                list(
                    beneficiarys_nonseeds_map.loc[
                        beneficiarys_nonseeds_map["map"] == kegg_map, "kegg_compound"
                    ]
                ),
                kegg_map,
                shortener,
            ),
        ]
        for kegg_map in maps_in
    ]

    merged_compl = ["^".join(gcompl) for gcompl in complements_verbose]

    return merged_compl


def seed_complement_edge(
    edge_id,
    beneficiary,
    donor,
    complement,
    competition,
    cooperation,
    node_names,
    update=False,
    cx_edges=None,
):
    """
    Conceptually, the donor is the source, since a compound would be secreted from it and drive to the beneficiary (target)
    his holds for the scores as well - for node_A in scores, we consider its seeds. Thus, the A of the score should be the beneficiary, i.e. target
    """
    if update:
        edge = cx_edges[edge_id]
    else:
        interacting = "".join(["(", COMPLEMENTING, ")"])

        edge = {
            "id": edge_id,
            "s": node_names.index(donor),
            "t": node_names.index(beneficiary),
            "v": {
                "interaction type": COMPLEMENTARITY_TYPE,
                "shared name": f"{beneficiary} {interacting} {donor}",
            },
        }

    # Edge attributes
    column = f"seedCompl::{beneficiary}:{donor}"
    edge["v"].update(
        {
            column: complement,
            "seed::competition": competition,
            "seed::cooperation": cooperation,
        }
    )

    return edge


def seed_complements(config, edgelist_df, node_names, cx_edges):

    shortener = pyshorteners.Shortener() if config.tinyurl else None
    kmap = load_seed_complement_files(config.kegg_mappings)

    seed_scores = pd.read_csv(
        config.phylomint_scores, sep="\t", header=None, skiprows=1
    )
    seed_scores.columns = ["A", "B", "Competition", "Complementarity"]

    # NOTE (Haris Zafeiropoulos, 2025-03-26):
    # Remember we only focus on the KEGG MODULES related seeds, nonseeds and compls
    with open(config.module_nonseeds, "rb") as f:
        non_seed_sets = pickle.load(f)

    with open(config.seed_complements, "rb") as f:
        seed_complements = pickle.load(f)

    seed_complements_dict = seed_complements.to_dict(orient="index")
    nodes_in_seed_compls_dict = set(seed_complements_dict.keys())

    for record in edgelist_df.iterrows():

        _, link = record
        node_a, node_b, _ = link

        if node_a in nodes_in_seed_compls_dict and node_b in nodes_in_seed_compls_dict:

            for beneficiary, donor in [(node_a, node_b), (node_b, node_a)]:

                edge_id, update = _get_edge_id(
                    beneficiary, donor, node_names, cx_edges, COMPLEMENTARITY_TYPE
                )

                scores = seed_scores.query("A == @beneficiary and B == @donor")
                if not scores.empty:
                    competAB, cooperAB = scores.iloc[0][
                        ["Competition", "Complementarity"]
                    ]
                else:
                    competAB, cooperAB = None, None  # or suitable defaults

                seed_complementAB = seed_complements_dict[beneficiary][donor]
                beneficiarys_nonseed = non_seed_sets.loc[beneficiary].to_list()[0]

                seed_complementAB = _verbose_seed_complement(
                    seed_complementAB, beneficiarys_nonseed, kmap, shortener
                )

                se = seed_complement_edge(
                    edge_id,
                    beneficiary,
                    donor,
                    seed_complementAB,
                    competAB,
                    cooperAB,
                    node_names,
                    update,
                    cx_edges,
                )

                _update_or_append(cx_edges, edge_id, se, update)


# -----------------------------
# NETWORK
# -----------------------------
def build_cx2(nodes, edges):

    # Create an empty net cx
    net_cx = ndex2.cx2.CX2Network()

    # Add nodes on the net_cx
    for node in nodes:
        node_attributes = node["v"]
        # Add node
        net_cx.add_node(attributes=node_attributes)

    # Add edges on the net_cx
    for edge in edges:

        source = edge["s"]
        target = edge["t"]
        attributes = edge["v"].copy()

        if attributes["interaction type"] in ["depletion", "cooccurrence"]:
            attributes["microbetag::weight"] = float(edge["v"]["microbetag::weight"])

        filtered_attributes = {
            key: value
            for key, value in attributes.items()
            if not (isinstance(value, list) and len(value) == 0)
        }

        # create an edge connecting the nodes, id of edge is returned
        _ = net_cx.add_edge(
            source=source, target=target, attributes=filtered_attributes
        )

    return net_cx


def mtg_annotate_network(config):

    edgelist_df = get_edgelist(config)

    # e.g.  'bin_31': 'd__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Reyranellales;f__Reyranellaceae;g__Reyranella;s__',
    seq_id_to_taxonomy_dic = config.seq_to_taxon_df.set_index(
        config.seq_to_taxon_df.columns[0]
    )[config.seq_to_taxon_df.columns[1]].to_dict()

    nodes, node_names, edges = init_nodes_and_edges(edgelist_df, seq_id_to_taxonomy_dic)

    if len(os.listdir(config.predictions_path)) > 0:
        update_with_phen_traits(config, nodes, node_names)

    if len(os.listdir(config.faprotax_sub_tables)) > 0:
        update_with_faprotax_traits(config, nodes, node_names)

    if config.net_cluster:  # os.path.exists(config.manta_net):
        manta_layout = update_with_manta(config, nodes, node_names)

    if config.path_compl:  # os.path.exists(config.compl_file):
        pathway_complements(config, edgelist_df, node_names, edges)

    if config.seed_compl:  # os.path.exists(config.seed_complements):
        seed_complements(config, edgelist_df, node_names, edges)

    net_cx = build_cx2(nodes, edges)

    timepoint = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    graphml_file = os.path.join(config.output_dir, f"mtag_net_{timepoint}.cx2")

    net_cx.set_network_attributes({"name": "microbetag annotated network"})
    net_cx.write_as_raw_cx2(graphml_file)
