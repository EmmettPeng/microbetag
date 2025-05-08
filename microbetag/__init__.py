from .config import Config

from .helpers import (
    MappingPaths,
    PathwayComplementarity,
    Faprotax,
    NetworkHandler,
    AbdTableHandler,
    SeedComplementarityHandler,
    BinsHandler,
    manta_input_net,
)  # GenresHandler

from .utils import (
    ko_list_parser,
    merge_ko,
    extend_complements,
    load_phenotypic_traits,
    extend_faprotax,
)

from .pathway_complementarity import (
    export_pathway_complementarities,
    all_complements,
    all_alternatives,
    build_kegg_url,
)

from .seed_complementarity import (
    ExportSeedComplementarities,
    load_seed_complement_files,
    build_url_with_seed_complements,
    kegg_module_related_intersect,
)
from .build_mtg_cx2 import mtg_annotate_network

from .db import (
    get_genomes_for_ncbi_tax_id,
    get_ncbi_tax_if_for_genome,
    get_patric_id_of_gc_accession_list,
    get_phen_traits,
    get_path_compls_for_ncbi_ids
)

# from .PhyloMint.PhyloMInt import PhylomintMGT

import os

_KEGG_MAPPINGS = os.path.join(
    os.path.dirname(__file__), "mtg_maps_models", "kegg_mappings"
)
KEGG_TERMS_PER_MODULE = os.path.join(_KEGG_MAPPINGS, "kegg_terms_per_module.tsv")
MODULE_DEFINITION_MAP = os.path.join(_KEGG_MAPPINGS, "module_definition_map.json")

# e.g cpd00020	C00022	M00001
KEGG_MODULES_TO_MAPS = os.path.join(_KEGG_MAPPINGS, "module_map_pairs.tsv")
