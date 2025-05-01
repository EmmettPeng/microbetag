# microbetag : a software suite to annotate microbial co-occurrence networks

# Copyright (c) 2025 Haris Zafeiropoulos

# Licensed under GNU LGPL.3, see LICENCE file

import os
import json
import pandas as pd

from .helpers import (
    Emojis,
    Faprotax,
    BinsHandler,
    MappingPaths,
    NetworkHandler,
    AbdTableHandler,
    PathwayComplementarity,
    SeedComplementarityHandler,
)
from .networks import get_edgelist
from .utils import (
    mtg_logger,
    detect_separator,
    resolve_file_path,
    convert_to_json_serializable
)


logger = mtg_logger(__name__)


def load_abundance(abd_file):
    """
    Load a tsv/csv format abundance table assuming the sequence id is procided in the first column
    and the taxonomy in the last one

    :return seq_id_to_taxonomy: A pd.DataFrame with the sequence id and their corresponding genome
    :return sequence_id_column_name: The name (``str``) of the column with the sequence identifier (e.g. ``seqId``)
    :return taxonomy_column_name: The name (``str``) of the column with the taxonomy
    """

    delimiter                  = detect_separator(abd_file)
    abd_tab_df                 = pd.read_csv(abd_file, sep=delimiter)
    sequence_id_column_name    = abd_tab_df.columns[0]
    taxonomy_column_name       = abd_tab_df.columns[-1]
    seq_id_to_taxonomy         = abd_tab_df[[sequence_id_column_name, taxonomy_column_name]]
    seq_id_to_taxonomy.columns = ["sequence_id", "taxonomy"]
    return seq_id_to_taxonomy, sequence_id_column_name, taxonomy_column_name, delimiter


def _required(conf, key1, key2):
    value = conf.get(key1, {}).get(key2)
    if value is None:
        raise ValueError(f"Missing required configuration: {key1}.{key2}")
    return value


class Config:
    """
    Parses a microbetag configuration file (yaml) to init a microbetag run.
    """

    def __init__(self, conf, config_file=None):

        # Pass loaded yaml object
        self.yaml = conf

        # Load emojis
        emojis = Emojis()

        # User's group and user id
        if config_file is not None:
            self.yaml_file = os.stat(config_file)
            self.user_id   = self.yaml_file.st_uid
            self.group_id  = self.yaml_file.st_gid

        # config_wd = os.path.dirname(os.path.realpath(__file__))
        # self.cwd  = os.path.dirname(config_wd)

        self.cwd = os.path.dirname(os.path.realpath(__file__))

        # Check if microbetag runs as a container
        self.mount = None
        if self.cwd == "/microbetag":
            self.mount    = "/data"
            self.base_dir = self.mount
        else:
            self.base_dir = os.path.dirname(config_file)

        self.onthefly = _required(conf, "onthefly", "value")

        # Threads to be used
        self.threads = conf["threads"]["value"] if conf["threads"]["value"] else 2

        # Output dir
        output_dir = conf.get("output_directory", {}).get("dir_path")
        if output_dir:
            self.output_dir = os.path.join(self.base_dir, output_dir)
        else:
            raise ValueError("Output directory needs to be specified.")

        # Steps
        self.faprotax    = conf.get("faprotax_annotation", {}).get("value")
        self.phen_traits = conf.get("phenotrex_traits", {}).get("value")
        self.path_compl  = conf.get("pathway_complementarity", {}).get("value")
        self.seed_compl  = conf.get("seed_complementarity", {}).get("value")
        self.net_cluster = conf.get("network_clustering", {}).get("value")

        # Bins/MAGs/genomes
        bins_fasta     = conf.get("bins_fasta", {}).get("dir_path")
        self.bins_path = resolve_file_path(self.base_dir, bins_fasta)

        # The abundance table is now optional, if no abundance table and no network provided, then it will only run pre-calculations
        abd_tbl_filename     = conf.get("abundance_table_file", {}).get("file_path")
        self.abundance_table = resolve_file_path(self.base_dir, abd_tbl_filename)

        # Edgelist of provided network
        edge_list    = conf.get("edge_list", {}).get("file_path")
        self.network = resolve_file_path(
            self.base_dir, edge_list
        )  # os.path.join(self.base_dir, edge_list) if edge_list else None

        # NOTE: Setting precalculations only as true allows to get a Config instance
        # without providing an abundance table or network, meaning you can use the Config instance
        # for partial/specific tasks of microbetag.
        precalc_only      = conf.get("precalculations_only").get("value")
        self.precalc_only = precalc_only if precalc_only in [0, 1] else False

        if (
            self.abundance_table is None
            and self.network is None
            and precalc_only is False
        ):
            raise ValueError(
                f"You need to provide at least one between an abundance table and a network's edgelist in 3-column format."
            )

        # IMPORTANT: Sequence to taxonomy map -- required if no abundance table is needed
        sequence_taxonomy_map = conf.get("sequence_id_taxonomy_map", {}).get(
            "file_path"
        )
        self.sequence_taxonomy_map = resolve_file_path(
            self.base_dir, sequence_taxonomy_map
        )

        if (
            self.abundance_table is None
            and self.sequence_taxonomy_map is None
            and precalc_only is False
        ):
            raise ValueError(
                f"Since an abundance table is not provided, you need to provide a 2-column file with the sequence id (e.g bin ids)"
                "and their corresponding taxonomy or taxon name."
            )

        # IMPORTANT: Sequence id to taxonomy map
        if self.network is None and self.abundance_table:
            (
                self.seq_to_taxon_df,
                self.sequence_id_column_name,
                self.taxonomy_column_name,
                self.delimiter,
            ) = load_abundance(self.abundance_table)

            self.seq_ids = self.seq_to_taxon_df["sequence_id"].unique().tolist()

        elif self.abundance_table is None and self.network:

            seq_to_taxon_df = pd.read_csv(
                self.sequence_taxonomy_map, sep=self.delimiter
            )
            seq_to_taxon_df.columns = ["sequence_id", "taxonomy"]
            self.seq_to_taxon_df = seq_to_taxon_df
            self.seq_ids = self.seq_to_taxon_df["sequence_id"].unique().tolist()

        elif self.abundance_table and self.network:

            # NOTE: Not all sequence ids in the seq_ids need to have a taxonomy in this case -- only those coming from the abundance table
            # Yet, in case that the network has taxa not present in the abundance table, apparently it will lead to errors.
            network_df = get_edgelist(self)
            net_seq_ids = (
                pd.concat([network_df.iloc[:, 0], network_df.iloc[:, 1]])
                .unique()
                .tolist()
            )

            (
                self.seq_to_taxon_df,
                self.sequence_id_column_name,
                self.taxonomy_column_name,
                self.delimiter,
            ) = load_abundance(self.abundance_table)

            abd_seq_ids = self.seq_to_taxon_df["sequence_id"].unique().tolist()
            self.seq_ids = net_seq_ids + abd_seq_ids

        else:
            logger.warning(
                "Neither an abundance table nor a network was procided.\n"
                "microbetag will only run some pre-calculations not requiring them.\n"
                f"This is only good to use if you are are quite familiar with microbetag and you know what you are doing. {emojis.WARNING_EMOJI}"
            )

        # Set bins --- NOTE: CHECK FOR CONFLICTS
        self.bins_ids = None
        if self.bins_path is not None:
            bn = BinsHandler(config=self)
            self.__dict__.update(vars(bn))

        # Get pathway complementarity related variables
        pcompl = conf.get("pathway_complementarity", {}).get("value")
        self.pathway_complementarity = pcompl if pcompl in [0, 1] else True
        pc = PathwayComplementarity(config=self)
        self.__dict__.update(vars(pc))

        # Load abundance table with taxonomy
        nsc = AbdTableHandler(config=self)
        self.__dict__.update(vars(nsc))

        # Check whether bin names are the same in both abundance and edgelist files
        nh = NetworkHandler(config=self)
        self.__dict__.update(vars(nh))

        # Build output dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.predictions_path = os.path.join(self.output_dir, "phen_predictions")
        os.makedirs(self.predictions_path, exist_ok=True)

        # Open Reading Frames
        orfs = conf.get("orfs", {}).get("path")
        if orfs is None:
            self.prodigal = os.path.join(self.output_dir, "ORFs")
            os.makedirs(self.prodigal, exist_ok=True)
        else:
            self.prodigal = os.path.join(self.base_dir, orfs)

        # ModelSEEDpy arguments
        self.gapfill_model = conf["gapfill_model"]["value"]
        self.gapfill_media = conf["gapfill_media"]["value"]

        # Flashweave arguments
        self.metadata_file   = conf.get("metadata_file").get("file_path")
        self.metadata        = "false" if self.metadata_file == "false" else "true"
        self.flashweave_args = conf["flashweave_args"]

        # Phenotrex
        self.phen_classes   = os.path.join(self.cwd, "mtg_maps_models/phenDB/classes/")
        self.genotypes_file = os.path.join(self.output_dir, "train.genotype")
        min_proba           = conf.get("min_proba", {}).get("value")
        self.min_proba      = min_proba if not None else 0.6

        # Mappings
        mappings = MappingPaths()
        self.__dict__.update(vars(mappings))

        # FAPROTAX
        if self.abundance_table is not None:
            faprotax = Faprotax(config=self)
            self.__dict__.update(vars(faprotax))

        # Manta
        # net_clust = conf.get("network_clustering").get("value")
        # self.net_cluster = net_clust if net_clust in [0, 1] else False
        if self.net_cluster:
            self.prev_manta_net = conf.get("prev_clustered_network").get("file_path")
            self.manta_net = (
                os.path.join(self.base_dir, self.prev_manta_net)
                if self.prev_manta_net is not None
                else os.path.join(self.output_dir, "manta_annotated.cyjs")
            )
            self.base_network_file = os.path.join(self.output_dir, "basenet.cyjs")

        # Seed complementarity
        sc = SeedComplementarityHandler(config=self)
        self.__dict__.update(vars(sc))

        # Intermediate annoteted network file name
        self.microbetag_annotated_network_file = os.path.join(
            self.output_dir, "pseudo_cx_annotated_net.cx"
        )
        self.tinyurl = (
            conf.get("tinyurl", {}).get("value")
            if conf.get("tinyurl", {}).get("value")
            else False
        )
        # ==========
        # Init torch -- machine learning library
        # ==========
        if not self.onthefly or not self.phen_traits:

            import torch
            from deepnog.utils import get_weights_path
            from deepnog.utils import set_device

            device = set_device("auto")
            try:
                weights_path = get_weights_path(
                    database="eggNOG5",
                    level=str(2),
                    architecture="deepencoding",
                )
                _ = torch.load(weights_path, map_location=device)
            except Exception:
                logger.warn(
                    "Could not load the deepnog weights. Please check the deepnog installation and setup."
                )
                pass

        logger.info("Configuration file loaded successfully.")

    def export_to_log(self, log_file="parameters.log"):
        args = convert_to_json_serializable(self.__dict__)
        with open(log_file, "w") as f:
            json.dump(args, f)
