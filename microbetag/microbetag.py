"""
Aim:
    Perform the microbetag () approach annotating user's bins instead of mapping taxonomies against the
    microbetagDB representative genomes.

Input:
    - A folder with .fasta files of the corresponding bins
    - The abundance table of the bins across the samples
    - (optional) a co-occurrence network in a three-columns format

Output:
    - An annotated network in .cx2 format
"""

__version__ = "1.0.3"
__author__ = "Haris Zafeiropoulos <haris.zafeiropoulos@kuleuven.be>"

import os
import sys
import yaml
import argparse

from .utils import (
    mtg_logger,
    merge_ko,
    load_merged_ko_file,
    bin_kos_to_file,
    ko_list_parser,

)
from .tools import (
    run_flashweave,
    run_faprotax,
    phenotrex_genotype,
    phenotrex_predict,
    run_prodigal,
    kegg_annotation,
    run_seed_complementarity,
    run_manta,
)
from .db import (
    get_phen_traits,
    # something for path compls
)
from .config import Config
from .genres import GEMSReconstruction
from .helpers import manta_input_net
from .build_mtg_cx2 import mtg_annotate_network
from .pathway_complementarity import export_pathway_complementarities


logger = mtg_logger(__name__)


def run_microbetag(config):

    # ----------------
    # Build network if not available
    # ----------------
    if config.precalc_only:
        logger.info(
            "microbetag is about to perform the precalculations for your list of bins/MAGs only."
            "No network will be built."
        )
    elif not os.path.exists(config.network) or os.path.getsize(config.network) == 0:
        logger.info(
            "[STEP] NETWORK INFERENCE WITH FLASHWEAVE. "
            "Using the abundance table provided, microbetag is about to build a co-occurrence network.\n"
        )
        run_flashweave(config)

    # ----------------
    # FAPROTAX
    # ----------------
    if config.abundance_table is not None and config.faprotax:
        logger.info("[STEP] LITERATURE ANNOTATION WITH FAPROTAX")
        run_faprotax(config)

    # ----------------
    # phen annotations
    # ----------------
    if config.phen_traits:
        if config.bins_ids is not None and not config.onthefly:
            logger.info("[STEP] PREDICTING PHENOTYPIC TRAITS")
            phenotrex_genotype(config=config)
            phenotrex_predict(config=config)

        elif config.onthefly:
            get_phen_traits(config=config)

    # ----------------
    # Prodigal - ORF prediction
    # ----------------
    if (
        config.pathway_complementarity or config.seed_complementarity
    ) and config.ko_merged is None and not config.onthefly:
        if len(os.listdir(config.prodigal)) != len(config.bins_ids):

            logger.info("[STEP  ] PREDICTING ORFs WITH PRODIGAL THROUGH DiTing")

            if config.bin_filenames is None:
                logger.error(
                    "Bins files have not been provided and they are required for the precalculation steps of microbetag."
                    "Provide the path to the directory with your bins/MAGs under the bins_fasta parameter of the config.yml file."
                )

            for bin_fa in config.bin_filenames:

                bin_filename = os.path.basename(bin_fa)

                bin_id, _ = os.path.splitext(bin_filename)
                bin_id = bin_id.split("/")[-1]

                bin_fa = os.path.join(config.bins_path, bin_fa)

                logger.info(f"Running Prodigal for {bin_id}")

                print(
                    bin_fa, bin_id, config.prodigal
                )  # TODO: check if bin_id is actually only the basename of the whole path until extension

                run_prodigal(bin_fa, bin_id, config.prodigal)

    # ----------------
    # Pathway complementarity
    # ----------------
    if config.pathway_complementarity:

        # ----------------
        # KEGG annotation - based on the DiTing implementation // required in case of pathway complementarities
        # ----------------

        if config.ko_merged is None and not config.onthefly:

            logger.info("[STEP ] KEGG ANNOTATION OF THE ORFs \n")

            ko_list = os.path.join(config.kegg_db_dir, "ko_list")
            ko_dic = ko_list_parser(ko_list)

            hmmout_dir = config.kegg_pieces_dir
            config.ko_merged = os.path.join(config.kegg_annotations, "ko_merged.txt")

            for bn in config.bin_filenames:
                bin_id, _ = os.path.splitext(bn)
                bin_kos_dir = os.path.join(hmmout_dir, bin_id)
                os.makedirs(bin_kos_dir, exist_ok=True)

                for afile in os.listdir(bin_kos_dir):
                    if afile.endswith(".hmmout.all"):
                        continue

                for bn in config.bin_filenames:
                    faa = os.path.join(config.prodigal, bin_id + ".faa")
                    # A folder with KO predictions (a single hmmout file for each KO) per bin
                    check = kegg_annotation(
                        faa,
                        bin_id,
                        config.kegg_pieces_dir,
                        config.kegg_db_dir,
                        ko_dic,
                        config.threads,
                    )
                    # Out of the 24K hmmout files, make a single one with the predictions as backup and one with the 3-columns
                    if check:
                        bin_kos_to_file(hmmout_dir=bin_kos_dir, bin_id=bin_id)

            # Make the 3-columns files with all bins and their KOs
            merge_ko(config.kegg_pieces_dir, config.ko_merged)

        elif not config.onthefly:

            logger.info("A 3-col KEGG annotation file already available.")
            print("A 3-col KEGG annotation file already available.")

        # ----------------
        # Extract pathway complementarities
        # ----------------

        if config.onthefly:
            # TODO (Haris Zafeiropoulos, 2025-05-01):
            print("needs to be buiilts")

        else:

            pivot_df = load_merged_ko_file(config.ko_merged)  # Load ko_merged.txt

            if not os.path.exists(config.alts_file) or not os.path.exists(
                config.compl_file
            ):

                print("[STEP ] GET PATHWAY COMPLEMENTS ")
                print(pivot_df)
                # bin_kos_per_module, alt_to_gapfill, complements =
                _, _ = export_pathway_complementarities(config, pivot_df)

    # ----------------
    # Seed complementarity
    # ----------------
    if config.seed_complementarity:

        # In case of seed complementarity:
        # 1. we need to make sure we have GEMs, then
        # 2. we need to extract the seed and non-seed sets and their compls

        # ----------------
        # Build GENREs
        # ----------------
        if not config.users_models:

            logger.info("[STEP] GENOME-SCALE METABOLIC NETWORK RECONSTRUCTIONS")

            # Init reconstruction class
            build_genres = GEMSReconstruction(config)

            # Annotate step
            if config.input_for_recon_type == "bins_fasta":

                if config.genre_reconstruction_with == "modelseedpy":
                    build_genres.rast_annotate_genomes()  # saves under config.reconstructions

                elif config.gene_predictor == "prodigal":
                    logger.info(
                        "DiTing .faa files will be used"
                    )  # go to the .faa case, i.e., the ORFs/

                elif config.gene_predictor == "fragGeneScan":
                    logger.info("Get annotations with FragGeneScan.")
                    build_genres.fgs_annotate_genomes()  # saves under config.reconstructions

            elif config.input_for_recon_type == "coding_regions":
                logger.info("CarveMe will be used with the users .ffn-like files.")

            else:
                logger.warning(
                    f"The combination of gene_predictor: {config.gene_predictor} \
                    \nand genre_reconstruction_with: {config.genre_reconstruction_with}, are not supported"
                )

            # Reconstruct step
            if config.genre_reconstruction_with == "modelseedpy":
                logger.info("Build draft reconstructions with ModelSEEDpy")
                build_genres.modelseed_reconstructions()

            elif config.genre_reconstruction_with == "carveme":
                logger.info("Build draft reconstructions with carveme")
                build_genres.carve_reconstructions()

            else:
                logger.info("User models to be used for the seed complementarity step.")

        # ----------------
        # microbetag implementation of Phylomint
        # ----------------
        logger.info("[STEP] COMPUTING SEED SETS AND SCORES")
        run_seed_complementarity(config)

    # ----------------
    # Network clustering
    # ----------------
    if config.net_cluster and config.prev_manta_net is None:

        logger.info(
            """[STEP]: network clustering using manta and the abundance table"""
        )
        # Build original input file in cyjs format
        manta_input_net(config)

        logger.info(
            "Base network has been built and saved."
            "manta is now clustering your network..."
        )
        # Run manta on the cyjs network
        run_manta(config)

        logger.info("Base network has been built and saved.")

    # ----------------
    # Annotate network in .cx format
    # ----------------
    if config.precalc_only is False:
        logger.info("[STEP] ANNOTATE NETWORK ")
        mtg_net = mtg_annotate_network(config)

    # ----------------
    # Keep arguments
    # ----------------
    config.export_to_log()
    logger.info("A parameters.log file with the parameters used in this run was built.")

    logger.info("microbetag completed.")

    return mtg_net


def print_help():
    help_message = """
    Usage: microbetag --config <path_to_config_yml>

    Other options:
    -h        Display this help message.
    -v        Display version.
    """
    print(help_message)


def print_version():
    print(__version__)


def print_config_message():
    conf_message = """
    The config file you provided cannot be loaded.
    Please make sure you follow the instructions on the documentation site:
    https://hariszaf.github.io/microbetag/docs/tutorials/local/#input-and-configyml-files
    """
    logger.error(conf_message)


def main():

    parser = argparse.ArgumentParser(description="Microbetag CLI")

    parser.add_argument("--config", "-c", help="Path to the configuration yaml file.")
    parser.add_argument(
        "-v", "--version", action="store_true", help="Show Microbetag version"
    )

    args = parser.parse_args()

    if args.version:
        print_version()
        sys.exit()

    elif args.config is None:
        print_help()
        sys.exit(0)

    try:
        with open(args.config, "r") as yaml_file:
            config = Config(yaml.safe_load(yaml_file), args.config)

    except yaml.YAMLError:
        print_config_message()
        sys.exit(1)

    # Run microbetag pipeline
    run_microbetag(config=config)


if __name__ == "__main__":

    main()
