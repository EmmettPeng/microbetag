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

__version__ = "v1.0.3"
__author__ = "Haris Zafeiropoulos <haris.zafeiropoulos@kuleuven.be>"

import os, sys
import yaml
import logging

from .utils import *
from .tools import *
from .config import Config
from .genres import GEMSReconstruction
from .helpers import manta_input_net
from .build_mtg_cx2 import build_pseudo_cx, build_ndex2_net
from .pathway_complementarity import export_pathway_complementarities
from .seed_complementarity import ExportSeedComplementarities


# Set up custom logging format
logging.basicConfig(
    format='%(levelname)s: %(message)s',  # Define the format without "root:"
    level=logging.INFO  # Set the logging level
)


def run_microbetag(config):

    # ----------------
    # Build network if not available
    # ----------------
    if config.precalc_only:
        logging.info(
            "microbetag is about to perform the precalculations for your list of bins/MAGs only."
            "No network will be built."
        )
    elif not os.path.exists(config.network) or os.path.getsize(config.network) == 0:
        logging.info(
            "[STEP] NETWORK INFERENCE WITH FLASHWEAVE"
            "Using the abundance table provided, microbetag is about to build a co-occurrence network.\n"
        )
        run_flashweave(config)

    # ----------------
    # FAPROTAX
    # ----------------
    if config.abundance_table is not None:
        logging.info("[STEP] LITERATURE ANNOTATION WITH FAPROTAX")
        run_faprotax(config)

    # ----------------
    # phen annotations
    # ----------------
    if config.bins_ids is not None:
        logging.info("[STEP] PREDICTING PHENOTYPIC TRAITS")
        phenotrex_genotype(config=config)
        phenotrex_predict(config=config)

    # ----------------
    # Prodigal - ORF prediction
    # ----------------
    if (config.pathway_complementarity or config.seed_complementarity) and config.ko_merged is None:
        if len(os.listdir(config.prodigal)) != len(config.bins_ids):
            logging.info("[STEP  ] PREDICTING ORFs WITH PRODIGAL THROUGH DiTing")

            if config.bin_filenames is None:
                logging.error(
                    "Bins files have not been provided and they are required for the precalculation steps of microbetag."
                    "Provide the path to the directory with your bins/MAGs under the bins_fasta parameter of the config.yml file."
                )

            for bin_fa in config.bin_filenames:

                bin_filename = os.path.basename(bin_fa)

                bin_id, extension = os.path.splitext(bin_filename)
                bin_id = bin_id.split("/")[-1]

                bin_fa = os.path.join(config.bins_path, bin_fa)

                logging.info(f"Running Prodigal for {bin_id}")

                print(bin_fa, bin_id, config.prodigal)          # TODO: check if bin_id is actually only the basename of the whole path until extension

                run_prodigal(bin_fa, bin_id, config.prodigal)

    # ----------------
    # Pathway complementarity
    # ----------------
    if config.pathway_complementarity:
        # ----------------
        # KEGG annotation - based on the DiTing implementation // required in case of pathway complementarities
        # ----------------

        if config.ko_merged is None:

            logging.info("[STEP ] KEGG ANNOTATION OF THE ORFs \n")

            ko_list = os.path.join(config.kegg_db_dir, 'ko_list')
            ko_dic = ko_list_parser(ko_list)

            hmmout_dir = config.kegg_pieces_dir
            config.ko_merged = os.path.join(config.kegg_annotations, 'ko_merged.txt')

            for bn in config.bin_filenames:
                bin_id, extension = os.path.splitext(bn)
                bin_kos_dir = os.path.join(hmmout_dir, bin_id)
                os.makedirs(bin_kos_dir, exist_ok=True)

                for afile in os.listdir(bin_kos_dir):
                    if afile.endswith(".hmmout.all"):
                        continue

                for bn in config.bin_filenames:
                    faa = os.path.join(config.prodigal, bin_id + '.faa')
                    # A folder with KO predictions (a single hmmout file for each KO) per bin
                    check = kegg_annotation(faa, bin_id, config.kegg_pieces_dir, config.kegg_db_dir, ko_dic, config.threads)
                    # Out of the 24K hmmout files, make a single one with the predictions as backup and one with the 3-columns
                    if check:
                        bin_kos_to_file(hmmout_dir=bin_kos_dir , bin_id=bin_id)

            # Make the 3-columns files with all bins and their KOs
            merge_ko(config.kegg_pieces_dir, config.ko_merged)

        else:
            logging.info("A 3-col KEGG annotation file already available.")

        # ----------------
        # Extract pathway complementarities
        # ----------------
        pivot_df = load_merged_ko_file(config.ko_merged)  # Load ko_merged.txt

        if not os.path.exists(config.alts_file) or not os.path.exists(config.compl_file):

            logging.info("[STEP ] GET PATHWAY COMPLEMENTS ")
            # bin_kos_per_module, alt_to_gapfill, complements =
            _, _ = export_pathway_complementarities(
                config,
                pivot_df
            )

    # ----------------
    # Build GENREs
    # ----------------
    if config.seed_complementarity:

        if not config.users_models:

            logging.info("[STEP] GENOME-SCALE METABOLIC NETWORK RECONSTRUCTIONS")

            # Init reconstruction class
            build_genres = GEMSReconstruction(config)

            # Annotate step
            if config.input_for_recon_type == "bins_fasta":

                if config.genre_reconstruction_with == "modelseedpy":
                    build_genres.rast_annotate_genomes()  # saves under config.reconstructions

                elif config.gene_predictor == "prodigal":
                    logging.info("DiTing .faa files will be used")  # go to the .faa case, i.e., the ORFs/

                elif config.gene_predictor == "fragGeneScan":
                    logging.info("Get annotations with FragGeneScan.")
                    build_genres.fgs_annotate_genomes()   # saves under config.reconstructions

            elif config.input_for_recon_type == "coding_regions":
                logging.info("CarveMe will be used with the users .ffn-like files.")

            else:
                logging.warning(f"The combination of gene_predictor: {config.gene_predictor} \
                    \nand genre_reconstruction_with: {config.genre_reconstruction_with}, are not supported")

            # Reconstruct step
            if config.genre_reconstruction_with == "modelseedpy":
                build_genres.modelseed_reconstructions()

            elif config.genre_reconstruction_with == "carveme":
                build_genres.carve_reconstructions()

            else:
                logging.info("User models to be used for the seed complementarity step.")

    # ----------------
    # Phylomint
    # ----------------
    if config.seed_complementarity:
        logging.info("[STEP] COMPUTING SEED SETS AND SCORES")
        if not os.path.exists(config.phylomint_scores):
            run_phylomint(config)
        else:
            logging.info("Seed scores already computed.")

    # ----------------
    # Export seed complementarities
    # ----------------
    if config.seed_complementarity:
        logging.info("[STEP] EXPORTING SEED COMPLEMENTS")
        seed_complements = ExportSeedComplementarities(config)
        """
        [NOTE]:consider running again "seed scores" (PhyloMint) using update seed sets
        in this case, we should also edit the ConfidenceScore dictionary
        by removing seeds that were removed in the update()
        """
        if config.genre_reconstruction_with == "carveme":
            logging.info("We will map the BIGG compounds to ModelSEED ones.\
                \nIn the future, we will map BiGG ids to KEGG so we do not have to go through ModelSEED in this scenario.")
            seed_complements.map_carveme_seeds()

        if not os.path.exists(seed_complements.module_seeds):
            seed_complements.module_related_seeds()
        else:
            logging.info("Seed and non seed sets with compounds related to KEGG modules already retrieved.")

        if not os.path.exists(seed_complements.seed_complements):
            seed_complements.export_seed_complements()
            logging.info("Seed complements were exported fine.")
        else:
            logging.info("Seed complements already exported.")

    # ----------------
    # Network clustering
    # ----------------
    if config.network_clustering and config.prev_manta_net is None:

        logging.info("""[STEP]: network clustering using manta and the abundance table""")
        # Build original input file in cyjs format
        manta_input_net(config)

        logging.info(
            "Base network has been built and saved."
            "manta is now clustering your network..."
        )
        # Run manta on the cyjs network
        run_manta(config)

        logging.info("Base network has been built and saved.")

    # ----------------
    # Annotate network in .cx format
    # ----------------
    if config.precalc_only is False:
        logging.info("[STEP] ANNOTATE NETWORK ")
        annotated_network = build_pseudo_cx(config)
        with open(config.microbetag_annotated_network_file, "w") as f:
                annotated_network2file = convert_to_json_serializable(annotated_network)
                json.dump(annotated_network2file, f)
                logging.info("A microbetag-annotated network in .cx format was built sucessfully.")

        # Build cx2 with ndex2 library
        if build_ndex2_net(config.microbetag_annotated_network_file):
            # os.remove(config.microbetag_annotated_network_file)
            logging.info("The pseudo .cx file was converted to CX2 through NDEx successfully.")

    config.export_to_log()
    logging.info("A parameters.log file with the parameters used in this run was built.")

    logging.info("microbetag completed.")


def print_help():
    help_message = """
    Usage: python microbetag.py <path_to_config_yml>

    Other options:
    h        Display this help message.
    v        Display version.
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
    logging.error(conf_message)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Set white background in an SVG file.")
    parser.add_argument("--config", "-c", required=True, help="Path to the configuration yaml file.")
    parser.add_argument("--version", "-v", required=True, help="Show microbetag stand-alone tool version")
    parser.add_argument("--help", "-h", required=True, help="Show help message.")

    args = parser.parse_args()

    if args.version:
        print_version(); sys.exit()

    if args.help:
        print_help(); sys.exit()

    if args.config:

        with open(args.config, 'r') as yaml_file:
            try:
                config = Config(yaml.safe_load(yaml_file), args.config)
            except yaml.YAMLError as exc:
                print_config_message() ; sys.exit(0)

        # Run microbetag pipeline
        run_microbetag(config=config)

