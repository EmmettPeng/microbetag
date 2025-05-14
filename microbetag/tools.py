import os
import sys
import glob
import shutil

import subprocess
import multiprocessing
from typing import List, TYPE_CHECKING

from .utils import (
    get_files_with_suffixes,
    get_library_version,
    get_tool_location,
    ensure_flashweave_format,
    ensure_same_namespace_after_fw,
    mtg_logger,
)

from .seed_complementarity import ExportSeedComplementarities


if TYPE_CHECKING:
    from .config import Config


_logger_ = mtg_logger(__name__)


def run_seed_complementarity(config: "Config"):
    """
    Invoke PhyloMInt as edited from microbetag team to support parallel calculation of the seed and non seed sets
    and save corresponding sets to json files.
    """

    if config.prev_conf is None or os.path.exists(config.prev_conf) is False:
        config.skip_sets = False
        if config.users_models:
            genre_files = [
                os.path.join(config.for_reconstructions, file)
                for file in os.listdir(config.for_reconstructions)
            ]
            for file in genre_files:
                dest_path = os.path.join(config.genres, os.path.basename(file))
                try:
                    shutil.copy(file, dest_path)
                except Exception:
                    pass
        else:
            all_files = [
                os.path.join(config.reconstructions, file)
                for file in os.listdir(config.reconstructions)
            ]
            genre_files = [file for file in all_files if file.startswith(".xml")]
            for file in genre_files:
                dest_path = os.path.join(config.genres, os.path.basename(file))
                shutil.move(file, dest_path)

    seeds = ExportSeedComplementarities(config)

    # If no seed and/or non-seed sets are missing, get them
    if seeds.skip_sets is False and seeds.api is False:
        seeds.get_sets()

    # If either the phylomint scores file or the one with the seed complementarities (pckl) is missing, exract them
    if seeds.get_scores or seeds.get_complements:
        seeds.get_scores_and_compls()


def hmmsearch(params: List):
    """
    Function to invoke hmmsearch software.

    Arguments:
        params: list of parameters to be passed to the hmmsearch function.

    Note:
        We use only 1 cpu since we use a muliprocessing.Pool in the kegg_annotation(). 
    """
    (threshold_method, threshold, outtype, output, hmm_db, faa) = params

    cmd_para = [
        "hmmsearch",
        threshold_method,
        threshold,
        "--cpu",
        "1",
        "-o /dev/null",
        outtype,
        output,
        hmm_db,
        faa,
    ]

    cmd = " ".join(cmd_para)

    try:
        os.system(cmd)
    except Exception:
        _logger_.warning("Something wrong with KEGG hmmsearch!")


def run_prodigal(fasta, basename, outdir):
    """
    Function to predict ORFs using Prodigal.
    By default outdir is the ORFs folder
    fna	FASTA nucleic acid	Used generically to specify nucleic acids
    ffn	FASTA nucleotide of gene regions	Contains coding regions for a genome
    """

    faa_file = os.path.join(outdir, basename + ".faa")
    ffn_file = os.path.join(outdir, basename + ".ffn")
    fna_file = os.path.join(outdir, basename + ".fna")
    gbk_file = os.path.join(outdir, basename + ".gbk")

    PRODIGAL = get_tool_location("prodigal")
    cmd_para = [
        PRODIGAL,
        "-q",
        "-i",
        fasta,
        "-p",
        "meta",
        "-a",
        faa_file,
        "-d",
        ffn_file,
        "-o",
        gbk_file,
    ]
    cmd = " ".join(cmd_para)
    if os.path.exists(faa_file) or os.path.exists(fna_file) or os.path.exists(ffn_file):
        _logger_.info("ORFs already predicted for bin: %s", basename)
    else:
        _logger_.info("ORFs to be predicted for bin: %s", basename)
        try:
            os.system(cmd)
        except Exception:
            _logger_.warning("Something wrong with prodigal annotation!")


def kegg_annotation(faa: str, basename: str, out_dir: str, db_dir: str, ko_dic: dict, threads: int) -> bool:
    """
    Function to perform KEGG annotation.
    The function invokes hmmsearch.

    Args:
        faa:      Filepath to the .faa file of the bin in process
        basename: Bin id
        out_dir:  Path to output directory where .hmmout files will be stored
        db_dir:   Path to KEGG database directory
        ko_dic:
        threads:  Number of threads to be used
    """
    _logger_.info("KEGG annotation for %s", basename)
    params = []

    if not os.path.exists(faa):
        _logger_.warn(
            f"A .faa file for {basename} is not availavle. microbetag will skip it."
        )
        return False

    for knum, info in ko_dic.items():

        # Check if hmmout for particular KO of a certain bin is already there
        hmmout_filename = ".".join([knum, str(basename), "hmmout"])
        output = os.path.join(out_dir, basename, hmmout_filename)
        if os.path.exists(output):
            continue

        # Get hmm for KO under consideration
        hmm_db = os.path.join(db_dir, "profiles", knum + ".hmm")
        if not os.path.exists(hmm_db):
            continue

        # Set params for hmmsearch based on the ko_list annotation
        if info[1] == "full":
            threshold_method = "-T"
            outtype = "--tblout"

        elif info[1] == "domain":
            threshold_method = "--domT"
            outtype = "--domtblout"

        elif info[1] == "custom":
            threshold_method = "-E"
            outtype = "--tblout"

        params.append((threshold_method, info[0], outtype, output, hmm_db, faa))

    _logger_.info("Number of KEGG processes to be performed: %s", str(len(params)))
    process = multiprocessing.Pool(threads)
    process.map(hmmsearch, params)

    return True


def phenotrex_genotype(config: "Config"):
    """
    Runs the `compute-genotype` program of phenotrex to get COGs present in the list of genomes under study.

    Cite:
    """
    if not os.path.exists(config.genotypes_file):

        _logger_.info("Get phenotrex predictions")
        suffixes = [".fa", ".fasta", ".gz"]
        bin_files = get_files_with_suffixes(config.bins_path, suffixes)
        bin_files_in_a_row = " ".join(bin_files)

        # Build genotypes
        compute_genotype_params = []
        if not config.cwd.startswith("/microbetag"):
            compute_genotype_params += ["conda run -n phendb"]

        compute_genotype_params += [
            "phenotrex",
            "compute-genotype",
            "--out",
            config.genotypes_file,
            "--threads",
            str(config.threads),
            bin_files_in_a_row,
        ]
        # Container - case
        if config.cwd.startswith("/microbetag"):

            sk_init_version = get_library_version("scikit-learn")
            pheno_init_version = get_library_version("phenotrex")

            sk_required_version = "1.3.2"
            if sk_init_version != sk_required_version:
                # os.system("python3 -m pip install scikit-learn==1.3.2")
                get_sklearn_v = "==".join(["scikit-learn", sk_required_version])
                subprocess.run([sys.executable, "-m", "pip", "install", get_sklearn_v])

            pheno_required_version = "0.6.0"
            if pheno_init_version != pheno_required_version:
                get_pheno_v = "==".join(["phenotrex[fasta]", pheno_required_version])
                subprocess.run([sys.executable, "-m", "pip", "install", get_pheno_v])

            # Run the command
            compute_genotype_command = " ".join(compute_genotype_params)
            if os.system(compute_genotype_command) != 0:
                _logger_.info("Try phenotrex genotype for the second time.")
                if os.system(compute_genotype_command) != 0:
                    _logger_.error("Phenotrex compute genotype command failed.")
                    sys.exit(0)
        # Local case
        else:
            compute_genotype_command = " ".join(compute_genotype_params)
            try:
                subprocess.run(compute_genotype_command, shell=True, check=True)
            except Exception as e:
                _logger_.error(e)
                sys.exit(1)


def phenotrex_predict(config: "Config"):
    """ 
    Runs the `predict` program of `phenotrex` to predict whether a genome does have a trait or not.

    Cite:
    """
    # Get predictions
    phen_models = [
        os.path.join(config.phen_classes, model)
        for model in os.listdir(config.phen_classes)
    ]

    # Parse through the phen models (classes)
    for model in phen_models:
        model_name = os.path.basename(model)
        model_predictions_output = "".join(
            [config.predictions_path, "/", model_name[:-4], ".prediction.tsv"]
        )

        # Skip if predictions already computed
        if os.path.exists(model_predictions_output):
            _logger_.info("Predictions already exist for model: %s", model_name)

        else:

            predict_traits_params = []
            # Local case
            if not config.cwd.startswith("/microbetag"):
                predict_traits_params = ["conda run -n phendb"]
            # All cases: local and container
            predict_traits_params += [
                "phenotrex",
                "predict",
                "--genotype",
                config.genotypes_file,
                "--classifier",
                model,
                "--min_proba",
                str(config.min_proba),
                "--verb >",
                model_predictions_output,
            ]

            predict_trait_command = " ".join(predict_traits_params)

            try:
                subprocess.run(predict_trait_command, shell=True, check=True)

            except subprocess.CalledProcessError as e:
                _logger_.warn(f"Command execution failed with return code {e}")

    # If no predictions file is present for any classes, probably something went off.
    if not any(glob.glob(os.path.join(config.predictions_path, "*.prediction.tsv"))):
        _logger_.error(
            "No prediction was able to be retrieved with phenotrex. Check your input files."
        )
        sys.exit(0)


def run_manta(config: "Config"):
    """
    Runs the manta package to perform network clustering.

    Cite:
    """
    # Build the manta command
    _logger_.info("Running manta clustering algorithm.")
    manta_output_file = "/".join([config.output_dir, "manta_annotated"])
    manta_params = [
        "manta",
        "-i",
        config.base_network_file,
        "-f",
        "cyjs",
        "-o",
        manta_output_file,
        "--layout",
    ]
    manta_command = " ".join(manta_params)

    # Run manta
    try:
        if os.system(manta_command) != 0:
            warning_msg = (
                "The manta clustering algorithm failed."
                "Most likely this is because clusters could not be grouped based on the provided network "
                "and the parameters setup of manta."
                "Yet, microbetag will continue to the following steps without considering for clusters."
            )
            _logger_.warning(warning_msg)
            # Changed cfg so the buld_cx_annotated_graph function will not fail.
            config.network_clustering = False
        else:
            _logger_.info("""manta ran fine.""")

    except Exception as e:
        _logger_.warning(e)
        config.network_clustering = False


def run_flashweave(config: "Config"):
    """ 
    Runs FlashWeave to infer co-occurrence network. 

    Cite:
    """
    # Run FlashWeave
    from julia.api import Julia

    _logger_.info("Fix FlashWeave arguments from config.")

    # Checking for format support.
    if not config.onthefly:
        ensure_flashweave_format(conf=config)

    pair_args = set()
    for arg, values in config.flashweave_args.items():
        if values["required"]:
            if isinstance(values["value"], bool):
                pair_args.add((arg, str(values["value"]).lower()))
            else:
                _logger_.error(
                    f'You need to provide values for "{arg}" argument of FlashWeave.'
                )
                sys.exit(0)
        else:
            if values["value"] is not None:
                if values["type"] == "Bool":
                    pair_args.add((arg, str(values["value"]).lower()))
                else:
                    pair_args.add((arg, values["value"]))

    pair_args.add(("transposed", "true"))

    learn_in = ",".join(f"{arg[0]}={arg[1]}" for arg in pair_args)

    _logger_.info("Init Julia through Python!")
    jl = Julia(compiled_modules=False)
    jl.using("FlashWeave")

    # Run FlashWeave based on presence/absence of a metadata file
    if config.metadata_file:
        _logger_.info("Running FlashWeaeve along with a metadata file.")
        _logger_.info(
            f'save_network("{config.network}", learn_network("{config.flashweave_abd_table}", "{config.metadata_file}", {learn_in}))'
        )
        jl.eval(
            f'save_network("{config.network}", learn_network("{config.flashweave_abd_table}", "{config.metadata_file}", {learn_in}))'
        )
    else:
        _logger_.info("Running FlashWeaeve.")
        _logger_.info(
            f'save_network("{config.network}", learn_network("{config.flashweave_abd_table}", {learn_in}))'
        )
        jl.eval(
            f'save_network("{config.network}", learn_network("{config.flashweave_abd_table}", {learn_in}))'
        )

    ensure_same_namespace_after_fw(config)


def run_faprotax(config: "Config"):
    """
    Runs FAPROTAX collapse_table.py script to annotate taxa based on their taxonomy using literature.

    Cite:
    """

    faprotax_params = [
        "python3",
        config.faprotax_script,
        "-i",
        config.abundance_table,
        "-o",
        config.faprotax_funct_table,
        "-g",
        config.faprotax_txt,
        "-c",
        '"' + "#" + '"',
        "-d",
        '"' + config.taxonomy_column_name + '"',
        "-v",
        "--force",
        "-s",
        config.faprotax_sub_tables,
    ]
    faprotax_command = " ".join(faprotax_params)

    try:
        process = subprocess.Popen(
            faprotax_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        _, _ = process.communicate()
    except Exception:
        raise TypeError(
            f"Something went wrong when running FAPROTAX with your abuandance table: {config.abundance_table}"
        )
