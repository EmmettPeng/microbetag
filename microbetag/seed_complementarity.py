import os
import sys
import json
import time

import pickle
import tarfile

import pandas as pd

from tqdm import tqdm
import multiprocessing

from .utils import mtg_logger
from .PhyloMint.lib import BuildGraphNetX
from .PhyloMint.lib.CalculateIndexes import calculate_scores, extract_complements


logger = mtg_logger(__name__)


class ExportSeedComplementarities:
    """
    Class to  export seed complements.
    Needs a config object to initiate it.

    # conda activate microbetag
    import yaml
    from config import Config

    config_file = "tests/dev_io_microbetag/config.yml"
    with open(config_file, 'r') as yaml_file:
        config = Config(yaml.safe_load(yaml_file), config_file)

    seed_complements = ExportSeedComplementarities(config)
    """

    def __init__(self, config):

        print("Initiating Export Seed Complementarities class..")

        self.dir_path       = config.genres
        self.outdir         = config.seeds
        self.scores_outfile = os.path.join(self.outdir, "phylomint_scores.tsv")
        self.save_dics      = True
        self.threads        = config.threads

        self.skip_sets     = config.skip_sets
        self.prev_conf     = config.prev_conf
        self.prev_nonseeds = config.prev_nonseeds

        self.perce_save          = 10
        self.seed_complements    = config.seed_complements
        self.module_seeds        = config.module_seeds
        self.module_nonseeds     = config.module_nonseeds

        self.get_scores          = not os.path.exists(self.scores_outfile)
        self.get_complements     = not os.path.exists(self.seed_complements)
        self.only_module_related = True

        self.namespace = "modelseed"
        if config.genre_reconstruction_with == "carveme":
            logger.info("Load bigg2seed map...")
            self.namespace = "BiGG"
            self.bigg2seed =  bigg_to_seed_mapping_df(config.metanetx_compounds)

        # Prefixes - suffixes
        self.ex_suffix = "e" if self.namespace == "BiGG" else "e0"
        self.int_suffix = "c" if self.namespace == "BiGG" else "c0"
        self.compound_prefix = "M"

        # self.seed_ko_mo = config.seed_ko_mo
        self.modelseed_compounds_of_interest = get_kegg_module_related(config.seed_ko_mo)

        if self.skip_sets:

            print("Load previously computed confidence scores and non-seed sets..")

            try:
                with open(self.prev_conf, "r") as f:
                    ConfidenceDic = json.load(f)
            except FileExistsError as e:
                raise e
            try:
                with open(self.prev_nonseeds, "r") as f:
                    nonSeedSetDic = json.load(f)
            except FileExistsError as e:
                raise e

            self.ConfidenceDic, self.nonSeedSetDic = {}, {}
            for k, v in ConfidenceDic.items():
                self.ConfidenceDic[k] = self._strip_pre_suff_from_dict(v)
            for k,v in nonSeedSetDic.items():
                self.nonSeedSetDic[k] = self._strip_pre_suff_from_list(v)


    def get_sets(self):
        """
        Get seed and non-seed sets for each model
        """
        # Build initial dictionaries
        SeedSetDic    = dict()
        nonSeedSetDic = dict()
        ConfidenceDic = dict()

        # Get all XML files in directory
        logger.info("Export seed and non seed sets....")

        sbml_files = [ os.path.join(self.dir_path, f) for f in os.listdir(self.dir_path) if f.endswith('.xml') ]
        num_threads = min(len(sbml_files), self.threads)
        with multiprocessing.Pool(processes=num_threads) as pool:
            with tqdm(total=len(sbml_files), desc="Processing SBML files to calculate seed and non-seed sets.") as pbar:
                results = []
                for result in pool.imap_unordered(self.process_sbml, sbml_files):
                    results.append(result)
                    pbar.update(1)  # Update progress bar as soon as a task completes

        # Unpack the results
        for result in results:
            try:
                sbml_base, SeedSet, nonSeedSet, SeedSetConfidence = result
            except:
                pass

            tmp                      = {key: None for key in SeedSet}
            SeedSetDic[sbml_base]    = tmp.keys()
            nonSeedSetDic[sbml_base] = nonSeedSet
            ConfidenceDic[sbml_base] = SeedSetConfidence

        pool.close()
        pool.join()

        self.SeedSetDic, self.nonSeedSetDic, self.ConfidenceDic = SeedSetDic, nonSeedSetDic, ConfidenceDic
        if self.save_dics:

            SeedSetDic_serial    = self._serialize_dic(SeedSetDic, os.path.join(self.outdir, 'SeedSetDic.json'))
            nonSeedSetDic_serial = self._serialize_dic(nonSeedSetDic, os.path.join(self.outdir, 'nonSeedSetDic.json'))
            _ = self._serialize_dic(ConfidenceDic, os.path.join(self.outdir, 'confidenceDic.json'))

            # NOTE (Haris Zafeiropoulos, 2025-03-26): In the pickle conversion we keep only the KEGG MODULE related
            self._dict_to_pickle(SeedSetDic_serial, self.module_seeds)
            self._dict_to_pickle(nonSeedSetDic_serial, self.module_nonseeds)

        logger.info("Seed and non seed sets have been exported.")


    def get_scores_and_compls(self):
        """
        Based on the seed and non-seed sets calculated, get all pairwise competition and cooperation scores and the seed complementarities
        between the modles under study.
        """

        logger.info("Exporting seed scores and complementarities.")

        total_species = self.ConfidenceDic.keys()
        self.perce_save = min(len(total_species), self.perce_save)  # Ensure perce_save does not exceed total species count

        # For multiprocessing
        lock    = multiprocessing.Lock()
        queue   = multiprocessing.Queue()
        manager = multiprocessing.Manager()
        shared_compls_dict = manager.dict()  # Shared dictionary for DataFrame
        processes = []

        # Start a separate process for tracking progress
        progress_process = multiprocessing.Process(target=progress_tracker, args=(queue, len(total_species)))
        progress_process.start()

        # Batch of processes..
        for i, species in enumerate(total_species):
            if i % self.threads == 0:
                for process in processes:
                    process.join()        # Wait for the batch to finish
                processes = []            # Clear completed processes

            # ..each running ther worker function
            process = multiprocessing.Process(target=self.worker_function, args=(lock, species, queue, shared_compls_dict))
            processes.append(process)
            process.start()

            # Periodically save progress in case of big data -- e.g. updating microbetagDB
            if len(total_species) > 500 and i % (len(total_species) // self.perce_save) == 0 and i != 0:

                print(f"Progress {i}/{len(total_species)} - Shared dict state: {len(shared_compls_dict)}", file=sys.stderr)
                for process in processes:
                    process.join()

                # Save temporary progress
                tmp_dict_serializable = {k: dict(v) for k, v in shared_compls_dict.items()}
                with open(f"tmp_cmpls_{i // (len(total_species) // self.perce_save)}.json", "w") as j:
                    json.dump(tmp_dict_serializable, j)
                shared_compls_dict = manager.dict()  # Reset shared dict after saving

        # Wait for any remaining processes to finish
        for process in processes:
            process.join()

        # Signal the progress tracker to stop
        queue.put(None)
        progress_process.join()

        # Finalize shared dictionary and convert to DataFrame
        final_compls_dict = {k: dict(v) for k, v in shared_compls_dict.items()}
        df = pd.DataFrame.from_dict(final_compls_dict)

        # Replace NaN with empty lists for species' complements with itself
        df = df.applymap(lambda x: [] if isinstance(x, float) and pd.isna(x) else x)

        # Save the final DataFrame
        with open(self.seed_complements, "wb") as f:
            pickle.dump(df, f)


    def scores_and_overlaps_for_a_species(self, species, lock, shared_dict):
        """
        [NEW] Get scores and complements for a specific model (species)
        """
        # Get species pairs and seedset confidence
        as_beneficiary, _ = generate_fixed_pairwise_comparisons(species, list(self.ConfidenceDic.keys()))
        species_seedset_confidence = self.ConfidenceDic[species]

        findings, compls = set(), {}
        species = species.replace(".PATRIC", "")  # Clean species name

        for partner in [pair[1] for pair in as_beneficiary if pair[1] != species]:
            partner_seedset_confidence, nonSeedB = self.ConfidenceDic[partner], self.nonSeedSetDic[partner]

            SeedA, SeedB, nonSeedB = set(species_seedset_confidence.keys()), set(partner_seedset_confidence.keys()), set(nonSeedB)

            if self.get_scores:
                MetabolicCooperationIdxAB, MetabolicCompetitionIdxAB = calculate_scores(SeedA, species_seedset_confidence, SeedB, nonSeedB)
                findings.add(f"{species}\t{partner.replace('.PATRIC', '')}\t{MetabolicCompetitionIdxAB}\t{MetabolicCooperationIdxAB}\n")

            if self.get_complements:
                B_complememts_to_A = extract_complements(SeedA, nonSeedB)
                if self.only_module_related:
                    B_complememts_to_A = kegg_module_related_intersect(B_complememts_to_A, self.modelseed_compounds_of_interest)
                compls[partner] = B_complememts_to_A

        # Update shared dictionary with complements if needed
        if self.get_complements:
            with lock:
                shared_dict[species] = compls

        # Write findings to file
        if self.get_scores:
            with open(self.scores_outfile, 'a') as f:
                f.writelines(findings)


    def worker_function(self, lock, species, queue, shared_dict):
        """Wrapper function to process a species and signal completion."""
        self.scores_and_overlaps_for_a_species(species, lock, shared_dict)
        with lock:
            queue.put(1)  # Signal that one task is completed


    def process_sbml(self, sbml_path, maxcc=2):
        """
        For each SBML model file (.xml) extract seeds, non-seeds and confidence scores
        using the PhyloMint adapted/refined approach of ours, i.e. building a directed graph
        with only the cytosol reactions, considering for the reversibility of a reaction.
        """
        filename = os.path.basename(sbml_path)
        sbml_base = filename.rstrip('.xml')

        # calculate SeedSets
        try:
            DG_sbml = BuildGraphNetX.buildDG(sbml_path)
        except:
            print("Failed to run PhyloMint for:", sbml_path)
            return

        # Get sets !
        # SeedSet: a dict_keys  |  nonSeedSet: a list already  |  SeedSetConfidence: a dict
        SeedSetConfidence, SeedSet, nonSeedSet = BuildGraphNetX.getSeedSet(
            DG_sbml,
            maxComponentSize=maxcc
        )

        # Remove any prefixes-suffixes
        SeedSetConfidence, SeedSet, nonSeedSet = (
            self._strip_pre_suff_from_dict(SeedSetConfidence),
            self._strip_pre_suff_from_list(SeedSet),
            self._strip_pre_suff_from_list(nonSeedSet)
        )

        # If carveme, map compounds to modelseed
        if self.namespace == "BiGG":

            seedSetBigg = SeedSet.copy()
            nonSeedSetBigg = nonSeedSet.copy()
            SeedSetConfidenceBigg = SeedSetConfidence.copy()

            SeedSet           = _bigg_to_modelseed(seedSetBigg, self.bigg2seed)
            nonSeedSet        = _bigg_to_modelseed(nonSeedSetBigg, self.bigg2seed)
            SeedSetConfidence = _bigg_to_modelseed(SeedSetConfidenceBigg, self.bigg2seed)

        return sbml_base, list(SeedSet), nonSeedSet, SeedSetConfidence


    def _strip_pre_suff_from_list(self, terms):
        return [
            term.split("_", 1)[-1] if term.startswith(self.compound_prefix) else term
            for term in (t.rsplit("_", 1)[0] if t.split("_")[-1] in {self.ex_suffix, self.int_suffix} else t for t in terms)
        ]

    # TODO (Haris Zafeiropoulos, 2025-03-28): check if this could be a static
    def _strip_pre_suff_from_dict(self, d):
        d_tmp = {}
        for k,v in d.items():
            new_k = self._strip_pre_suff_from_list([k])[0]
            d_tmp[new_k] = v
        return d_tmp

    # TODO (Haris Zafeiropoulos, 2025-03-28): like above
    def _serialize_dic(self, dict, json_file):

        dict_serial    = {k: list(v) for k, v in dict.items()}
        with open(json_file, 'w') as out_file:
            json.dump(dict_serial, out_file)
        return dict_serial


    def _dict_to_pickle(self, dict, pickle_file):

        dict_tmp = {}
        for k,v in dict.items():
            dict_tmp[k] = [
                kegg_module_related_intersect(v, self.modelseed_compounds_of_interest)
            ]
        df = pd.DataFrame.from_dict(dict_tmp)
        with open(pickle_file, "wb") as f:
            pickle.dump(df.T, f)




def kegg_module_related_intersect(intersect, modelseed_compounds_of_interest):
    """Check if KEGG MODULE related"""
    intersect = list(intersect)
    tmp_intersect = intersect.copy()
    for compl in tmp_intersect:
        if compl not in modelseed_compounds_of_interest:
            intersect.remove(compl)
    return intersect


def get_kegg_module_related(seed_ko_mo):
    """Load map file with KEGG modules and their terms and return a set with all the KOs there"""
    modules_compounds = pd.read_csv(seed_ko_mo, sep="\t")
    modules_compounds.columns = ["modelseed", "kegg", "module"]
    return set(modules_compounds["modelseed"].unique().tolist())


def _bigg_to_modelseed(bigg_obj, bigg2seed):
    """
    bigg2seed (pd.DataFrame)
    """
    if isinstance(bigg_obj, list):
        modelseed_seed_list = []
        for seed_BiggId in bigg_obj:
            if seed_BiggId not in bigg2seed:
                # print("not found: %s", seed_BiggId)  # This should be almost impossible..
                continue
            else:
                # NOTE (Haris Zafeiropoulos, 2025-03-24):
                # The bigg2seed dictionary has lists for values; if more than 1 modelseed compounds hit to the same BiGG
                # we keep the one with the lowest cpd since it's probably more involved to key processes
                mapped_ids = bigg2seed[seed_BiggId]
                mapped_ids = [x for x in mapped_ids if not isinstance(x, float)]
                if len(mapped_ids) > 0:
                    modelseed_ids = sorted(mapped_ids)
                    modelseed_seed_list.append(modelseed_ids[0])
                else:
                    # print("not found: %s", seed_BiggId)
                    modelseed_seed_list.append(seed_BiggId)
        return modelseed_seed_list

    elif isinstance(bigg_obj, dict):
        modelseed_seed_dict = {}
        for seed_BiggId, value in bigg_obj.items():
            if seed_BiggId not in bigg2seed:
                # logging.warning("not found in dictionary case: %s", seed_BiggId)
                continue
            else:
                mapped_ids = bigg2seed[seed_BiggId]
                mapped_ids = [x for x in mapped_ids if not isinstance(x, float)]
                if len(mapped_ids) > 0:
                    modelseed_ids = sorted(mapped_ids)
                    modelseed_seed_dict[modelseed_ids[0]] = value
                else:
                    # print("not found in dictionary case: %s", seed_BiggId)
                    modelseed_seed_dict[seed_BiggId] = value
        return modelseed_seed_dict


def progress_tracker(queue, total):
    """Progress bar updater."""
    with tqdm(total=total, desc="Calculate cooperation and competition scores as well as complementarities.") as pbar:
        for _ in range(total):
            queue.get()  # Wait for a task to finish
            pbar.update(1)


def generate_fixed_pairwise_comparisons(fixed_item, reconstruction_filenames):
    """Generate and return two lists: one with the fixed item in the first position and one with it in the second."""

    fixed_seedset_as_A = set()
    fixed_nonseedset_as_A = set()

    # Generate pairs where fixed_item is in the first position
    for B in reconstruction_filenames:
        fixed_seedset_as_A.add((fixed_item, B))

    # Generate pairs where fixed_item is in the second position
    for A in reconstruction_filenames:
        fixed_nonseedset_as_A.add((A, fixed_item))

    return list(fixed_seedset_as_A), list(fixed_nonseedset_as_A)


def bigg_to_seed_mapping_df(metanetx_compounds):

    # Open the tar.gz file
    with tarfile.open(metanetx_compounds, "r:gz") as tar:

        # List files in the archive to identify the one you want to read
        file_names = tar.getnames()  # Returns a list of files in the tar.gz

        # Extract the file of interest as a file-like object
        file_to_read = tar.extractfile(file_names[0])

        if file_to_read:
            metanetx = pd.read_csv(file_to_read, delimiter="\t", skiprows=353, header=None)  # Adjust delimiter as needed

    metanetx.columns = ["source", "id", "description"]

    metanetx[['source_namespace', 'source_id']] = metanetx['source'].str.split(pat=':', n=1, expand=True)
    metanetx.drop(columns=['source'], inplace=True)
    metanetx = metanetx[metanetx['source_namespace'].isin(["bigg.metabolite", "seed.compound"])]

    metanetx_bigg_metabolite = metanetx[metanetx["source_namespace"] == "bigg.metabolite"]
    metanetx_seed_compound   = metanetx[metanetx["source_namespace"] == "seed.compound"]

    merged_df = pd.merge(metanetx_bigg_metabolite, metanetx_seed_compound, on="id", how="left")
    bigg2seed = merged_df.groupby("source_id_x")["source_id_y"].apply(list).to_dict()

    return bigg2seed


# ---- Utils not related to the extaction but to the assignment of the complements or scores to the net

def load_seed_complement_files(path_to_kegg_seed_mappings):
    """
    Loads mapping files to be used for the building of the cx2 network.

    """
    kmap = pd.read_csv(os.path.join(path_to_kegg_seed_mappings, "seedId_keggId_module.tsv"), sep="\t", header=None)
    kmap.columns = ["modelseed", "kegg_compound", "kegg_module"]

    module_to_map           = pd.read_csv(os.path.join(path_to_kegg_seed_mappings, "module_map_pairs.tsv"), sep="\t", header=None)
    module_to_map.columns   = ["module", "map"]
    module_to_map['module'] = module_to_map['module'].str.replace("md:", '')
    module_to_map['map']    = module_to_map["map"].str.strip()

    module_map_dict = module_to_map.set_index('module')['map'].to_dict()
    kmap['map']     = kmap['kegg_module'].map(module_map_dict)
    maps_cat_descrs = pd.read_csv(
        os.path.join(path_to_kegg_seed_mappings, "related_kegg_maps_descriptions.tsv"),
        sep="\t",
        header=None
    )
    maps_cat_descrs.columns = ["map", "description", "category"]
    kmap                    = pd.merge(kmap, maps_cat_descrs, on='map', how='left')

    return kmap


def build_url_with_seed_complements(seed_complements, nonseeds, kmap, shortener=None):
    """

    """
    base_url = "https://www.kegg.jp/kegg-bin/show_pathway?"
    url      = "".join([base_url, kmap]) + "/"

    present_compounds_color   = "%20skyblue%2Cblue/"
    complemet_compounds_color = "%09%23ff0000/"

    for compound in nonseeds:
        url += compound + present_compounds_color
    for compound in seed_complements:
        url += compound + complemet_compounds_color
    if shortener is not None:
        logger.info("Shortening the URL.")
        url =  shortener.tinyurl.short(url)
    return url
