import os
import json
import shutil
import pickle
import h5py
import logging
import tarfile
import multiprocessing
import pandas as pd
from tqdm import tqdm
# from joblib import Parallel, delayed

from .PhyloMint.lib import BuildGraphNetX
from .PhyloMint.lib.CalculateIndexes import microbetagPI



class ExportSeedComplementarities():
    """
    Class to  export seed complements.
    Needs a config object to initiate it.

    # conda activate microbetag
    import yaml
    from config import Config
    from utils import ExportSeedComplementarities

    config_file = "tests/dev_io_microbetag/config.yml"
    with open(config_file, 'r') as yaml_file:
        config = Config(yaml.safe_load(yaml_file), config_file)

    seed_complements = ExportSeedComplementarities(config)
    """
    def __init__(self, config):
        self.seeds = config.seeds
        self.seed_sets = os.path.join(config.seeds, "SeedSetDic.json")
        self.non_seed_sets = os.path.join(config.seeds, "nonSeedSetDic.json")
        self.genres = config.genres
        self.logfile = os.path.join(config.seeds, "log.tsv")
        self.seed_ko_mo = config.seed_ko_mo
        self.module_seeds = os.path.join(self.seeds, "module_related_seeds.pckl")
        self.module_non_seeds = os.path.join(self.seeds, "module_related_non_seeds.pckl")
        self.seed_complements = os.path.join(self.seeds, "seed_complements.pckl")
        self.metanetx_compounds = config.metanetx_compounds
        self.genre_reconstruction_with = config.genre_reconstruction_with
        self.ex_suffix = "e" if self.genre_reconstruction_with == "carveme" else "e0"
        self.int_suffix = "c" if self.genre_reconstruction_with == "carveme" else "c0"
        # In case of carveme
        self.compound_prefix = "M"






    def export_seed_complements(self):
        """
        Export pairwise seed complmenents.
        Returns a df where beneficiary species are in the rows and potential donors in the columns.

        example:
        BIN                                             bin101-contigs                                     bin151-contigs                                      bin19-contigs                                     bin189-contigs
        BIN
        bin101-contigs                                                 []  [cpd02678, cpd00094, cpd02893, cpd00641, cpd00...  [cpd00259, cpd02678, cpd00641, cpd00200, cpd00...  [cpd02678, cpd00641, cpd00200, cpd00142, cpd00...
        bin151-contigs  [cpd03049, cpd00239, cpd03831, cpd11466, cpd00...                                                 []  [cpd03049, cpd00239, cpd03831, cpd11466, cpd00...  [cpd03049, cpd00239, cpd03831, cpd11466, cpd00...
        bin19-contigs   [cpd01777, cpd00055, cpd00121, cpd00482, cpd00...  [cpd01777, cpd00055, cpd00121, cpd00338, cpd00...                                                 []  [cpd00145, cpd21480, cpd01777, cpd02160, cpd00...
        """



class PhylomintMGT:

    def __init__(self, config):

        print("Initiating PhylomintMGT..")

        self.dir_path  = config.genres
        self.outdir    = config.seeds
        self.outfile   = os.path.join(self.outdir, "phylomint_scores.tsv")
        self.save_dics =  True
        self.threads   = config.threads
        self.sets_only = config.sets_only
        self.skip_sets = config.skip_sets
        self.prev_conf = config.prev_conf
        self.prev_nonseeds = config.prev_nonseeds
        self.genre_reconstruction_with = config.genre_reconstruction_with
        self.module_seeds = os.path.join(self.outdir, "module_related_seeds.pckl")

        # If we are using
        if config.users_models and self.skip_sets is False:
            if len(os.listdir(config.genres)) != len(os.listdir(config.for_reconstructions)):
                genre_files = [
                os.path.join(config.for_reconstructions, file)
                        for file in os.listdir(config.for_reconstructions)
                ]
                for file in genre_files:
                    dest_path = os.path.join(config.genres, os.path.basename(file))
                    shutil.copy(file, dest_path)

       #
        if self.genre_reconstruction_with == "carveme":
            print("Load bigg2seed map...")
            # Open the tar.gz file
            with tarfile.open(config.metanetx_compounds, "r:gz") as tar:
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
            metanetx_seed_compound = metanetx[metanetx["source_namespace"] == "seed.compound"]
            merged_df = pd.merge(metanetx_bigg_metabolite, metanetx_seed_compound, on="id", how="left")
            self.bigg2seed = merged_df.groupby("source_id_x")["source_id_y"].apply(list).to_dict()

        print("Load KEGG MODULE related terms..")
        self.seed_ko_mo = config.seed_ko_mo
        modules_compounds = pd.read_csv(self.seed_ko_mo, sep="\t")
        modules_compounds.columns = ["modelseed", "kegg", "module"]
        self.modelseed_compounds_of_interest = set(modules_compounds["modelseed"].unique().tolist())

        if self.skip_sets:
            print("Load previously computed confidence scores and non-seed sets..")
            try:
                with open(self.prev_conf, "r") as f:
                    self.ConfidenceDic = json.load(f)
            except FileExistsError as e:
                raise e
            try:
                with open(self.prev_nonseeds, "r") as f:
                    self.nonSeedSetDic = json.load(f)
            except FileExistsError as e:
                raise e


        self.ex_suffix = "e" if self.genre_reconstruction_with == "carveme" else "e0"
        self.int_suffix = "c" if self.genre_reconstruction_with == "carveme" else "c0"
        self.compound_prefix = "M"


    def get_sets(self):

        """
        Get seed and non-seed sets for each model
        """

        # Build initial dictionaries
        SeedSetDic = dict()
        nonSeedSetDic = dict()
        ConfidenceDic = dict()

        # Get all XML files in directory
        print("Export seed and non seed sets....")

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

            tmp = {key: None for key in SeedSet}
            SeedSetDic[sbml_base] = tmp.keys()
            nonSeedSetDic[sbml_base] = nonSeedSet
            ConfidenceDic[sbml_base] = SeedSetConfidence

        pool.close()
        pool.join()

        print("Seed and non seed sets have been exported.")

        self.SeedSetDic, self.nonSeedSetDic, self.ConfidenceDic = SeedSetDic, nonSeedSetDic, ConfidenceDic

        if self.save_dics:
            SeedSetDic_serializable = {k: list(v) for k, v in SeedSetDic.items()}
            with open(f'{self.outdir}/SeedSetDic.json', 'w') as out_file:
                json.dump(SeedSetDic_serializable, out_file)
            with open(f'{self.outdir}/nonSeedSetDic.json', 'w') as out_file:
                json.dump(nonSeedSetDic, out_file)
            with open(f'{self.outdir}/confidenceDic.json', 'w') as out_file:
                json.dump(ConfidenceDic, out_file)


    def get_scores(self):
        """
        Based on the seed and non-seed sets calculated, get all pairwise competition and cooperation scores and the seed complementarities
        between the modles under study.
        """
        # TODO (Haris Zafeiropoulos, 2025-03-24):
        # what if I do have the scores and i only need the complements... (funny but you never know)

        total_species = self.ConfidenceDic.keys()

        lock = multiprocessing.Lock()
        queue = multiprocessing.Queue()
        processes = []

        # for shared dict ----------------------------
        manager = multiprocessing.Manager()
        shared_dict = manager.dict()  # Shared dictionary for DataFrame

        # Start a separate process for tracking progress
        progress_process = multiprocessing.Process(target=progress_tracker, args=(queue, len(total_species)))
        progress_process.start()

        # Manually spawn processes with a limited number of concurrent threads
        for i, species in enumerate(total_species):
            if i % self.threads == 0:
                for process in processes:
                    process.join()  # Wait for the batch to finish
                processes = []  # Clear completed processes

            process = multiprocessing.Process(target=self.worker_function, args=(lock, species, queue, shared_dict))
            processes.append(process)
            process.start()

        # Wait for the remaining processes to finish
        for process in processes:
            process.join()

        # Signal the progress tracker to stop
        queue.put(None)
        progress_process.join()

        # Ensure inner dictionaries are not lost
        shared_dict_serializable = {k: dict(v) for k, v in shared_dict.items()}

        with open("asd.json", "w") as f:
            json.dump(shared_dict_serializable, f)  # Pretty print for readability

        df = pd.DataFrame.from_dict(shared_dict_serializable)
        with open("asd.pckl", "wb") as f:
            pickle.dump(df, f)

        """
        import pandas as pd
        # Load only rows where "Source" is "GPB:bin_000082"
        df_partial = pd.read_hdf("module_related_seeds.h5", key="seeds", where='Source == "GPB:bin_000082"')

        print(df_partial)
        """



    def scores_and_overlaps_for_a_species(self, species, lock, shared_dict):
        """
        Get scores and complements for a specific model (species)
        """
        print(species)

        # Get pairs with species as A and species as B
        as_beneficiary, _ = generate_fixed_pairwise_comparisons(species, list(self.ConfidenceDic.keys()))

        # Get species seed and non-seed sets
        species_seedset = self.ConfidenceDic[species]

        # Species as A
        compls = {}
        findings = set()
        for pair in as_beneficiary:
            partner = pair[1]
            if partner == species:
                continue

            species = species.replace(".PATRIC", "")
            partner = partner.replace(".PATRIC", "")

            SeedSetBConfidence, nonSeedB = self.ConfidenceDic[partner], self.nonSeedSetDic[partner]
            MetabolicCooperationIdxAB, MetabolicCompetitionIdxAB, B_complememts_to_A = microbetagPI(species_seedset, SeedSetBConfidence, nonSeedB)
            # Get only KEGG module - related
            kegg_module_related_B_complememts_to_A = self.kegg_module_related_intersect(B_complememts_to_A)
            # Line to print
            output = f"{species}\t{partner}\t{MetabolicCompetitionIdxAB}\t{MetabolicCooperationIdxAB}\n"
            compls[partner] = kegg_module_related_B_complememts_to_A
            findings.add(output)

        # Safely write to the file with a lock
        with lock:
            with open(self.outfile, 'a') as f:
                for r in findings:
                    f.write(r)
            # Load KEGG related complements to shared dict
            shared_dict[species] = compls


    def worker_function(self, lock, species, queue, shared_dict):
        """Wrapper function to process a species and signal completion."""
        self.scores_and_overlaps_for_a_species(species, lock, shared_dict)
        with lock:
            queue.put(1)  # Signal that one task is completed


    def kegg_module_related_intersect(self, intersect):
        """Check if KEGG MODULE related"""
        intersect = list(intersect)
        tmp_intersect = intersect.copy()
        for compl in tmp_intersect:
            if compl not in self.modelseed_compounds_of_interest:
                intersect.remove(compl)
        return intersect


    def process_sbml(self, sbml_path, maxcc=2):

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
        if self.genre_reconstruction_with == "carveme":

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

    def _strip_pre_suff_from_dict(self, d):
        d_tmp = {}
        for k,v in d.items():
            new_k = self._strip_pre_suff_from_list([k])[0]
            d_tmp[new_k] = v
        return d_tmp


def _bigg_to_modelseed(bigg_obj, bigg2seed):
    """
    bigg2seed (pd.DataFrame)
    """
    if isinstance(bigg_obj, list):
        modelseed_seed_list = []
        for seed_BiggId in bigg_obj:
            if seed_BiggId not in bigg2seed:
                print("not found: %s", seed_BiggId)  # This should be almost impossible..
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
                    print("not found: %s", seed_BiggId)
                    modelseed_seed_list.append(seed_BiggId)
        return modelseed_seed_list

    elif isinstance(bigg_obj, dict):
        modelseed_seed_dict = {}
        for seed_BiggId, value in bigg_obj.items():
            if seed_BiggId not in bigg2seed:
                logging.warning("not found in dictionary case: %s", seed_BiggId)
            else:
                mapped_ids = bigg2seed[seed_BiggId]
                mapped_ids = [x for x in mapped_ids if not isinstance(x, float)]
                if len(mapped_ids) > 0:
                    modelseed_ids = sorted(mapped_ids)
                    modelseed_seed_dict[modelseed_ids[0]] = value
                else:
                    print("not found in dictionary case: %s", seed_BiggId)
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



def load_seed_complement_files(path_to_kegg_seed_mappings):
    """

    """
    kmap = pd.read_csv(os.path.join(path_to_kegg_seed_mappings, "seedId_keggId_module.tsv"), sep="\t", header=None)
    kmap.columns = ["modelseed", "kegg_compound", "kegg_module"]

    module_to_map = pd.read_csv(os.path.join(path_to_kegg_seed_mappings, "module_map_pairs.tsv"), sep="\t", header=None)
    module_to_map.columns = ["module", "map"]
    module_to_map['module'] = module_to_map['module'].str.replace("md:", '')
    module_to_map['map'] = module_to_map["map"].str.strip()

    module_map_dict = module_to_map.set_index('module')['map'].to_dict()
    kmap['map'] = kmap['kegg_module'].map(module_map_dict)

    maps_categories_and_descrs = pd.read_csv(os.path.join(path_to_kegg_seed_mappings, "related_kegg_maps_descriptions.tsv"), sep="\t", header=None)
    maps_categories_and_descrs.columns = ["map", "description", "category"]

    kmap = pd.merge(kmap, maps_categories_and_descrs, on='map', how='left')

    return kmap


def order_seed_complements(r):
    """
    Order seed complements so they display based on their metabolism category
    which have been ranked according to what metabolic interactions we believe most common.
    """
    # Define a custom sorting function
    def custom_sort(item):
        return category_index.get(item[0], len(order_list))

    # Create a dictionary to map each category to its corresponding index in the order_list
    order_list = [
        'Amino acid metabolism',
        'Metabolism of cofactors and vitamins',
        'Energy metabolism',
        'Carbohydrate metabolism',
        'Nucleotide metabolism',
        'Biosynthesis of other secondary metabolites',
        'Biosynthesis of terpenoids and polyketides',
        'Lipid metabolism',
        'Glycan metabolism',
        'Xenobiotics biodegradation'
    ]
    category_index = {category: index for index, category in enumerate(order_list)}

    # Sort the data using the custom sorting function
    sorted_data = sorted(r, key=custom_sort)
    return sorted_data


def build_url_with_seed_complements(seed_complements, nonseeds, kmap, shortener=None):
    """

    """
    base_url = "https://www.kegg.jp/kegg-bin/show_pathway?"
    url = "".join([base_url, kmap]) + "/"

    present_compounds_color = "%20skyblue%2Cblue/"
    complemet_compounds_color = "%09%23ff0000/"

    for compound in nonseeds:
        url += compound + present_compounds_color
    for compound in seed_complements:
        url += compound + complemet_compounds_color
    if shortener is not None:
        logging.info("Shortening the URL.")
        url =  shortener.tinyurl.short(url)
    return url

