import os
import json
import time
import cobra
import shutil
import pickle
import logging
import tarfile
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
from .utils import convert_to_json_serializable


def process_sets(file_path, compound_prefix, ex_suffix, int_suffix):
    with open(file_path) as f:
        sets = json.load(f)
    updated_sets = sets.copy()
    for k, v in sets.items():
        tmp_v = [x.split("_", 1)[-1] if x.split("_")[0] == compound_prefix else x for x in v]
        new_v = [x.rsplit("_", 1)[0] if x.split("_")[-1] in [ex_suffix, int_suffix] else x for x in tmp_v]
        updated_sets[k] = new_v
    return updated_sets


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

        if config.users_models:
            if len(os.listdir(config.genres)) != len(os.listdir(config.for_reconstructions)):
                genre_files = [
                os.path.join(config.for_reconstructions, file)
                        for file in os.listdir(config.for_reconstructions)
                ]
                for file in genre_files:
                    dest_path = os.path.join(config.genres, os.path.basename(file))
                    shutil.copy(file, dest_path)
        self.updated_seeds     = process_sets(self.seed_sets, self.compound_prefix, self.ex_suffix, self.int_suffix)
        self.updated_non_seeds = process_sets(self.non_seed_sets, self.compound_prefix, self.ex_suffix, self.int_suffix)

    def module_related_seeds(self):
        """
        Get seed and non seed sets with terms related to KEGG MODULES.
                                                       NonSeedSet
        BIN
        bin101-contigs  [cpd02817, cpd00344, cpd03123, cpd00482, cpd11...
        """
        modules_compounds = pd.read_csv(self.seed_ko_mo, sep="\t")
        modules_compounds.columns = ["modelseed", "kegg", "module"]
        modelseed_compounds_of_interest = set(modules_compounds["modelseed"].unique().tolist())

        number_of_models = 0
        patricId_to_seeds_of_interest = {}
        patricId_to_non_seeds_of_interest = {}
        mean_non_seedset_length = 0 ; mean_non_seedset_of_interest = 0
        # non_seedset_file = json.load(open(self.non_seed_sets,"r"))
        # model_names = list(non_seedset_file.keys())
        model_names = list(self.updated_non_seeds.keys())

        for smodel_name in model_names:
            non_seedset = set( [x[2:] for x in self.updated_non_seeds[smodel_name]] )
            non_seedset_no_compartments = set([
                x[2:].rsplit("_", 1)[0]
                for x in self.updated_non_seeds[smodel_name]
            ])
            non_seeds_of_interest = non_seedset_no_compartments.intersection(modelseed_compounds_of_interest)
            mean_non_seedset_length += len(non_seedset)
            patricId_to_non_seeds_of_interest[smodel_name] = non_seeds_of_interest
            mean_non_seedset_of_interest += len(non_seeds_of_interest)

        mean_seedset_length = 0
        mean_seedset_of_interest = 0
        # seedset_file = json.load(open(self.seed_sets, "r"))
        # model_names = list(seedset_file.keys())
        model_names = list(self.updated_seeds.keys())
        for model_name in model_names:
            number_of_models += 1
            # Get seeds with and without their compartment specific part
            seedset = set([x[2:]  for x in self.updated_seeds[model_name]])
            seedset_no_compartments = set([
                x[2:].rsplit("_", 1)[0]  for x in self.updated_seeds[model_name]
            ])
            mean_seedset_length += len(seedset)
            seeds_of_interest = seedset_no_compartments.intersection(modelseed_compounds_of_interest)
            seeds_of_interest_tmp = list(seeds_of_interest.copy())
            for pot_seed in seeds_of_interest:
                if pot_seed in patricId_to_non_seeds_of_interest[model_name]:
                    seeds_of_interest_tmp.remove(pot_seed)
            patricId_to_seeds_of_interest[model_name] = set(seeds_of_interest_tmp)
            mean_seedset_of_interest += len(set(seeds_of_interest_tmp))

        logging.info("Mean length of initial seedset: %s", str(mean_seedset_length/number_of_models))
        logging.info("Mean length of seedsets of interest: %s", str(mean_seedset_of_interest/number_of_models))
        logging.info("Mean of initial length of non seed sets: %s", str(mean_non_seedset_length/number_of_models))
        logging.info("Mean of non seed sets of interest: %s", str(mean_non_seedset_of_interest/number_of_models))

        tmp_dict = {key: list(value) for key, value in patricId_to_seeds_of_interest.items()}
        df1 = pd.DataFrame(list(tmp_dict.items()), columns=['BIN', 'SeedSet'])
        df1['BIN'] = df1['BIN'].str.replace('.BIN', '')
        df1.set_index('BIN', inplace=True)

        with open(self.module_seeds,"wb") as f:
            pickle.dump(df1, f)

        tmp_dict = {key: list(value) for key, value in patricId_to_non_seeds_of_interest.items()}
        df2 = pd.DataFrame(list(tmp_dict.items()), columns=['BIN', 'NonSeedSet'])
        df2['BIN'] = df2['BIN'].str.replace('.BIN', '')
        df2.set_index('BIN', inplace=True)

        with open( self.module_non_seeds,"wb") as f:
            pickle.dump(df2, f)

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
        # Function to calculate overlap
        def calculate_overlap(seed_set, non_seed_set):
            return list(set(seed_set) & set(non_seed_set))

        # Parallelized function to calculate overlap for one row in df1 with all rows in df2
        def calculate_overlap_parallel(row1, df2):
            return [calculate_overlap(row1['SeedSet'], row2['NonSeedSet']) for _, row2 in df2.iterrows()]

        # Load my case
        with open(self.module_seeds, "rb") as f:
            df1 = pickle.load(f)
        with open(self.module_non_seeds, "rb") as f:
            df2 = pickle.load(f)

        # Create a new DataFrame for overlaps
        overlaps_df = pd.DataFrame(index=df1.index, columns=df2.index)

        # Parallelize the overlap calculation using joblib with tqdm for progress tracking
        num_cores = -1  # Use all available cores

        results = Parallel(n_jobs=num_cores)(
            delayed(calculate_overlap_parallel)(row1, df2)
            for _, row1 in tqdm(df1.iterrows(),
                                total=len(df1)
                                )
            )

        # Fill in the overlaps DataFrame with calculated values
        for i, row in enumerate(results):
            overlaps_df.iloc[i] = row

        with open(self.seed_complements,"wb") as f:
            pickle.dump(overlaps_df, f)

        del results

    def map_carveme_seeds(self):
        """
        https://www.metanetx.org/mnxdoc/mnxref.html
        map to modelseed and/or kegg..
        [REMEMBER] MGG points to modelseed ids
        We will map
        """

        # Open the tar.gz file
        with tarfile.open(self.metanetx_compounds, "r:gz") as tar:
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
        bigg2seed = merged_df.groupby("source_id_x")["source_id_y"].apply(list).to_dict()

        seeds_biggIds, _     = process_carve_seeds(self.updated_seeds, bigg2seed, self.int_suffix)
        non_seeds_biggIds, _ = process_carve_seeds(self.updated_non_seeds, bigg2seed, self.int_suffix)

        bigg_seeds = os.path.join(self.seeds, "updatedBiggSeedsDic.json")
        bigg_non_seeds = os.path.join(self.seeds, "updatedBiggNonSeedsDic.json")

        with open(bigg_seeds, "w") as f:
            json.dump(seeds_biggIds, f)
        with open(bigg_non_seeds, "w") as f:
            json.dump(non_seeds_biggIds, f)

        self.updated_seeds = seeds_biggIds
        self.updated_non_seeds = non_seeds_biggIds


def process_carve_seeds(seeds_dict, bigg2seed, int_suffix):
    """

    bigg2seed (pd.DataFrame)
    """
    updated_biggIds = {}
    bigg_ids_not_mapped_to_seed = {}
    for bin_id, seeds in seeds_dict.items():
        updated_biggIds[bin_id] = []
        for seed_id in seeds:
            # seed_id_part = "_".join(seed_id.split("_")[1:-1])
            if seed_id not in bigg2seed:
                logging.warning("not found: %s", seed_id)
                bigg_ids_not_mapped_to_seed.setdefault(bin_id, []).append(seed_id)
                continue
            updated_biggIds[bin_id].append(bigg2seed[seed_id])
        flat_list = [
            "_".join(["M", item, int_suffix])
            for sublist in updated_biggIds[bin_id]
            if sublist
            for item in sublist
            if not pd.isna(item)
        ]
        updated_biggIds[bin_id] = flat_list
    return updated_biggIds, bigg_ids_not_mapped_to_seed


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

