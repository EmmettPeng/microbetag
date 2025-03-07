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
    seed_complements.update()
    """
    def __init__(self, config):
        self.seeds = config.seeds
        self.seed_sets = os.path.join(config.seeds, "SeedSetDic.json")
        self.nonseed_sets = os.path.join(config.seeds, "nonSeedSetDic.json")
        self.genres = config.genres
        self.logfile = os.path.join(config.seeds, "log.tsv")
        self.updated_seed_sets = os.path.join(config.seeds, "updated_SeedsDic.json")
        self.updated_non_seed_sets = os.path.join(config.seeds, "updated_nonSeedsDic.json")
        self.seed_ko_mo = config.seed_ko_mo
        self.module_seeds = os.path.join(self.seeds, "module_related_seeds.pckl")
        self.module_non_seeds = os.path.join(self.seeds, "module_related_non_seeds.pckl")
        self.seed_complements = os.path.join(self.seeds, "seed_complements.pckl")
        self.metanetx_compounds = config.metanetx_compounds
        self.genre_reconstruction_with = config.genre_reconstruction_with

        self.ex_suffix = "_e" if self.genre_reconstruction_with == "carveme" else "_e0"
        self.int_suffix = "_c" if self.genre_reconstruction_with == "carveme" else "_c0"

        # In case of carveme
        self.periplasm = "_p"
        self.compound_prefix = "M_"

        if config.users_models:
            if len(os.listdir(config.genres)) != len(os.listdir(config.for_reconstructions)):
                genre_files = [
                os.path.join(config.for_reconstructions, file)
                        for file in os.listdir(config.for_reconstructions)
                ]
                for file in genre_files:
                    dest_path = os.path.join(config.genres, os.path.basename(file))
                    shutil.copy(file, dest_path)

    def update(self):
        """
        PhyloMInt does not consider the
        Update seed and non-seed sets returned by PhyloMint by:
        - removing compounds from seed sets that are related to environmental metabolites that can be produced in several ways within the cell.
        - removing from non-seed sets compounds that cannot be produced in any other way than from entering the cell from the environment.

        If both _c0 and _e0 exist in the seed list, _e0 is kept.
        If _c0 is missing from the model, _e0 is kept as a seed.
        If _c0 can be produced without _e0, _e0 is not a seed.
        Otherwise, _e0 is kept as a seed and _c0 is removed from non-seeds.

        """
        updated_seeds = {}
        updated_nonSeeds = {}
        current_seeds = json.load(open(self.seed_sets, "r"))
        current_nonSeeds = json.load(open(self.nonseed_sets, "r"))

        for xml in os.listdir(self.genres):

            s1 = time.time()

            model_id, _ = os.path.splitext(xml)
            model       = cobra.io.read_sbml_model(os.path.join(self.genres, xml))

            tmp_seeds     = [x.lstrip(self.compound_prefix) for x in current_seeds[model_id]]
            tmp_seeds_main = set(["_".join(pot_seed.split("_")[:-1]) for pot_seed in tmp_seeds])

            tmp_non_seeds = [x.lstrip(self.compound_prefix) for x in current_nonSeeds[model_id]]
            tmp_non_seeds_main = set(["_".join(pot_nonseed.split("_")[:-1]) for pot_nonseed in tmp_non_seeds])

            models_non_seeds = tmp_non_seeds_main.copy()
            model_seeds = set()

            # Check if we keep intra- or extracellular compound
            for seed_id in tmp_seeds_main:

                ex = seed_id + self.ex_suffix
                cel = seed_id + self.int_suffix
                per = seed_id + self.periplasm

                # NOTE (2025-03-07):
                # PhyloMInt does not consider the fact that in a metabolic model, you may have the same metabolites across different compartments
                reactions = model.reactions - model.exchanges
                check = [
                    rxn.id for rxn in reactions
                    if (
                        # TODO (2025-03-07): Should we remove this the reversibility check ?
                        # not rxn.reversibility and
                        any(m.id == cel for m in rxn.products) and
                        all(m.id not in {ex, per} for m in rxn.reactants + rxn.products)
                    )
                ]

                if len(check) == 0:
                    model_seeds.add(seed_id)

            for nonseed_id in tmp_non_seeds_main:
                if nonseed_id in  model_seeds:
                    models_non_seeds.discard(nonseed_id)

            print("original # of seeds:", len(tmp_seeds), "updated # of seeds:", len(model_seeds) )
            print("original # of non-seeds:", len(tmp_non_seeds), "updated # of non-seeds:", len(models_non_seeds) )

            # Update dics
            updated_seeds[model_id] = list(model_seeds)
            updated_nonSeeds[model_id] = list(models_non_seeds)

            s2 = time.time()
            logging.info(f"{s2-s1} seconds for a .xm_tags load")

        logging.info("Update function is done and about to save updated json files.")

        with open(self.updated_seed_sets, "w") as f:
            json.dump(convert_to_json_serializable(updated_seeds), f)

        with open(self.updated_non_seed_sets, "w") as f:
            json.dump(convert_to_json_serializable(updated_nonSeeds), f)

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
        non_seedset_file = json.load(open(self.updated_non_seed_sets,"r"))
        model_names = list(non_seedset_file.keys())

        for smodel_name in model_names:
            non_seedset = set( [x[2:] for x in non_seedset_file[smodel_name]] )
            non_seedset_no_compartments = set( [x[2:].rsplit("_", 1)[0] for x in non_seedset_file[smodel_name] ])
            non_seeds_of_interest = non_seedset_no_compartments.intersection(modelseed_compounds_of_interest)
            mean_non_seedset_length += len(non_seedset)
            patricId_to_non_seeds_of_interest[smodel_name] = non_seeds_of_interest
            mean_non_seedset_of_interest += len(non_seeds_of_interest)

        mean_seedset_length = 0
        mean_seedset_of_interest = 0
        seedset_file = json.load(open(self.updated_seed_sets, "r"))
        model_names = list(seedset_file.keys())

        for model_name in model_names:
            number_of_models += 1
            # Get seeds with and without their compartment specific part
            seedset = set([x[2:]  for x in seedset_file[model_name]])
            seedset_no_compartments = set([
                x[2:].rsplit("_", 1)[0]  for x in seedset_file[model_name]
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
        with open(self.updated_seed_sets) as f:
            updated_seeds = json.load(f)
        with open(self.updated_non_seed_sets) as f:
            updated_non_seeds = json.load(f)
        bigg_seeds = os.path.join(self.seeds, "updatedBiggSeedsDic.json")
        bigg_non_seeds = os.path.join(self.seeds, "updatedBiggNonSeedsDic.json")

        shutil.move(self.updated_seed_sets, bigg_seeds)
        shutil.move(self.updated_non_seed_sets, bigg_non_seeds)

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
        print("\n\n>>>>", bigg2seed)

        updated_seeds_biggIds, _     = process_carve_seeds(updated_seeds, bigg2seed, self.int_suffix)
        updated_non_seeds_biggIds, _ = process_carve_seeds(updated_non_seeds, bigg2seed, self.int_suffix)

        with open(self.updated_seed_sets, "w") as f:
            json.dump(updated_seeds_biggIds, f)
        with open(self.updated_non_seed_sets, "w") as f:
            json.dump(updated_non_seeds_biggIds, f)



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
            "".join(["M_", item, int_suffix])
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

