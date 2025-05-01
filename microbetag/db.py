import os
import json
import pandas as pd
import mysql.connector
from mysql.connector import pooling

from .utils import mtg_logger

logger = mtg_logger(__file__)
_KEGG_MAPPINGS = os.path.join(
    os.path.dirname(__file__), "mtg_maps_models", "kegg_mappings"
)
# -----------
# Database
# -----------

script_dir = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(script_dir, ".env_dev.json")) as f:
    c = json.load(f)
    config = {
        "user": c["Luna"]["USER_NAME"],
        "password": c["Luna"]["PASSWORD"],
        "host": c["Luna"]["HOST"],      # 'localhost'
        "port": c["Luna"]["PORT"],
        "database": c["Luna"]["DB_NAME"],
    }


def execute(phrase):
    """
    Establish a database connection and perform an action
    """
    # Database connection configuration
    # [TODO] Switch to "db" when running on the container
    cnx    = mysql.connector.connect(**config)
    cursor = cnx.cursor()
    cursor.execute(phrase)
    rows = cursor.fetchall()
    cursor.close()
    cnx.commit()
    cnx.close()
    return rows


def execute_in_a_pool(cursor, query):
    """
    Executes a query using a connection from the connection pool.
    """
    try:
        cursor.execute(query)
        result = [row for row in cursor]
        return result
    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))
        print(query)


def init_connection_pool():
    """
    Initiates a connection pool.
    A pool opens a number of connections and handles thread safety when providing connections to requesters.
    For more see: https://dev.mysql.com/doc/connector-python/en/connector-python-connection-pooling.html
    """
    connection_pool = pooling.MySQLConnectionPool(
        pool_name="microbetagDB_pool",
        pool_size=5,
        user=config["user"],
        password=config["password"],
        host=config["host"],
        database=config["database"],
    )
    return connection_pool


# --------
# Mapping
# --------
def gc_unify(gc_list):
    """
    Removes duplicates of a genome that has entries both as GCA and GCF in the db.
    """
    return list(set([x.replace("GCA", "GCF") for x in gc_list]))


def alt_genome_prefix(gc):
    """
    Switches GCA prefic of a genome accession id to GCF and vice-versa.
    """
    if gc.startswith("GCA_"):
        gc_alt = gc.replace("GCA_", "GCF_")

    elif gc.startswith("GCF_"):
        gc_alt = gc.replace("GCF_", "GCA_")

    else:
        return 0
    return gc_alt


def get_genomes_for_ncbi_tax_id(ncbi_tax_id=1281578):
    """
    Get the genome IDs corresponding to a given NCBI Taxonomy ID from the microbetagDB.
    """
    query      = f"SELECT genomeId FROM genome2taxNcbiId WHERE ncbiTaxId = {ncbi_tax_id};"
    genome_ids = execute(query)

    return {ncbi_tax_id: gc_unify([x[0] for x in genome_ids])}


def get_ncbi_tax_if_for_genome(gc_id="GCA_018819265.1"):
    #
    query      = f"SELECT ncbiTaxId FROM genome2taxNcbiId WHERE genomeId = '{gc_id}';"
    ncbi_tax_ids = execute(query)
    #
    return {gc_id: list(ncbi_tax_ids)}


def get_patric_id_of_gc_accession_list(gc_accession_list=["GCA_003184265.1"]):
    """
    Gets a list of GC accession ids and returns a dictionary where the GC ids are the keys
    and their corresponding PATRIC ids are the values.
    """
    gc_to_patric_dict = {}
    for gc in gc_accession_list:
        gc_alt = alt_genome_prefix(gc)
        if gc_alt != 0:
            query = f"SELECT patricId FROM patricId2genomeId WHERE gtdbGenomeAccession = '{gc}' OR gtdbGenomeAccession = '{gc_alt}';"
            patricId = execute(query)
            gc_to_patric_dict[gc] = patricId[0][0] if patricId and patricId[0] else None
    logger.info(gc_to_patric_dict)
    return gc_to_patric_dict


# --------
# Phen related
# --------
def get_phen_traits(config):

    traits_per_genome = {
        genome: get_phendb_traits(genome)
        for genome in config.repr_genomes_present
    }

    df = pd.DataFrame(traits_per_genome)

    # Drop 'gtdbId' only if it exists in the index to avoid errors
    df = df.drop("gtdbId", errors="ignore")  # drop to make it safe if "gtdbId" isn't present.

    # all_traits_file = os.path.join(config.predictions_path, "all_predictions.tsv")
    # df.to_csv(all_traits_file, sep="\t")

    write_trait_file(df, config.predictions_path)


def phen_query(gtdb_id):

    return "".join([
        "SELECT * FROM phenDB WHERE SUBSTRING_INDEX(gtdbId, '.', 1) = SUBSTRING_INDEX('",
        gtdb_id,
        "', '.', 1);"
    ])


def get_phendb_traits(gtdb_genome_id="GCA_018819265.1"):
    """
    Get phenotypical traits based on phenDB classes based on its GTDB representative genome
    """
    gtdb_id = gtdb_genome_id.strip()

    rows = execute(phen_query(gtdb_id))

    if len(rows) == 0:
        alt_gtdb_id = alt_genome_prefix(gtdb_genome_id)
        rows = execute(phen_query(alt_gtdb_id))

    if len(rows) == 0:
        logger.info(f"Genome {gtdb_genome_id} is not a NCBI accession id. It could be a MGnify or a KEGG one.")
        return 0

    query_colnames = "SHOW COLUMNS FROM phenDB;"
    colnames = [list(x)[0] for x in execute(query_colnames)]
    genomes_traits = {i: j for i, j in zip(colnames, rows[0])}

    return genomes_traits  # gtdb_genome_id


def write_trait_file(df, output_dir=None):

    if output_dir is None:
        output_dir = os.getcwd()

    # Loop through traits (every 2 rows)
    for i in range(0, df.shape[0], 2):
        trait    = df.index[i]
        presence = df.iloc[i]
        score    = df.iloc[i + 1]

        output = pd.DataFrame({
            'Identifier': presence.index,
            'Trait present': presence.values,
            'Confidence': score.values
        })

        # Write to .tsv
        trait    = f"{trait}.prediction.tsv"
        filename = os.path.join(output_dir, trait)
        with open(filename, "w") as f:
            f.write(f"# Trait: {trait}\n")
            f.write("Identifier\tTrait present\tConfidence\n")
            output.to_csv(f, sep="\t", index=False, header=False)


# --------
# Pathway complementarity
# --------

def get_complements_of_list_of_pair_of_ncbiIds(relative_genomes, pairs_of_interest):
    """
    pairs_of_interest={("1260918", "1819566")}
    relative_genomes={
        "1260918": {"GCF_002102185.1"},
        "1819566": {"GCF_009711525.1"}
    }
    """

    logger.info("===> Building queries for genome pairs...")
    complements_ids_queries = build_complement_queries(relative_genomes, pairs_of_interest)

    logger.info("===> Executing unique queries...")
    unique_queries = {q for genome_pairs in complements_ids_queries.values() for q in genome_pairs.values()}
    unique_queries2comples = get_complement_ids(unique_queries)

    logger.info("===> Mapping complement IDs to genome pairs...")
    pairs_to_compl_ids = map_queries_to_pairs(complements_ids_queries, unique_queries2comples)

    logger.info("===> Fetching full complement metadata...")
    all_compl_ids2coloured_compls = get_coloured_complements(pairs_to_compl_ids)

    logger.info("===> Assembling final result...")
    return build_pairs_complements(pairs_to_compl_ids, all_compl_ids2coloured_compls)


def query_for_getting_compl_ids(beneficiary="GCA_003184265.1", donor="GCA_000015645.1"):
    """
    Gets 2 gc accession ids and returns a query for their pathway complementarities
    """
    beneficiary_alt = alt_genome_prefix(beneficiary)
    donor_alt       = alt_genome_prefix(donor)

    return "".join(
        [
            'SELECT complmentId FROM pathwayComplementarity WHERE (beneficiaryGenome = "',
            str(beneficiary),
            '" or beneficiaryGenome = "',
            str(beneficiary_alt),
            '") AND (donorGenome = "',
            str(donor),
            '" OR donorGenome = "',
            str(donor_alt),
            '");',
        ]
    )


def build_complement_queries(relative_genomes, pairs_of_interest):

    complements_ids_queries = {}

    for ncbi_a, ncbi_b in pairs_of_interest:
        for genome_a in relative_genomes[ncbi_a]:
            for genome_b in relative_genomes[ncbi_b]:
                query = query_for_getting_compl_ids(genome_a, genome_b)
                complements_ids_queries.setdefault((ncbi_a, ncbi_b), {})[(genome_a, genome_b)] = query

    return complements_ids_queries


def get_complement_ids(unique_queries):
    cnx_pool      = init_connection_pool()
    my_connection = cnx_pool.get_connection()
    cursor        = my_connection.cursor()

    query_results = {}
    for query in unique_queries:
        result = execute_in_a_pool(cursor, query)
        if result:
            query_results[query] = result[0][0].split(",")
    return query_results


def map_queries_to_pairs(complements_ids_queries, unique_queries2comples):
    pairs_to_compl_ids = {}
    for ncbi_pair, genome_query_map in complements_ids_queries.items():
        for genome_pair, query in genome_query_map.items():
            if query in unique_queries2comples:
                pairs_to_compl_ids.setdefault(ncbi_pair, {})[genome_pair] = unique_queries2comples[query]
    return pairs_to_compl_ids


def get_coloured_complements(pairs_to_compl_ids):
    unique_compl_ids = {
        compl_id
        for genome_map in pairs_to_compl_ids.values()
        for compl_list in genome_map.values()
        for compl_id in compl_list
    }
    unique_compl_ids = list(unique_compl_ids)
    if not unique_compl_ids:
        return {}

    query = """SELECT KoModuleId, complement, pathway FROM uniqueComplements
               WHERE complementId IN ('{}') ORDER BY FIELD(complementId, '{}');""".format(
        "','".join(unique_compl_ids), "','".join(unique_compl_ids)
    )

    cnx_pool                 = init_connection_pool()
    my_connection            = cnx_pool.get_connection()
    cursor                   = my_connection.cursor()
    query_result             = execute_in_a_pool(cursor, query)
    query_result_list        = [[r] for r in query_result]
    colored_complements_list = build_kegg_urls(query_result_list)

    return dict(zip(unique_compl_ids, colored_complements_list))


def build_pairs_complements(pairs_to_compl_ids, all_compl_ids2coloured_compls):
    pairs_complements = {}
    for ncbi_pair, genome_map in pairs_to_compl_ids.items():
        for genome_pair, compl_ids in genome_map.items():
            for compl_id in compl_ids:
                colored = all_compl_ids2coloured_compls.get(compl_id)
                if colored:
                    pairs_complements.setdefault(ncbi_pair, {}).setdefault(genome_pair, []).append(colored)
    return pairs_complements


# ====================================================

def build_kegg_urls(genome_pair_compls):
    """
    Takes as input the complements list between two genomes and
    build urls to colorify the related to the module kegg map based on the KO terms of the beneficiary (pink)
    and those it gets from the donor (green).
    NOTE: some modules do not belong to any map, e.g. https://www.kegg.jp/module/M00705. In these cases, we will have a N/A value in the url.
    """
    # Constants
    color_map_base_url   = "https://www.kegg.jp/kegg-bin/show_pathway?"
    present_kos_color    = "%09%23EAD1DC/"
    complement_kos_color = "%09%2300A898/"

    # Load module-to-map mappings
    module_map_file = os.path.join(_KEGG_MAPPINGS, "module_map_pairs.tsv")

    logger.info(module_map_file)

    with open(module_map_file, "r") as f:
        module_map = {
            line.split("\t")[0][3:]: line.split("\t")[1].strip()[1:-1]
            for line in f
        }

    updated_complements = []

    for [module_id, complement_str, ko_terms_str] in [list(c[0]) for c in genome_pair_compls]:

        beneficiary_kos    = []
        complement_kos     = []
        complement_kos_set = set(complement_str.split(";"))
        ko_terms           = ko_terms_str.split(";")

        for ko in ko_terms:
            if ko in complement_kos_set:
                complement_kos.append(ko + complement_kos_color)
            else:
                beneficiary_kos.append(ko + present_kos_color)

        # Construct URL or fallback
        url    = "N/A"
        map_id = module_map.get(module_id)
        if map_id:
            url = f"{color_map_base_url}{map_id}/" + "".join(beneficiary_kos + complement_kos)

        # Append new colored URL to the record
        updated_complements.append([[module_id, complement_str, ko_terms_str, url]])

    return updated_complements

# ====================================================
