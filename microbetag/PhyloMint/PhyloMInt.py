#!/usr/bin/env python3

# Built-in/Generic Imports
import sys
import os
import json
import itertools
import functools
import subprocess
from multiprocessing import Manager
from multiprocessing import Process
from argparse import ArgumentParser
from pathlib import Path

# Libraries
from libsbml import readSBML
import multiprocessing
from tqdm import tqdm

# Local modules
from .lib import BuildGraphNetX
from .lib.CalculateIndexes import microbetagPI

__author__ = "Haris Zafeiropoulos"
__credits__ = "Tony J. Lam, Moses Stamboulian, Wontack Han, Yuzhen Ye"
__version__ = "0.1.1"
__maintainer__ = "Haris Zafeiropoulos"
__email__ = "haris.zafeiropoulos@kuleuven.be"
__status__ = "Development"


class PhylomintMGT:

    def __init__(self, config):

        self.dir_path  = config.genres
        self.outdir    = config.seeds
        self.outfile   = "phylomint_scores.tsv"
        self.save_dics =  True
        self.threads   = config.threads
        self.sets_only = config.sets_only
        self.skip_sets = config.skip_sets
        self.prev_conf = config.prev_conf
        self.prev_nonseeds = config.prev_nonseeds

        if self.skip_sets:
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

    def get_sets(self):
        # build initial dictionaries
        SeedSetDic = dict()
        nonSeedSetDic = dict()
        ConfidenceDic = dict()

        # get all XML files in directory
        print("Export seed and non seed sets....")
        sbml_files = [ os.path.join(self.dir_path, f) for f in os.listdir(self.dir_path) if f.endswith('.xml') ]
        # reconstruction_filenames = [ f for f in os.listdir(self.dir_path) if f.endswith('.xml') ]

        num_threads = min(len(sbml_files), self.threads)

        with multiprocessing.Pool(processes=num_threads) as pool:
            with tqdm(total=len(sbml_files), desc="Processing SBML files") as pbar:
                results = []
                for result in pool.imap_unordered(process_sbml, sbml_files):
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

        total_species = self.ConfidenceDic.keys()

        lock = multiprocessing.Lock()
        queue = multiprocessing.Queue()
        processes = []

        # Start a separate process for tracking progress
        progress_process = multiprocessing.Process(target=progress_tracker, args=(queue, len(total_species)))
        progress_process.start()

        # Manually spawn processes with a limited number of concurrent threads
        for i, species in enumerate(total_species):
            if i % self.threads == 0:
                for process in processes:
                    process.join()  # Wait for the batch to finish
                processes = []  # Clear completed processes

            process = multiprocessing.Process(target=self.worker_function, args=(lock, species, queue))
            processes.append(process)
            process.start()

        # Wait for the remaining processes to finish
        for process in processes:
            process.join()

        # Signal the progress tracker to stop
        queue.put(None)
        progress_process.join()


    def scores_and_overlaps_for_a_species(self, species, lock):

        print(species)

        # Get pairs with species as A and species as B
        as_beneficiary, as_donor = generate_fixed_pairwise_comparisons(species, list(self.ConfidenceDic.keys()))

        # Get species seed and non-seed sets
        species_seedset, species_nonseed_set = self.ConfidenceDic[species], self.nonSeedSetDic[species]

        findings = set()

        # Species as A
        for pair in as_beneficiary:
            partner = pair[1]
            SeedSetBConfidence, nonSeedB = self.ConfidenceDic[partner], self.nonSeedSetDic[partner]
            MetabolicCooperationIdxAB, MetabolicCompetitionIdxAB, intersectAB = microbetagPI(species_seedset, SeedSetBConfidence, nonSeedB)
            output = f"{species}\t{partner}\t{MetabolicCompetitionIdxAB}\t{MetabolicCooperationIdxAB}\n"
            findings.add(output)

        # Safely write to the file with a lock
        with lock:
            with open(self.outfile, 'a') as f:
                for r in findings:
                    f.write(r)

        findings = set()
        # Species as B
        for pair in as_donor:
            partner = pair[0]
            SeedSetBConfidence, nonSeedB = self.ConfidenceDic[partner], self.nonSeedSetDic[partner]
            MetabolicCooperationIdxBA, MetabolicCompetitionIdxBA, intersectBA = microbetagPI(SeedSetBConfidence, species_seedset, species_nonseed_set)
            output = f"{partner}\t{species}\t{MetabolicCompetitionIdxBA}\t{MetabolicCooperationIdxBA}\n"
            findings.add(output)

        # Safely write to the file with a lock
        with lock:
            with open(self.outfile, 'a') as f:
                for r in findings:
                    f.write(r)




    def worker_function(self, lock, species, queue):
        """Wrapper function to process a species and signal completion."""
        self.scores_and_overlaps_for_a_species(species, lock)
        with lock:
            queue.put(1)  # Signal that one task is completed

def progress_tracker(queue, total):
    """Progress bar updater."""
    with tqdm(total=total, desc="Processing SBML files") as pbar:
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


def process_sbml(sbml_path, maxcc=2):

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

    return sbml_base, list(SeedSet), nonSeedSet, SeedSetConfidence







