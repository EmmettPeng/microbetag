import os
import unittest

from microbetag.tools import run_seed_complementarity
from microbetag.helpers import MappingPaths


root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
test_data = os.path.join(root_dir, "test_data/test_seed_compl")

models_dir = os.path.join(
    root_dir, "test_data/test_carve/output_files/reconstructions/GENREs/"
)

output_dir = os.path.join(test_data, "output_files")
previous_confidence = os.path.join(output_dir, "confidenceDic.json")
previous_nonseeds = os.path.join(output_dir, "nonSeedSetDic.json")

maps = MappingPaths()


class Config:
    def __init__(self):
        # If False, then microbetag should have built GENREs before running run_seed_complementarity()
        # To use the run_seed_complementarity() on its own, this would have to be always True
        self.users_models = True

        # Directory where GENREs are stored
        self.genres = self.for_reconstructions = models_dir

        # Directory where to save seed complementarity - related files
        self.seeds = output_dir
        self.seed_complements = os.path.join(output_dir, "seed_complements.pckl")
        self.module_seeds = os.path.join(self.seeds, "kegg_module_related_seeds.pckl")
        self.module_nonseeds = os.path.join(
            self.seeds, "kegg_module_related_nonseeds.pckl"
        )

        self.threads = 2
        self.seed_ko_mo = maps.seed_ko_mo
        self.genre_reconstruction_with = "carveme"
        self.metanetx_compounds = maps.metanetx_compounds

        self.skip_sets = False
        self.prev_conf = None
        self.prev_nonseeds = None


class TestSeedComplementarity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = Config()

    def test1WihtoutSets(self):

        run_seed_complementarity(config=self.config)

    def test2WithUsersModelsAndSets(self):

        self.config.skip_sets = True
        self.config.prev_conf = previous_confidence
        self.config.prev_nonseeds = previous_nonseeds

        run_seed_complementarity(config=self.config)


if __name__ == "__main__":

    unittest.main()
