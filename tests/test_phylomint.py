
import os
import unittest

from microbetag.tools import run_seed_complementarity
from microbetag.helpers import MappingPaths


root_dir   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
test_data  = os.path.join(root_dir, "test_data")
output_dir = os.path.join(test_data, "test_phylomint/output_files")
models_dir = os.path.join(test_data, "test_carve/output_files/reconstructions/GENREs/")

maps = MappingPaths()


class Config:

    def __init__(self):
        self.users_models = True
        self.genres       = models_dir
        self.seeds        = output_dir

        self.for_reconstructions =  models_dir  #[os.path.join(models_dir, x) for x in os.listdir(models_dir)]

        self.seed_complements = os.path.join(output_dir, "seed_complements.pckl")
        self.module_seeds     = os.path.join(self.seeds, "kegg_module_related_seeds.pckl")
        self.module_nonseeds  = os.path.join(self.seeds, "kegg_module_related_nonseeds.pckl")

        self.threads = 2
        # New
        self.sets_only = False
        self.skip_sets = True
        self.prev_conf = "../test_data/test_phylomint/output_files/confidenceDic.json"
        self.prev_nonseeds = "../test_data/test_phylomint/output_files/nonSeedSetDic.json"
        self.seed_ko_mo = maps.seed_ko_mo
        self.genre_reconstruction_with = "carveme"
        self.metanetx_compounds = maps.metanetx_compounds


class TestPhylomint(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.config = Config()

    def testWithUsersModels(self):

        run_seed_complementarity(config=self.config)



if __name__ == "__main__":

    unittest.main()

