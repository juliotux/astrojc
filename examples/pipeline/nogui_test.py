import os
path = os.path.abspath(os.path.dirname('./'))
import sys
sys.path.append(path)
from shutil import rmtree

import numpy as np
from astropy import log

from pipeline.core.utils import empty_value
from pipeline.core.pipeline import Pipeline
import glob

log.setLevel('DEBUG')

ppl = Pipeline('./testing_pipeline.json', nodes_dir='./nodes')
ppl.load_products_file('./test.json')

ppl.run()
