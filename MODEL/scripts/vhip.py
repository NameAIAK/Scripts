import os
import sys
from vhip.predict_interactions import PredictInteractions
from vhip.mlmodel.build import BuildModel

ml_training = '/home/meta/metagenomic/softs/VirusHostInteractionPredictor/data/hostrange/ml_input.csv'

model = BuildModel(ml_training)
VHIP = model.build()

batch = sys.argv[1]

virus_directory_path = batch +'/virus_sequences/'
host_directory_path = batch+'/host_sequences/'

blastn_path = batch +'/StaphStudy_virusvhosts.tsv'
spacer_path = batch +'/StaphStudy_virusvspacers_blastn.tsv'

output_filename = batch +'/predictions.tsv'
CPU_CORES = 6


predictions = PredictInteractions(virus_directory_path, host_directory_path)
predictions.model = VHIP
predictions.add_blastn_files(blastn_path, spacer_path)
predictions.do_setup()
predictions.run_parallel(CPU_CORES)
predictions.predict()
predictions.save_predictions(output_filename)