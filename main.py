import layergraph
import model
import os

# PRACTITIONER DEFINED EXPERIMENTAL INPUTS
# REMEMBER TO CHANGE RUNNAME OR YOU WILL OVERWRITE THE RESULTS!!!
NUM_LAYERS = 10
NUM_PER_LAYER = 10
ARCS_PER_NODE = 3
MAX_EVADERS = 10
SAMPLES = 20
MU = 200
SIGMA = 20
R_0 = 4
DATA_DIR = "dat"
RUNNAME = "run2"
FILENAME = RUNNAME + ".csv"
LOGNAME = RUNNAME + ".log"
FILEPATH = os.path.join(DATA_DIR, FILENAME)
LOGPATH = os.path.join(DATA_DIR, LOGNAME)

ExpBed = layergraph.TestBed(
    NUM_LAYERS, NUM_PER_LAYER, ARCS_PER_NODE, MAX_EVADERS, SAMPLES, MU, SIGMA, R_0)

results = {}
for evaders in range(1, MAX_EVADERS+1):
    current_results = []
    for sample in range(1, SAMPLES+1):
        sample_results = model.M2Model(
            ExpBed.cc[evaders][sample], ExpBed.d, ExpBed.r_0, ExpBed.G.arcs, ExpBed.G.n)
        current_results.append(sample_results)
    results[evaders] = current_results

with open(FILEPATH, "w") as file:
    for evaders in range(1, MAX_EVADERS+1):
        line = ""
        for sample in range(SAMPLES):
            line = line + str(results[evaders][sample][2][0]) + ","
        file.write(line + "\n")

ExpBed.writeBed(LOGPATH)


print(results)
