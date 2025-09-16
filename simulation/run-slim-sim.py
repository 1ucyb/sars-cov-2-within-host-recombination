import subprocess
import random
import os

##### Run simulation
# Master program for calling all the constituent parts of a SLiM simulation
# including analysis. Formatted this way so it can be run on the cluster.

individual = 1
#individual = os.environ["SLURM_ARRAY_TASK_ID"]
rhos = ["1e-5", "1e-6", "1e-7", "1e-8", "0"]
random.seed(individual)
randints = []
for i in range(0,len(rhos)):
    randints.append(random.randint(0, 99999))

for rho, randint in zip(rhos, randints):
    # Sorts out a formatting issue with the way SLiM handles scientific notation
    if rho == "0":
        longRho = rho
    else:
        longRho = str(rho)[0] + ".0e-0" + str(rho)[3]
      
    # Creates necessary directories for SLiM
    if not os.path.exists("output/" + longRho + "/" + individual):
        os.makedirs("output/" + longRho + "/" + individual)

    subprocess.run(["run-slim-sim", rho, str(randint), individual])

    for individual in range (1, 51):
        subprocess.run(["python3", "slim-output-parser.py", rho, str(individual)])
        subprocess.run(["python3", "slim-residuals-calc.py", rho, str(individual)])
