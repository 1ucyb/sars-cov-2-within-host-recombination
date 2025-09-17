import os
import random
import pandas as pd
import numpy as np

##### Empirical null distribution
#
# This script resamples the original dataset to produce an empirical null distribution.
# NOTE: You'll need to set the number of residuals to use this script correctly.

noResiduals = 2624

# Sets random seed(s). Some functions use the native random seed, others use numpy.
np.random.seed(6373627)
random.seed(6373627)

workingList = []
cols1 = ["person_id", "mut_id", "1", "2", "3"]
cols2 = ["person_id", "mut_id", "1", "2"]

# Creates a list of mutations in each individual
for personID in os.listdir("sequences"):
    if personID == ".DS_Store": # Mac moment
        continue
    else:
        workingDF = pd.read_csv("sequences/" + personID + "/" + personID + "_working_subset.csv")
        workingDF.insert(0, "person_id", personID)
        if workingDF.shape[1] == 4:
            workingList.append(pd.DataFrame(workingDF.values, columns=cols2))
        elif workingDF.shape[1] == 5:
            workingList.append(pd.DataFrame(workingDF.values, columns=cols1))
        else:
            print("Error with " + personID)

# Combines into one dataframe
df = pd.concat(workingList, ignore_index = True)

# Generates a list of all person_ids
people = df["person_id"].unique().tolist()

doneList = []
residuals = []
distances = []

i = 0

# Generates the same number of pairs as in the real sample
while i <= int(noResiduals/2):
    # Picks two mutations which have not yet been picked
    sample = random.sample(people, 2)
    sample.sort()
    muts1 = df[df["person_id"].isin([sample[0]])]
    muts2 = df[df["person_id"].isin([sample[1]])]
    mut1 = muts1.sample()
    mut2 = muts2.sample()
    choiceID = [mut1["mut_id"].to_string(), mut2["mut_id"].to_string(), sample[0], sample[1]]
    if choiceID in doneList:
        continue
    else:
        i += 1

        # Carries out the same residual calculation process
        doneList.append(choiceID)

        res1 = abs(float(mut1["1"].iloc[0]) - float(mut2["1"].iloc[0]))
        res2 = abs(float(mut1["2"].iloc[0]) - float(mut2["2"].iloc[0]))

        dist = abs(int(mut1["mut_id"].iloc[0][:-1]) - int(mut2["mut_id"].iloc[0][:-1]))

        residuals.append(res1)
        residuals.append(res2)
        distances.append(dist)
        distances.append(dist)

# Creates and outputs a dataframe of residuals and distances
outputDF = pd.DataFrame(columns = ["residuals", "distance"])
outputDF["residuals"] = residuals
outputDF["distance"] = distances
outputDF.to_csv("empirical-control_residuals_subset.csv", index = False)
