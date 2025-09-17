import pandas as pd
import numpy as np
from itertools import combinations
import os

##### Calculating residuals
# The aim here is to calculate the residuals between the differences and a line y=x which is what would
# happen if there was no recombination.

workingColNames = ["mut-id"]

def diffcalc(row):
    # Creates variables allowing me to reference each mutation
    pos = int(row["mut-id"][:-1])
    base = row["mut-id"][-1]
    freqs = []

    # Creates a list of frequencies for that mutation
    for timePoint in timePoints:
        seriesthing = usefulMutations.query("position == @pos & time_point == @timePoint & to_base == @base")["frequency"]
        freqs.append(float(seriesthing.iloc[0]))

    # Creates a list of differences in frequency over time
    for x in range(1, len(timePoints)):
        row.iloc[x] = freqs[x - 1] - freqs[x]

    # Appends this list to the list of dicts
    dictList.append(row.to_dict())

listOfResiduals = []

for personID in os.listdir("sequences"):
    if personID == ".DS_Store":
        continue
    else:
        usefulMutations = pd.read_csv("sequences/"+ personID +"/" + personID + "_useful_mutations_subset.csv")
        print(personID)

        # Setting up working dataframe
        workingDF = pd.DataFrame(columns = workingColNames)

        # Creates a new dataframe which has a row dedicated to each unique mutation
        uniqueCombos = usefulMutations.loc[:, ["position", "to_base"]].drop_duplicates()
        uniqueCombos["position"] = uniqueCombos["position"].astype(np.int64)

        # Generates a unique ID for each mutation based on its position and base
        uniqueCombos["mut_id"] = uniqueCombos["position"].astype(str) + uniqueCombos["to_base"]

        # Creates a list of all the timepoints
        timePoints = list(usefulMutations["time_point"].unique().astype(int))

        timePoints.sort()

        # Creates a new dataframe to be working with, and a list of dicts to append to it
        workingDF["mut-id"] = uniqueCombos["mut_id"]
        dictList = []

        # Creates useful columns
        for i in range(len(timePoints) - 1):
            workingDF["interval_" + str(timePoints[i]) + "-" + str(timePoints[i + 1])] = ""

        workingDF.reset_index(inplace = True, drop = True)

        # Runs function
        workingDF.apply(diffcalc, axis = 1)

        # Converts to dataframe
        workingDF = pd.DataFrame.from_dict(dictList)

        workingDF.to_csv("sequences/"+ personID +"/" + personID + "_working_subset.csv", index = False)

        # Creates a list of all possible pairs of mutations
        pairs = list(combinations(workingDF.index, 2))

        # Creates holding list
        soManyDifferences = []
        soManyDistances = []

        # Calculates the "residuals" of each pair and appends to list
        for pair in pairs:
            row1 = workingDF.iloc[pair[0], 1:]
            row2 = workingDF.iloc[pair[1], 1:]
            diffOfDiffs = (row2 - row1).abs()
            pos1 = workingDF.iloc[pair[0], 0][:-1]
            pos2 = workingDF.iloc[pair[1], 0][:-1]
            distance = abs(int(pos2) - int(pos1))
            listHack = [distance] * len(row1)
            soManyDifferences.append(diffOfDiffs.to_list())
            soManyDistances.append(listHack)

        soManyDifferences = [x
                             for xs in soManyDifferences
                             for x in xs]

        soManyDistances = [x
                           for xs in soManyDistances
                           for x in xs]

        residuals = pd.DataFrame(soManyDifferences, columns = ["residuals"])
        residuals["distance"] = soManyDistances

        listOfResiduals.append(residuals)

        # Saves out the correlation coefficients
        residuals.to_csv("sequences/"+ personID +"/" + personID + "_residuals_subset.csv", index = False)
        print("Saved")

allResiduals = pd.concat(listOfResiduals)
allResiduals.to_csv("all_residuals.csv", index = False)