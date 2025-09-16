import pandas as pd
pd.options.mode.chained_assignment = None
import sys

#### Parsing the SLiM output
# This script takes the SLiM output to the "useful mutations" stage from the real data analysis.

# Setting parameters
rho = sys.argv[1]
individual = sys.argv[2]
minimumFreq = 0.12
maximumFreq = 1

# Sorts out a formatting issue with the way SLiM handles scientific notation
if rho == "0":
    longRho = rho
else:
    longRho = str(rho)[0] + ".0e-0" + str(rho)[3]

headers = ["out", "tick", "cycle", "tracked", "population", "id", "type", "position", "selection_coef", "dom_coef", "origin_pop", "origin_tick", "prevalence", "to_base"]

# Gets population sizes
log = pd.read_csv("output/" + longRho + "/" + individual + "/sim_log.txt")

# Sets up Pandas dataframe to store the information we will use later
colNames = ["position", "time_point", "to_base", "frequency", "depth", "sample_name"]
df = pd.DataFrame(columns=colNames)

# List for storing dataframes so concatenation is more efficient later
dfList = []

for x in range(1, 4):
    # Converts the main part of the output to a Pandas dataframe
    output = pd.read_csv("output/" + longRho + "/" + individual + "/output_" + str(x) + ".txt", sep=" ", names=headers)

    # Creates a working dataframe for this timepoint
    workingDF = pd.DataFrame(columns = colNames)

    # Works out all columns
    workingDF["position"] = output["position"]
    workingDF["time_point"] = [x] * len(workingDF)
    workingDF["to_base"] = output["to_base"]

    # Calculates the mutation frequencies
    workingDF["frequency"] = output["prevalence"] / int(log.iloc[x-1, 0])

    dfList.append(workingDF)

df = pd.concat(dfList)
df.to_csv("output/" + longRho + "/" + individual + "/output_parsed.csv", index = False)
# This part is the equivalent to mutation-filter.py

# Collects all mutations that are within the frequency requirements at time point 1
midFreqMuts = df.query("@maximumFreq > frequency > @minimumFreq and time_point == 1").drop_duplicates(subset = ["position"])

timePoints = [1, 2, 3]
usefulMuts = []

# This basically collects all the information we have into a usefully organised dataframe
for index, mutation in midFreqMuts.iterrows():
    pos = int(mutation["position"])
    base = mutation["to_base"]

    allTimePoints = df.loc[(df["position"] == pos) & (df["to_base"] == base)]
    allTimePoints.reset_index(drop=True, inplace=True)

    # Aggregates duplicate mutations
    aggregation_functions = {'position': 'first', 'time_point': "first", "to_base": "first", "frequency": "sum", "depth": "first", "sample_name": "first"}
    allTimePoints = allTimePoints.groupby(["position", "time_point", "to_base"], as_index = False).agg(aggregation_functions)

    # Adds zeroes in where necessary
    if len(allTimePoints) < 3:
        for time in timePoints:
            if str(time) not in str(allTimePoints["time_point"].to_list()):
                allTimePoints.loc[len(allTimePoints)] = [pos,
                                                         time,
                                                         base,
                                                         0,
                                                         "NaN",
                                                         "NaN"]

    usefulMuts.append(allTimePoints)

# Saves this information out
dfMuts = pd.concat(usefulMuts).reset_index(drop = True)
dfMuts.to_csv("output/" + longRho + "/" + individual + "/useful_mutations.csv", index = False)
