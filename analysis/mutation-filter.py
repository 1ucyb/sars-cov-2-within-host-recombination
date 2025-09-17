import os
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

##### Calculating and filtering mutation trajectories
#
# The aim of this script is to detect useful mutations and use them to create
# frequency trajectories from BaseFreq files.

# Sets the minimum frequency and minimum read depth to accept as a mutation.
minFreq = 0.12
minDepth = 10
firstTimePointOnly = True # Searches only mutations present in the first time point if true
maxFreq = 1
maf = 10 # The MAF threshold to choose from the artefacts spreadsheet
artefacts = pd.read_csv("artefacts.csv")
samplePrimers = pd.read_csv("sample_list_primers.csv")
# A list of bases which the program can iterate over later.
BASES = ["A", "C", "G", "T"]
GAPBASES = ["A", "C", "G", "T", "gap"]

# Reads in the useful columns from the ONS-CIS summary spreadsheet
ONSSummary = pd.read_csv("filtered_ons_summary.csv",
                         usecols=["person_id", "sequence_name", "Day", "sample_name"],
                         dtype={"person_id": str, "sequence_name": str, "Day": np.int16, "sample_name": str})

# Sets up columns for the information we want on each mutation
colNames = ["position", "time_point", "to_base", "frequency", "depth", "sample_name"]

# Function to apply to dataframe later (it's speedier!)
def append_muts(basefreqs_row):
    # Iterates over each base, checking to see if any available mutations meet the given criteria
    for BASE in BASES:
        depth = basefreqs_row["total"]
        frequency = basefreqs_row[BASE + " freq"]
        if (frequency > minFreq and
                basefreqs_row.iloc[1] != BASE and
                depth >= minDepth):

            # Appends a new row to the summary document with all our new information
            thisDict = {
                "position": basefreqs_row.iloc[0],
                "time_point": timePoint,
                "to_base": BASE,
                "frequency": frequency,
                "depth": depth,
                "sample_name": sampleName
            }
            dictList.append(thisDict)

# Iterates over each individual
print("Searching for mutations and calculating frequencies")
for personID in os.listdir("sequences"):
    if personID == ".DS_Store": # Mac problems :(
        continue
    else:
        # Reads in the relevant summary document rows for the current person
        rows = ONSSummary.query("person_id == @personID")

        # Creates empty list to hold dicts later
        dictList = []

        # Searches for all the BaseFreq files in this person's directory
        for file in os.listdir("sequences/" + personID):
            if file.endswith("BaseFreqs.csv"):

                # Searches for the correct sample ID
                fileStart = file[:11]
                sampleName = rows[rows["sample_name"].str.contains(fileStart)]["sample_name"].iloc[0]

                # Reads in the BaseFreq file
                baseFreqs = pd.read_csv("sequences/" + personID + "/" + file)

                # Adds together all the columns with read counts to give the depth
                baseFreqs["total"] = baseFreqs[[BASE + " count" for BASE in GAPBASES]].sum(axis=1)
                # Calculates the frequency of each variant
                for BASE in BASES:
                    baseFreqs[BASE + " freq"] = baseFreqs[BASE + " count"]/baseFreqs["total"]

                # Works out what the sequence name is so we can find the time point
                mask = rows["sequence_name"].str.contains(file[:11])
                relevantRows = rows.query("@mask")
                timePoint = int(relevantRows["Day"].iloc[0])

                # Applies earlier function to process individual mutations
                baseFreqs.apply(append_muts, axis = 1)

        personSummary = pd.DataFrame.from_dict(dictList)
        # Saves file out
        print("Saving...")
        personSummary.to_csv("sequences/" + personID + "/" + personID + "_summary.csv", index = False)
        print("Saved.")

##### Filtering out mutations
#
# This half of the script filters out mutations even further, and identifies artefacts

print("Filtering mutations...")
for roots, dirs, files in os.walk("sequences"):
    for f in files:
        if f.endswith("_summary.csv"):
            personID = f[:10]

            # Reads in mutations
            mutations = pd.read_csv("sequences/" + personID + "/" + f)

            # Collects all the mutations that at some point are within the frequency requirements
            midFreqMuts = mutations.query("@maxFreq > frequency > @minFreq").drop_duplicates(subset = ["position"])
            midFreqMuts.dropna(axis = 0, inplace = True)

            # Filters out any that don't meet the minimum depth requirements at any point
            midFreqMuts = midFreqMuts.query("depth >= @minDepth")

            # Searches for artefacts
            sampleNames = mutations["sample_name"].unique()
            for sampleName in sampleNames:
                ref = sampleName[:4]
                primer = "NA"
                if ref == "NORT" or ref == "NORW":
                    primer = samplePrimers[samplePrimers.iloc[:, 0] == sampleName].iloc[:, 2]
                    primer = float(primer.iloc[0])
                    artefactList = artefacts.query("centre == @ref & primers == @primer & `MAF (%)` == @maf")["Artefact iSNVs"].to_string(index = False)
                else:
                    artefactList = artefacts.query("centre == @ref & `MAF (%)` == @maf")["Artefact iSNVs"].to_string(index = False)

            # Sets up a dataframe for useful mutations to be stored in
            colNames = ["position", "time_point", "to_base", "frequency", "depth"]
            usefulMuts = []

            # Works out what time points this individual spans
            timePoints = mutations["time_point"].unique()

            # Filters out those that aren't present for the first time point
            if firstTimePointOnly:
                firstTimePoint = timePoints.min()
                midFreqMuts = midFreqMuts.query("time_point == @firstTimePoint")

            # Iterates over all the selected mutations
            for index, mutation in midFreqMuts.iterrows():
                pos = int(mutation["position"])
                base = mutation["to_base"]
                sampleName = mutation["sample_name"]

                # Removes artefacts
                mutID = str(pos) + base
                if mutID in artefactList:
                    continue

                # Collects the information on this mutation for each time point
                allTimePoints = mutations.loc[(mutations["position"] == pos) &
                                              (mutations["to_base"] == base)]
                allTimePoints.reset_index(drop = True, inplace = True)

                # Fills in zeroes where necessary
                if len(allTimePoints) < len(timePoints):
                    for time in timePoints:
                        if str(time) not in str(allTimePoints["time_point"].to_list()):
                            allTimePoints.loc[len(allTimePoints)] = [pos,
                                                                     time,
                                                                     base,
                                                                     0,
                                                                     np.nan,
                                                                     sampleName]

                # Appends to that dataframe
                usefulMuts.append(allTimePoints)

            # Saves out the mutations
            print("Saving...")
            if len(usefulMuts) == 0:
                print("No mutations found for " + personID)
                continue
            else:
                dfMuts = pd.concat(usefulMuts).reset_index(drop = True)
                if firstTimePointOnly:
                    dfMuts.to_csv("sequences/" + f[:10] + "/" + f[:10] + "_useful_mutations_subset.csv", index=False)
                else:
                    dfMuts.to_csv("sequences/" + f[:10] + "/" + f[:10] + "_useful_mutations.csv", index = False)
            print("Saved.")
