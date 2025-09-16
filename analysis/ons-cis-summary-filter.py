import pandas as pd
import numpy as np

##### Filtering the ONS summary spreadsheet
#
# This script selects individuals from the ONS-CIS summary spreadsheet which meet
# the following criteria:
#
# - At least three samples on different dates
# - All samples within 30 days
# - All samples are the same variant
#
# It should be noted that this doesn't check for multiple stretches of consecutive samples,
# which isn't a problem for my dataset but could be an issue for future readers.

print("Reading in summary spreadsheet...")
df = pd.read_csv("ONS_summary_spreadsheet.csv")

# Counts how many samples there are per person, and removes anybody for whom there's less than three samples
sampleCounts = df.person_id.value_counts()
df = df[df.person_id.isin(sampleCounts.index[sampleCounts > 2])]
df = df.drop_duplicates(subset = ["person_id", "collection_date"], keep = False)

# Adds the number of samples as an extra column on the main dataframe
sampleCounts = sampleCounts[sampleCounts > 2].to_frame() # Excludes everything with less than three samples
df = pd.merge(df, sampleCounts, how="left", on="person_id")
df["count"] = df["count"].astype(int) # Converts the number of samples to int

# Creates a new dataframe with one entry per person for easier iteration
people = df[["person_id", "count"]].drop_duplicates(subset = ["person_id"])
# Creates a new blank column on that dataframe to store date ranges
people["date_range"] = np.nan

# Works out how far apart samples are.
# The author would like to apologise to any actual data scientists reading this for how
# inefficient it is, but she couldn't think of any other way.
print("Iterating... please wait")
for i, row in people.iterrows():

    # Collects all samples for the person in question and creates a list of the days they were taken
    matches = df.query("person_id == @row.person_id")
    listOfDays = matches["Day"].to_list()

    # Sorts the samples and checks if any three items are within 30 days of each other
    # and then saves the date range if so.
    listOfDays.sort(reverse = True)
    x = 2
    while x < len(listOfDays):
        dateRange = listOfDays[x-2] - listOfDays[x]
        if dateRange < 30:
            people.loc[i, "date_range"] = dateRange
            break
        x += 1

    # Checks if all three samples are the same variant
    listOfVariants = matches["paper_lineage"].to_list()
    if len(set(listOfVariants)) != 1:
        people.loc[i, "date_range"] = np.nan

# This way, any samples which are irrelevant to my study remain as NaN and can be filtered out later.
# This merges the two dataframes again
people.drop(columns=["count"], inplace=True) # Drops the count column to avoid duplicating it
df = df.merge(people, how="left", on="person_id")
df = df.dropna(subset = "date_range") # And we can drop the NaN rows here!
print("Saving...")
df.to_csv("filtered_ons_summary.csv", index = False)
print("Saved!")
