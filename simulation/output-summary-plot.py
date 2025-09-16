import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os

##### Output Summary Plot
# Creates some plots to examine the output of my SLiM simulations.

# List of recombination rates for future reference
rhos = ["1e-5", "1e-6", "1e-7", "1e-8", "0"]
longRhos = ["1.0e-05", "1.0e-06", "1.0e-07", "1.0e-08", "0"]

bigResidualsList = []

# Iterates over recombination rates
for rho, longRho in zip(rhos, longRhos):

    print(rho)

    # Collects all lists of residuals into one big dataframe
    residualsList = []
    for x in range(1, 51, 4):
        try:
            df = pd.read_csv("output/" + longRho + "/" + str(x) + "/residuals.csv")
            residualsList.append(df)
        except:
            print("There was an issue with run " + str(x) + " at recombination rate " + rho)
    
    residualDF = pd.concat(residualsList)
    residualDF = pd.to_numeric(residualDF["residuals"], downcast = "float").to_frame()

    #Plots and saves histograms describing this dataframe
    print("Creating histogram")
    hist = sns.displot(data = residualDF,
                        x = "residuals").set(title = rho)
    print("Saving histogram")
    hist.savefig("output/" + longRho + "/histogram.png")
    print("Saved")
    
    print("Creating ECDF")
    ecdf = sns.displot(data = residualDF,
                        x = "residuals",
                        kind = "ecdf").set(title = rho)
    print("Saving ECDF")
    ecdf.savefig("output/" + longRho + "/ecdf.png")
    print("Saved")

    # Adds this DF into a big main one along with an identifier for its recombination rate
    print("Saving recombination rate")
    residualDF["rho"] = rho
    print("Appending to main DF")
    bigResidualsList.append(residualDF)

bigResidualsDF = pd.concat(bigResidualsList)
bigResidualsDF = bigResidualsDF.sample(frac = 0.5)

# Plots and saves histograms describing this dataframe
hist = sns.displot(data = bigResidualsDF,
                    x = "residuals",
                    hue = "rho",
                    element = "step",
                    stat = "density",
                    common_norm = False,
                    binwidth = 0.01).set(title = "Histogram of residuals")
hist.savefig("output/histogram.png")

ecdf = sns.displot(data = bigResidualsDF,
                    x = "residuals",
                    kind = "ecdf",
                    hue = "rho").set(title = "ECDF of residuals")
ecdf.savefig("output/ecdf.png")
