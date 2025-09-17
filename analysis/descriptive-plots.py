import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
from statsmodels.graphics.gofplots import qqplot_2samples

##### Descriptive plots
#
# This script creates several plots illustrating various aspects of the results.

realResiduals = pd.read_csv("all_residuals.csv")
realResiduals["origin"] = "real"
empiricalNull = pd.read_csv("empirical-control_residuals_subset.csv")
empiricalNull["origin"] = "null"

# Creates a QQ plot illustrating that the two distributions are different
qqplot_2samples(realResiduals["residuals"], empiricalNull["residuals"],
                xlabel = "Real residual quantiles", ylabel = "Empirical null residual quantiles",
                line = "45")

plt.show()

# Joins the two dataframes
df = pd.concat([realResiduals, empiricalNull], ignore_index = True)

# Creates bins for easier display of distance later
bins = range(0, 30000, 1000)
df["binned"] = pd.cut(df["distance"], bins = bins, right = True)

sns.set_theme(style = "whitegrid")

# Plots the distributions of the empirical null and real datasets against each other
sns.histplot(data = df,
             x = "residuals",
             hue = "origin",
             alpha = 0.5,
             binwidth = 0.05,
             legend = False).set(title = "Distribution of residuals", xlabel = "Residuals")

plt.legend(title = "Dataset", labels = ["Empirical null", "Real dataset"])
plt.yscale("log")
plt.show()

# Same but cumulative
sns.histplot(data = df,
             x = "residuals",
             hue = "origin",
             element = "step",
             cumulative = True).set(title = "Cumulative distribution of residuals", xlabel = "Residuals")

plt.legend(title = "Dataset", labels = ["Empirical null", "Real dataset"])
plt.show()

# Shows the relationship of distance to the distribution
g = sns.displot(data = df,
            x = "residuals",
            hue = "binned",
            multiple = "stack",
            binwidth = 0.05,
            palette = "crest",
            linewidth = 0,
                legend = False)
plt.yscale("log")
plt.show()

closeMuts = df.query("distance <= 5000")
farMuts = df.query("distance >= 20000")

# Plots only the nearby mutations
g = sns.displot(data = closeMuts,
            x = "residuals",
            hue = "binned",
            multiple = "stack",
            binwidth = 0.05,
            palette = "crest",
            linewidth = 0,
                legend = False)
plt.yscale("log")
plt.ylim(None,10000)
plt.show()

# Plots only the far apart mutations
g = sns.displot(data = farMuts,
            x = "residuals",
            hue = "binned",
            multiple = "stack",
            binwidth = 0.05,
            palette = "crest",
            linewidth = 0)
plt.yscale("log")
sns.move_legend(g, loc = "center right", labelspacing = 0.1)
plt.ylim(None,10000)
plt.show()
