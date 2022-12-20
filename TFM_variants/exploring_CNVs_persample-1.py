import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import zscore
from matplotlib.ticker import ScalarFormatter
plt.style.use('seaborn-colorblind')
sns.set_style('darkgrid')

# Read coverage per base file and add column names to the file
header = ["chr", "start", "end", "feature", "base", "17NR2122", "17NR2123", 
          "17NR2124", "17NR2125", "17NR2126", "17NR2127", "17NR2128", "17NR2129", 
          "17NR2130", "17NR2131", "17NR2132", "17NR2133", "17NR2134", "17NR2135", 
          "17NR2136", "17NR2137", "17NR2138", "17NR2139", "17NR2140", "17NR2141", 
          "17NR2142", "17NR2143", "17NR2144", "17NR2145", "17NR2146", "17NR2147", 
          "17NR2148", "17NR2149", "17NR2150", "17NR2151", "17NR2152", "17NR2153", 
          "17NR2154"]
df = pd.read_csv("./pasted_coverages-chr7-2.tsv", sep="\t", names=header, index_col=False)

dic_total_reads = {"17NR2122":6794137, "17NR2123":3870331, "17NR2124":8761807, "17NR2125":7956425, 
                "17NR2126":4976172, "17NR2127":2399678, "17NR2128":5845975, "17NR2129":7520553, "17NR2130":7544035, 
                "17NR2131":8224731, "17NR2132":6843457, "17NR2133":6604556, "17NR2134":7510350, "17NR2135":7016579, 
                "17NR2136":8359478, "17NR2137":7663546, "17NR2138":9244365, "17NR2139":7487638, "17NR2140":8640616, 
                "17NR2141":7406296, "17NR2142":5602032, "17NR2143":7164379, "17NR2144":7068970, "17NR2145":9307389,
                "17NR2146":7300131, "17NR2147":7013634, "17NR2148":7278219, 
                   "17NR2149":6378808, "17NR2150":8210828, "17NR2151":8707031, 
                   "17NR2152":7102156, "17NR2153":8520552, "17NR2154":7290042}

list_samples = ["17NR2122", "17NR2123", "17NR2124", "17NR2125", 
                "17NR2126", "17NR2127", "17NR2128", "17NR2129", "17NR2130", 
                "17NR2131", "17NR2132", "17NR2133", "17NR2134", "17NR2135", 
                "17NR2136", "17NR2137", "17NR2138", "17NR2139", "17NR2140", 
                "17NR2141", "17NR2142", "17NR2143", "17NR2144", "17NR2145", 
                "17NR2146", "17NR2147", "17NR2148", "17NR2149", "17NR2150", 
                "17NR2151", "17NR2152", "17NR2153", "17NR2154"]

dff = df.copy()

for elem in list_samples:
    dff[elem] = (dff[elem]/dic_total_reads[elem])*10**5

dff["mean"] = dff.loc[:, list_samples].mean(axis=1)
list_samples = list_samples + ["mean"]

c = dff[list_samples].apply(zscore)
for el in list_samples:
    dff[el] = c[el]

for e in list_samples:
    dff[e] = dff[e] - dff["mean"]

## BY ADDING THE START BASE WE GET X AXIS WITH CHROMOSOMIC POSSITIONS
for samplename in list_samples:

    dff["index"] = dff.index

    # Create a plot per gene
    fig = plt.figure(figsize=(40,10))
    ax1 = fig.add_subplot(111)

    # Plot limits
    ax1.set_ylim(-2, 2)
    ax1.set_xlim(dff.loc[0, "start"], dff['index'].max() + dff.loc[0, "start"])
    #ax1.set_xlim(dff.loc[0, "start"] + 350000, dff.loc[0, "start"] + 450000)
    #ax1.set_xlim(57250000, 57260000)
    ax1.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax1.xaxis.offsetText.set_fontsize(22)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    ax1.hlines(y=1, xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[-1], color='r', linestyle='--')
    ax1.hlines(y=-1, xmin=ax1.get_xlim()[0], xmax=ax1.get_xlim()[-1], color='r', linestyle='--')


    for el in list_samples:
        ax1.plot(dff['index'] + dff.loc[0, "start"], dff[el], alpha=0.5, color='gray')
    ax1.plot(dff['index'] + dff.loc[0, "start"], dff[samplename], alpha=0.5, color='blue')
    ax1.plot(dff["index"] + dff.loc[0, "start"], dff["mean"], alpha=0.5, color='r', linestyle='--')


    # Getting positions to plot vertical lines separating exons
#index_positions = [2683191, 2721224]
#names = ['E11-KCNQ1', 'START-KCNQ10T1\nantisense']
        
    # Set limits for the vertical lines, it is a plot on top of a plot
#ax2 = ax1.twiny()
#ax2.set_ylim(ax1.get_ylim())
#ax2.set_xlim(ax1.get_xlim())

    # Paint vertical lines, we use [1:] because the first point is always zero
#ax2.vlines(x=index_positions, ymin=-2, ymax=ax2.get_ylim()[1], linestyle='--', 
#            alpha=0.5, color='black')
#ax2.grid(b=False)
#ax2.set_xticks(index_positions)
#ax2.set_xticklabels(names, rotation=60, minor=False, fontsize=20)

    ax1.set_xlabel('Chromosomal position (bp)', fontsize=24)
    ax1.set_ylabel('Coverage', fontsize=24)


    fig.tight_layout()
#figname = "{}_pbcov.png".format(gene)
    figname = "chr7-2-{}.png".format(samplename)
    fig.savefig(os.path.join("./", figname))
#plt.show()
    plt.close(fig)
