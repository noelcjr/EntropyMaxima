#!/usr/bin/python

"""
Created on Fri Jun 24 16:49:07 2016
@author: noel
For a description of the program, type:
python flower.py
"""
import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

######## histogram 1 ###############
np.random.seed(0)

n_bins = 10
x = np.random.randn(1000, 3)

fig, axes = plt.subplots(nrows=2, ncols=2)
ax0, ax1, ax2, ax3 = axes.flatten()

colors = ['red', 'tan', 'lime']
ax0.hist(x, n_bins, normed=1, histtype='bar', color=colors, label=colors)
ax0.legend(prop={'size': 10})
ax0.set_title('bars with legend')

ax1.hist(x, n_bins, normed=1, histtype='bar', stacked=True)
ax1.set_title('stacked bar')

ax2.hist(x, n_bins, histtype='step', stacked=True, fill=False)
ax2.set_title('stack step (unfilled)')

# Make a multiple-histogram of data-sets with different length.
x_multi = [np.random.randn(n) for n in [10000, 5000, 2000]]
ax3.hist(x_multi, n_bins, histtype='bar')
ax3.set_title('different sample sizes')

fig.tight_layout()
plt.show()
######## histogram 2 ###############
np.random.seed(0)

# example data
mu = 100    # mean of distribution
sigma = 15  # standard deviation of distribution
x = mu + sigma * np.random.randn(437)

num_bins = 50

fig, ax = plt.subplots()

# the histogram of the data
n, bins, patches = ax.hist(x, num_bins, normed=1)

# add a 'best fit' line
y = mlab.normpdf(bins, mu, sigma)
ax.plot(bins, y, '-')
ax.set_xlabel('Smarts')
ax.set_ylabel('Probability density')
ax.set_title(r'Histogram of IQ: $\mu=100$, $\sigma=15$')

# Tweak spacing to prevent clipping of ylabel
fig.tight_layout()
plt.show()
########################################################################################################################
ro = 0.315061
ep = 0.6
#= (2^(1/6))*TM.ro;
#wallCutOff = 1.122462*ro
regionX = 1.8662
x = [float(i)/100 for i in range(95,150,1)]
fcValx = []
uSum = []
uSum2 = []
for i in x:
    dr = i
    rri = 1/(dr*dr)
    rri3 = rri*rri*rri
    fcValx.append(48 * rri3 * (rri3 - 0.5) *rri *dr)
    uSum.append(4 * rri3 * (rri3 - 1))
    uSum2.append(4 * rri3 * (rri3 - 1) + 1)

DF = pd.DataFrame(columns=['d','E','E2','F'])
DF['d'] = pd.Series(x)
DF['E'] = pd.Series(uSum)
DF['E2'] = pd.Series(uSum2)
DF['F'] = pd.Series(fcValx)
plt.figure()
plt.plot(DF['d'],DF['E'])
plt.plot(DF['d'],DF['E2'],color='red')
plt.axvline(x=1.12)
plt.axhline(y=0.0)
plt.title("Lennard-Jones and Wall Potentials")
plt.xlabel('ro (ro=0.315061 nm)')
plt.ylabel('e (e=0.6364 kj/mol)')
plt.show()
########################################################################################################################
# Circles for GB derivation
circle1 = plt.Circle((3,5),1.0, color='r')
circle2 = plt.Circle((7,5),1.0, color='b')

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()
ax.set_xlim((1, 9))
ax.set_ylim((3, 7))

# major ticks every 20, minor ticks every 5
#major_ticks = np.arange(0, 11, 1)
minor_ticks = np.arange(0, 11, 0.5)
#ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
#ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)

# and a corresponding grid
ax.grid(which='both')
# or if you want differnet settings for the grids:
ax.grid(which='minor', alpha=0.2)
#ax.grid(which='major', alpha=0.8)

ax.plot((7), (5), '.', color='k')
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.plot((3), (5), '.', color='k')
#fig.tight_layout()
plt.rc('grid', linestyle="-", color='k')
plt.grid()
plt.gca().set_aspect('equal', adjustable='box')
#plt.draw()
fig.set_size_inches(10.0, 10.0)
#fig.savefig('/home/noel/plotcircles2.png')
plt.show()
########################################################################################################################
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cwd = os.getcwd()
folder = '/home/noel/Projects/Protein_design/Cas9/aln_pdb/4OO8/GBenergy/'
file_name = '4oo8_abc_swiss_min3_1_self.mmgbsa'
mmgbsa = pd.read_csv(folder+file_name)
mmgbsa.sort(['MMGBSA'],inplace=True)
mmgbsa_dict = {}
for i in mmgbsa.index:
    mmgbsa_dict[mmgbsa.comp1[i]+mmgbsa.comp2[i]] = i  
top = mmgbsa.head(1000000)
#bot = mmgbsa.tail(1000000)
count1 = 0
count2 = 0
init = 0 
final = 100
inter = {}
# Get first 100 components
# For some reasom MMGBSA gives a few components between B and B component interaction.
for i in top.index:
    count2 += 1
    if count1 == final:
        break
    if count1 >= init:
        if mmgbsa.comp1[i] in inter:
            pass
        else:
            # This is because there were no 1-4 interactions and some B-B and C-C interactions where too high
            # and populated the first 100 top energies. The ifs here somehow filter the desired interactions 
            if top.comp1[i][0] == 'A' and top.comp2[i][0] == 'A':
                inter[mmgbsa.comp1[i]] = count1
                count1 += 1
            elif top.comp1[i][0] == 'A' and top.comp2[i][0] == 'B':
                inter[mmgbsa.comp1[i]] = count1
                count1 += 1
            elif top.comp1[i][0] == 'A' and top.comp2[i][0] == 'C':
                inter[mmgbsa.comp1[i]] = count1
                count1 += 1
        if count1 == final:
            break
        else:
            if mmgbsa.comp2[i] in inter:
                pass
            else:
                if top.comp1[i][0] == 'A' and top.comp2[i][0] == 'A':
                    inter[mmgbsa.comp2[i]] = count1
                    count1 += 1
                elif top.comp1[i][0] == 'A' and top.comp2[i][0] == 'B':
                    inter[mmgbsa.comp2[i]] = count1
                    count1 += 1
                elif top.comp1[i][0] == 'A' and top.comp2[i][0] == 'C':
                    inter[mmgbsa.comp2[i]] = count1
                    count1 += 1
#df_top = pd.DataFrame(top_10000,columns=mmgbsa.columns)

mmgbsa_arr = np.zeros((100,100))
polar_arr = np.zeros((100,100))
nonpolar_arr = np.zeros((100,100))
not_found = []
for i in inter:
    for j in inter:
        if i+j in mmgbsa_dict:
            if (i[0] != 'B') and (j[0] != 'B'):
                mmgbsa_arr[inter[i]][inter[j]] = mmgbsa.MMGBSA[mmgbsa_dict[i+j]]
                mmgbsa_arr[inter[j]][inter[i]] = mmgbsa.MMGBSA[mmgbsa_dict[i+j]]
                polar_arr[inter[i]][inter[j]] = mmgbsa.Polar[mmgbsa_dict[i+j]]
                polar_arr[inter[j]][inter[i]] = mmgbsa.Polar[mmgbsa_dict[i+j]]
                nonpolar_arr[inter[i]][inter[j]] = mmgbsa.Nonpolar[mmgbsa_dict[i+j]]
                nonpolar_arr[inter[j]][inter[i]] = mmgbsa.Nonpolar[mmgbsa_dict[i+j]]
        elif j+i in mmgbsa_dict:
            if (i[0] != 'B') and (j[0] != 'B'):
                mmgbsa_arr[inter[i]][inter[j]] = mmgbsa.MMGBSA[mmgbsa_dict[j+i]]
                mmgbsa_arr[inter[j]][inter[i]] = mmgbsa.MMGBSA[mmgbsa_dict[j+i]]
                polar_arr[inter[i]][inter[j]] = mmgbsa.Polar[mmgbsa_dict[j+i]]
                polar_arr[inter[j]][inter[i]] = mmgbsa.Polar[mmgbsa_dict[j+i]]
                nonpolar_arr[inter[i]][inter[j]] = mmgbsa.Nonpolar[mmgbsa_dict[j+i]]
                nonpolar_arr[inter[j]][inter[i]] = mmgbsa.Nonpolar[mmgbsa_dict[j+i]]
        else:
            if i != j:
                not_found.append(i+j)
# http://seaborn.pydata.org/generated/seaborn.heatmap.html
ax = sns.heatmap(mmgbsa_arr)
fig = ax.get_figure()
fig.savefig("/home/noel/Projects/Protein_design/Cas9/aln_pdb/4OO8/GBenergy/mmgbsa.png")
ax2 = sns.heatmap(polar_arr)
ax3 = sns.heatmap(nonpolar_arr)
top.to_csv(folder+'top.csv')
bot.to_csv(folder+'bot.csv')

