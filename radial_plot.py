# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 15:56:17 2019

@author: Katarzyna Luszczak
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# loading the data and setting main variables
#------------------------------------------------------------------------------
sample_name = 'SAMPLE'        # type in the sample name
data = pd.read_csv("input_data_EDM.csv", delimiter=";")

if (data.columns[0] == 'Age' or data.columns[0] == 'age'):
    age = data.iloc[:,0]
    unc = data.iloc[:,1]
    comp = data.iloc[:,2]
    comp_type = data.columns[2]
        
elif (data.columns[0] == 'Ns'):
    Ns = data.iloc[:,0]
    Ni = data.iloc[:,1]
    A = data.iloc[:,2]
    comp = data.iloc[:,3]
    comp_type = data.columns[3]
    
    lambd = 1.55125E-10      # radioactive decay rate of 238U in 1/yr
    c = 0.5
    rho_d = data.iloc[0,4]
    Nd = data.iloc[0,5]
    zeta = data.iloc[0,6]
    zeta_se = data.iloc[0,7]
    
    age = round(((1/1000000) * (1/lambd) * np.log(1 + (lambd * zeta * c * (Ns/A) * rho_d) / (Ni/A))),2)
    unc = age * ((1/Ns) + (1/Ni) + (1/Nd) + (zeta_se/zeta)**2 ) ** 0.5

#------------------------------------------------------------------------------
z = age.apply(np.log)       # using logarthmic transformation (after Vermeesch, 2008)
sigma_z = unc / age
central_value = ((z/(sigma_z**2)).sum()) / ((1/(sigma_z**2)).sum())
central_age = np.exp(central_value)

x_corr = 1/sigma_z
y_corr = (z - central_value)/sigma_z
x_max = round(max(x_corr)+0.5)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Figure
#------------------------------------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
radial_plot = ax.scatter(x_corr, y_corr, c=comp, s=50, cmap='autumn_r', edgecolor='black')
ax.set_xlim(0,(x_max+2))
ax.set_ylim(-3,3)
ax.set_xlabel('precision (t/$\sigma$)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_position(('data',-3))
ax.spines['bottom'].set_bounds(0,x_max+1)

ax.spines['left'].set_bounds(-2,2)
ax.set_xticks(np.arange(0,x_max+2,1))
ax.set_yticks([-2, -1, 0, 1, 2])

ax2 = ax.twiny()
ax2.set_xlabel('% relative error', labelpad=-30)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_position(('data',-2.9))
ax2.spines['bottom'].set_bounds(0,x_max+1)
ax2.spines['left'].set_bounds(-2,2)
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")

x_ticks = np.arange(2,x_max+2,1)
ax2.set_xticks(x_ticks)
ax2.set_xticklabels((100/x_ticks).astype(int))
ax2.set_xlim(ax.get_xlim())
ax2.tick_params(direction='in', pad=-15)

plt.colorbar(radial_plot, orientation='horizontal', aspect=30, pad=0.1, label=comp_type)
fig.set_size_inches(8, 10)

R = x_max + 1
R2 = R + 0.1
h = (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])

age_diff = round((max(age)+5),-1) - round((min(age)+5),-1)
if age_diff > 100:
    ticks_span = 50
elif age_diff > 50:
    ticks_span = 20
elif age_diff <= 50:
    ticks_span = 10

tick_ages_0 = np.arange((round((min(age)+5),-1)), (round(max(age),-1)), ticks_span)
tick_ages = np.append(tick_ages_0, ([round(min(age),1), round(max(age),1)]))
ticks = np.log(tick_ages)

arc_point_ages = np.arange(min(age),max(age),0.5)
arc_points = np.log(arc_point_ages)

x_t = R / (1 + (h ** 2) * ((ticks - central_value) ** 2)) ** 0.5
y_t = (ticks - central_value) * x_t

x_t2 = R2 / (1 + (h ** 2) * ((ticks - central_value) ** 2)) ** 0.5
y_t2 = (ticks - central_value) * x_t2

x_arc = R / (1 + (h ** 2) * ((arc_points - central_value) ** 2)) ** 0.5
y_arc = (arc_points - central_value) * x_arc

arc = plt.plot(x_arc,y_arc, c='black', linewidth=1.0)

plt.plot([0,R],[0,0],linestyle=':',c='gray')

for i, value in enumerate(tick_ages):
    x_1 = x_t[i]
    y_1 = y_t[i]
    x_2 = x_t2[i]
    y_2 = y_t2[i]
    tick_plot = plt.plot([x_1,x_2],[y_1,y_2], c='black')
    plt.annotate(value, (x_t2[i]+0.1, (y_t2[i]-0.1)), transform=ax.transAxes)

plt.text(0, 2.7, sample_name, fontsize=12)
plt.text(0, 2.4, ('Central age: '+ str(round((central_age),2))+' Ma'), fontsize=11)

plt.savefig("RadialPlot.png")

plt.show()