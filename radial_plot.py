# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 15:56:17 2019

@author: Katarzyna Luszczak
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
sample_name = 'SAMPLE'                  # type in the sample name
file_name = 'input_data_EDM.csv'        # type in the input file name
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
def load_data(input_file):
    data = pd.read_csv(input_file, delimiter=';')
    
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
        
        age = round(((1 / 1000000) * (1/lambd) * np.log(1 + (lambd * zeta * c * (Ns / A) * rho_d) / (Ni / A))), 2)
        unc = age * ((1 / Ns) + (1 / Ni) + (1 / Nd) + (zeta_se / zeta) ** 2) ** 0.5
    
    return age, unc, comp, comp_type

def draw_radial_plot(age, unc, comp, comp_type):
    
    def xy_coordinates():  
        z = age.apply(np.log)       # using logarthmic transformation (after Vermeesch, 2008)
        sigma_z = unc / age
        central_value = (z / sigma_z ** 2).sum() / (1 /sigma_z ** 2).sum()
        central_age = np.exp(central_value)
        
        x_corr = 1 / sigma_z
        y_corr = (z - central_value) / sigma_z
    
        return x_corr, y_corr, central_value, central_age
    
    def arc_ticks_span():
        age_diff = round((max(age) + 5), -1) - round((min(age) + 5), -1)
        
        if age_diff > 100:
            ticks_span = 50
        elif age_diff > 50:
            ticks_span = 20
        elif age_diff <= 50:
            ticks_span = 10
            
        return ticks_span
    
    def axes_setup(ax, x_corr):
        x_max = round(max(x_corr) + 0.5)       
        ax.set_xlim(0, x_max + 2)
        ax.set_ylim(-3, 3)
        
        ax.set_xlabel('precision (t/$\sigma$)')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_position(('data', -3))
        ax.spines['bottom'].set_bounds(0, x_max + 1) 
        ax.spines['left'].set_bounds(-2, 2)
        ax.set_xticks(np.arange(0, x_max + 2, 1))
        ax.set_yticks([-2, -1, 0, 1, 2])
        
        ax2 = ax.twiny()
        ax2.set_xlabel('% relative error', labelpad=-30)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_position(('data', -2.9))
        ax2.spines['bottom'].set_bounds(0, x_max + 1)
        ax2.spines['left'].set_bounds(-2, 2)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.xaxis.set_label_position('bottom')
        
        x_ticks = np.arange(2, x_max + 2, 1)
        ax2.set_xticks(x_ticks)
        ax2.set_xticklabels((100 / x_ticks).astype(int))
        ax2.set_xlim(ax.get_xlim())
        ax2.tick_params(direction='in', pad=-15)
        
        return x_max
    
    def add_colorbar(plot, comp_type):
        plt.colorbar(plot, orientation='horizontal', aspect=30, pad=0.1, label=comp_type)

    def z_axis_arc(points, radius, central_value, h):
        x_arc = radius / (1 + h**2 * ((points - central_value) ** 2)) ** 0.5
        y_arc = (points - central_value) * x_arc
        return x_arc, y_arc
    
    def arc_ticks_and_labels(tick_ages, x_t, y_t, x_t2, y_t2, ax): 
        for i, value in enumerate(tick_ages):
            x_1 = x_t[i]
            y_1 = y_t[i]
            x_2 = x_t2[i]
            y_2 = y_t2[i]
            plt.plot([x_1, x_2], [y_1, y_2], c='black')
            plt.annotate(value, (x_t2[i] + 0.1, y_t2[i] - 0.1), transform=ax.transAxes)
        
    def draw_arc(plot, ax, x_max, age, ticks_span, central_value):
        R = x_max + 1
        R2 = R + 0.1
        h = (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
        
        tick_ages_0 = np.arange((round((min(age) + 5), -1)), (round(max(age), -1)), ticks_span)
        tick_ages = np.append(tick_ages_0, ([round(min(age), 1), round(max(age), 1)]))
        ticks = np.log(tick_ages)
        
        arc_point_ages = np.arange(min(age),max(age),0.5)
        arc_points = np.log(arc_point_ages)
            
        x_t, y_t = z_axis_arc(ticks, R, central_value, h)
        x_t2, y_t2 = z_axis_arc(ticks, R2, central_value, h)
        x_arc, y_arc = z_axis_arc(arc_points, R, central_value, h)
            
        plt.plot(x_arc, y_arc, c='black', linewidth=1.0)
        arc_ticks_and_labels(tick_ages, x_t, y_t, x_t2, y_t2, ax)
        return R
        
    def draw_central_age_line(R):    
        plt.plot([0,R], [0,0], linestyle=':', c='gray')
        
    def add_text_labels(central_age):    
        plt.text(0, 2.7, sample_name, fontsize=12)
        plt.text(0, 2.4, ('Central age: ' + str(round((central_age), 2)) + ' Ma'), fontsize=11)
                
    def save_plot(sample_name):   
        save_name = sample_name + '_RadialPlot.png'
        plt.savefig(save_name)
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    fig.set_size_inches(8, 10)
    
    x_corr, y_corr, central_value, central_age = xy_coordinates()
    radial_plot = ax.scatter(x_corr, y_corr, c=comp, s=50, cmap='autumn_r', edgecolor='black')
    
    x_max = axes_setup(ax, x_corr)
    add_colorbar(radial_plot, comp_type)
    ticks_span = arc_ticks_span()
    R = draw_arc(radial_plot, ax, x_max, age, ticks_span, central_value)
    draw_central_age_line(R)
    add_text_labels(central_age)
    
    save_plot(sample_name)
    plt.show()

        
def main():
    age, unc, comp, comp_type = load_data(file_name)
    draw_radial_plot(age, unc, comp, comp_type)
       
if __name__ == "__main__":
    main()