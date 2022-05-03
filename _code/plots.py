#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:30:41 2020

@author: GD
"""

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

## Just making the plots look better
mpl.style.use('ggplot')
mpl.rcParams['figure.figsize'] = (8,6)
mpl.rcParams['font.size'] = 12

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 140)

path_root = "/Users/GD/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion"

db = pd.read_csv(path_root+"/data/dataset/CPdiff_dataset_v2.csv")
countries_oecd = ['Sweden', 'United Kingdom', 'Finland','Slovenia', 'Norway', 'Switzerland',
                  'Denmark', 'New Zealand', 'United States', 'Canada']
countries_noecd = ['China', 'South Africa', 'Brazil', 'Russia', 'Philippines',
                   'Indonesia', 'Argentina']

markers_ctry = ['H','x','*','o','s','h','d','3','>','<']
lines_ctry = ['-', '--', '-', '--', '-', '-', '-', '-', '-', '--']
colors_ctry_oecd = ['darksalmon','chocolate','cadetblue','darkorange','darkgreen','navy',
                    'blueviolet','darkcyan','firebrick','darkslategray']
colors_ctry_noecd = ['lightseagreen', 'darkslategray', 'darkkhaki',
                     'olivedrab', 'tan', 'royalblue', 'seagreen']
colors_ctry = [colors_ctry_oecd, colors_ctry_noecd]

labels = ["carbon taxes", "emissions trading", "FiT", "RPS", "all pricing mechanisms", "all tech. pol."]
colors = ['lightseagreen', 'teal', "indianred", "brown", 'darkkhaki', 'dodgerblue']

## GDP per cap vs Adoption

fig, axs = plt.subplots(2, 2, sharey=True, sharex=True, figsize=(18, 14))

i = 0
j = 0
k = 0

for dummy in ["tax", "ets", "fit", "rps"]:
    plots_db = db.loc[db[dummy]==1, :]
    plots_db = plots_db.drop_duplicates(["Country", dummy])
    axs[i, j].scatter(plots_db.Year, plots_db.gdp_pc_ppp, color=colors[k])
    
    axs[i, j].set_title(r"GDP per cap. PPP $\times$ Year - "+labels[k], fontsize=22)
    axs[i, j].set_ylabel("USD", fontsize=18)
    axs[i, j].tick_params(axis="x", labelsize=18)
    axs[i, j].tick_params(axis="y", labelsize=18)
    
    axs[i, j].get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    
    if j == 0:
        j += 1
    elif j > 0:
        i += 1
        j = 0

    k += 1

plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/adoption_scatterII.pdf")


#Adoption and survival - time series

survival = db[["Country", "Year", "tax", "ets", "fit", 
               "rps", "pricing", "techpol"]]

survival = survival.groupby(["Year"]).sum()/len(survival.loc[survival.Year==2016,:])

for col in survival.columns:
    survival[col+"_surv"] = 1 - survival[col]


#Cumulative adoption and survival

markers = ['d','3','>','<', '+', '*']

fig = plt.figure(figsize=(24,14))

fig.patch.set_alpha(0.0)
i = 0

ax1 = fig.add_subplot(121)

for dummy in ["tax", "ets", "fit", "rps", "pricing", "techpol"]:
    #ax1.plot(cum_count.index, cum_count[dummy], label=labels[i], color=colors[i])
    ax1.plot(survival.index, survival[dummy], label=labels[i], color=colors[i], marker=markers[i])
    i += 1

ax1.set_ylabel("% of countries with policy", fontsize=28)
ax1.tick_params(axis="x", labelsize=24)
ax1.tick_params(axis="y", labelsize=24)
ax1.set_xlim([1990, 2016])
ax1.set_ylim(0, 1.01)

ax2 = fig.add_subplot(122)
i=0

for dummy in ["tax_surv", "ets_surv", "fit_surv", "rps_surv", "pricing_surv", "techpol_surv"]:
    ax2.plot(survival.index, survival[dummy], label=labels[i], color=colors[i], marker=markers[i])
    i += 1

ax2.set_ylabel("% of countries without policy", fontsize=28)
ax2.tick_params(axis="x", labelsize=24)
ax2.tick_params(axis="y", labelsize=24)
ax2.set_xlim([1990, 2016])
ax2.set_ylim(0,1.01)

ax1.legend(labels=["carbon tax", "emissions trading", "FiT", "RPS", 
                   "all pricing mechanisms", "all energy policies"],
           ncol=3, loc="lower center", bbox_to_anchor=(1.1, -0.2),
           fontsize=28)

plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/AdoptionIII.pdf")
plt.close()


# Diffusion regressor: leakage index 

titles = ["Carbon tax", "Emissions trading", "FiT", "RPS", 
          "All pricing mechanisms", "All energy policies"]
figname = ["oecd", "noecd"]

def diffusion_visuals(features, category):
    
    y_labels = {"stringency":"% of total imports", "leakage":"Leakage risk index"}
    
    m = 0
    
    for countries in [countries_oecd, countries_noecd]:
        
        fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(18, 18))
        
        i = 0
        j = 0
        l = 0
        
        for feature in features:
            temp = db.loc[:, ["Country", "Year", feature]]
            
            k = 0
            
            for ctry in countries:
    
                temp_ctry = temp.loc[temp.Country==ctry, :]
                axs[i, j].plot(temp_ctry.Year, temp_ctry[feature], color=colors_ctry[m][k], marker=markers_ctry[k], linestyle=lines_ctry[k], label=ctry)
            
                k += 1
                
            axs[i, j].set_title(titles[l], fontsize=24)
            axs[i, j].set_ylabel(y_labels[category], fontsize=20)
            axs[i, j].tick_params(axis="x", labelsize=18)
            axs[i, j].tick_params(axis="y", labelsize=18)
            plt.xlim([1990,2016])
            
            
            if j == 0:
                j += 1
            elif j > 0:
                i += 1
                j = 0
        
            l += 1
            
        fig.legend( ncol=5,
               labels=countries,   # The labels for each line
               loc="lower center",   # Position of legend
               borderaxespad=0.1,    # Small spacing around legend box
               fontsize=20,
               )
        
        plt.savefig(path_root+"/manuscript/Tex_file/figures/"+category+"_"+figname[m]+"II.pdf")
        m +=1 

#diffusion_visuals(["tax_imp", "ets_imp", "fit_imp", "rps_imp", "pricing_imp", "techpol_imp"], "stringency")
diffusion_visuals(["tax_imp", "ets_imp", "fit_imp", "rps_imp"], "stringency")

#diffusion_visuals(['tax_leakage_index', 'ets_leakage_index', 'fit_leakage_index',
#                   'rps_leakage_index', 'pricing_leakage_index', 'techpol_leakage_index'], 
#                  "leakage")

diffusion_visuals(['tax_leakage_index', 'ets_leakage_index', 'fit_leakage_index',
                   'rps_leakage_index'], "leakage")



# Diffusion regressors: technology

titles = ["OECD", "Non-OECD"]          
figname = ["oecd", "noecd"]

fig_xlabels = {'patents_ds_imp':'number of patents (x 100)', 'patents_ds_exp':'number of patents (x 100)',
               'capacity_ws_ds_imp':'Mw', 'capacity_ws_ds_exp':'Mw'}

for feature in ['patents_ds_imp', 'patents_ds_exp', 'capacity_ws_ds_imp', 'capacity_ws_ds_exp']:
    
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(24, 12))
    
    i = 0
    l = 0

    for countries in [countries_oecd, countries_noecd]:
        
        temp = db.loc[:, ["Country", "Year", feature]]
        
        k = 0
        
        for ctry in countries:
            temp_ctry = temp.loc[temp.Country==ctry, :]
            axs[i].plot(temp_ctry.Year, temp_ctry[feature], color=colors_ctry[l][k], marker=markers_ctry[k], linestyle=lines_ctry[k], label=ctry)
        
            k += 1
        
        axs[i].set_ylabel(fig_xlabels[feature], fontsize= 27)
        axs[i].set_xticklabels([1990, 1995, 2000, 2005, 2010, 2015], fontsize= 20)
        axs[i].set_title(titles[l], fontsize=30)
        axs[i].tick_params(axis="y", labelsize=20)
        axs[i].legend(fontsize=26)
        
        if i == 0:
            i += 1
        
        l += 1
        
    plt.xticks([1990, 1995, 2000, 2005, 2010, 2015], fontsize=20)
    plt.xlim([1990,2016])       
    plt.tight_layout()
    plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/diffusion_tech_"+feature+"II.pdf")
     


# Diffusion regressors: communication

m = 0

for countries in [countries_oecd, countries_noecd]:
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(18, 18))
    
    i = 0
    j = 0
    l = 0

    for feature in ['tax_impexp', 'ets_impexp', 'fit_impexp', 'rps_impexp']: #, 'pricing_impexp', 'techpol_impexp'
        temp = db.loc[:, ["Country", "Year", feature]]
        
        k = 0
        
        for ctry in countries:
            temp_ctry = temp.loc[temp.Country==ctry, :]
            axs[i, j].plot(temp_ctry.Year, temp_ctry[feature], color=colors_ctry[m][k], marker=markers_ctry[k], linestyle=lines_ctry[k], label=ctry)
        
            k += 1

        axs[i, j].set_xticklabels([1990, 1995, 2000, 2005, 2010, 2015], fontsize= 18)
        axs[i, j].set_ylabel("% of total bilateral trade", fontsize=20)
        axs[i, j].tick_params(axis="y", labelsize= 18)            
        axs[i, j].set_title(titles[l], fontsize=24)
        plt.xlim([1990,2016])
        
        if j == 0:
            j += 1
        elif j > 0:
            i += 1
            j = 0
    
        l += 1
    
    fig.legend(labels=countries,
               loc='lower center',
               borderaxespad=0.1,    # Small spacing around legend box
               fontsize=20,
               ncol=5)
    
    plt.savefig(path_root+"/manuscript/Tex_file/figures/description/diffusion_com_"+figname[m]+"_v2.pdf")
    m +=1 




