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

db = pd.read_csv(path_root+"/data/dataset/CPdiff_datasetII.csv")
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

labels = ["carbon taxes", "emissions trading", "FiTs", "RPS", "all pricing mechanisms", "all tech. pol."]
colors = ['lightseagreen', 'teal', "indianred", "brown", 'darkkhaki', 'dodgerblue']

## GDP per cap vs Adoption

fig, axs = plt.subplots(2, 2, sharey=True, sharex=True, figsize=(18, 14))

i = 0
j = 0
k = 0

for dummy in ["Tax_dummy", "ETS_dummy", "FiT_dummy", "RPS_dummy"]:
    plots_db = db.loc[db[dummy]==1, :]
    plots_db = plots_db.drop_duplicates(["Country", dummy])
    axs[i, j].scatter(plots_db.Year, plots_db.gdp_pc_ppp, color=colors[k])
    
    axs[i, j].set_title("GDP per cap. PPP x Year - "+labels[k], fontsize=22)
    axs[i, j].set_ylabel("USD", fontsize=18)
    axs[i, j].tick_params(axis="x", labelsize=18)
    axs[i, j].tick_params(axis="y", labelsize=18)
    
    if j == 0:
        j += 1
    elif j > 0:
        i += 1
        j = 0

    k += 1

plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/adoption_scatterII.pdf")


#Adoption and survival - time series


survival = db[["Country", "Year", "Tax_dummy", "ETS_dummy", "FiT_dummy", 
               "RPS_dummy", "Pricing_dummy", "Tech_dummy"]]

survival = survival.groupby(["Year"]).sum()/len(survival.loc[survival.Year==2016,:])

for col in survival.columns:
    survival[col+"_surv"] = 1 - survival[col]


#Cumulative adoption and survival

markers = ['d','3','>','<', '+', '*']

fig = plt.figure(figsize=(20,12))

fig.patch.set_alpha(0.0)
i = 0

ax1 = fig.add_subplot(121)

for dummy in ["Tax_dummy", "ETS_dummy", "FiT_dummy", "RPS_dummy", "Pricing_dummy", "Tech_dummy"]:
    #ax1.plot(cum_count.index, cum_count[dummy], label=labels[i], color=colors[i])
    ax1.plot(survival.index, survival[dummy], label=labels[i], color=colors[i], marker=markers[i])
    i += 1

ax1.set_ylabel("% of countries with policy", fontsize=24)
ax1.tick_params(axis="x", labelsize=20)
ax1.tick_params(axis="y", labelsize=20)
ax1.set_xlim([1990, 2016])
ax1.set_ylim(0, 1.01)

ax2 = fig.add_subplot(122)
i=0

for dummy in ["Tax_dummy_surv", "ETS_dummy_surv", "FiT_dummy_surv", "RPS_dummy_surv", "Pricing_dummy_surv", "Tech_dummy_surv"]:
    ax2.plot(survival.index, survival[dummy], label=labels[i], color=colors[i], marker=markers[i])
    i += 1

ax2.set_ylabel("% of countries without policy", fontsize=24)
ax2.tick_params(axis="x", labelsize=20)
ax2.tick_params(axis="y", labelsize=20)
ax2.set_xlim([1990, 2016])
ax2.set_ylim(0,1.01)

ax1.legend(labels=["carbon tax", "emissions trading", "FiT", "RPS", "all pricing", "all tech. pol."],
           ncol=4, loc="lower center", bbox_to_anchor=(1.1, -0.2),
           fontsize=24)

plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/AdoptionIII.pdf")
plt.close()


# Diffusion regressor: stringency

titles = ["Tax", "ETS", "FiT", "RPS", "Pricing", "Tech. pol."]
figname = ["oecd", "noecd"]


def diffusion_visuals(features, category):
    
    m = 0
    
    for countries in [countries_oecd, countries_noecd]:
        
        fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(18, 18))
        
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
            axs[i, j].set_ylabel("% of total imports", fontsize=20)
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
        
        plt.savefig(path_root+"/Manuscript/Tex_file/Figures/"+category+"_"+figname[m]+"II.pdf")
        m +=1 

diffusion_visuals(["tax_imp", "ets_imp", "fit_imp", "rps_imp", "pricing_imp", "techpol_imp"], "stringency")

# Diffusion: leakage

diffusion_visuals(['tax_leakage_index', 'ets_leakage_index', 'fit_leakage_index',
                   'rps_leakage_index', 'pricing_leakage_index', 'techpol_leakage_index'], 
                  "leakage")

# Diffusion regressors: technology

titles = ["OECD", "non OECD"]          
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

titles = ["Carbon tax", "Emissions trading", "FiT", 
          "RPS"]
figname = ["oecd", "noecd"]

m = 0

for countries in [countries_oecd, countries_noecd]:
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(18, 14))
    
    i = 0
    j = 0
    l = 0

    for feature in ['Tax_impexp', 'ETS_impexp', 'FiT_impexp', 'RPS_impexp']:
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
    
    plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/diffusion_com_"+figname[m]+".pdf")
    m +=1 




