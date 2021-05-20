#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 09:30:41 2020

@author: GD
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

from lifelines.fitters import ParametricRegressionFitter
from lifelines import CoxPHFitter
from lifelines import CoxTimeVaryingFitter
from lifelines import WeibullAFTFitter

from sklearn.preprocessing import MinMaxScaler

## Just making the plots look better
mpl.style.use('ggplot')
mpl.rcParams['figure.figsize'] = (8,6)
mpl.rcParams['font.size'] = 12

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 140)

path_root = "/Users/GD/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion"

db = pd.read_csv(path_root+"/data/dataset/CPdiff_dataset.csv")
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

labels = ["carbon taxes", "emissions trading", "FiTs", "RPS", "All pricing mechanisms"]
colors = ['lightseagreen', 'teal', "indianred", "brown", 'darkkhaki']

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


count = db[["Country", "Year", "Tax_dummy", "ETS_dummy", "FiT_dummy", "RPS_dummy"]]

# Forward fill loop to account for cumulative adoption (and disregard potential repeal)
for ctry in count.Country.unique():
    value_tax = 0
    value_ets = 0

    for yr in count.loc[count.Country==ctry, :].Year.unique():
        if count.loc[(count.Country==ctry) & (count.Year==yr), "Tax_dummy"].item() == 1:
            value_tax = 1
        count.loc[(count.Country==ctry) & (count.Year==yr), "Tax_dummy"] = value_tax
        
        if count.loc[(count.Country==ctry) & (count.Year==yr), "ETS_dummy"].item() == 1:
            value_ets = 1
        count.loc[(count.Country==ctry) & (count.Year==yr), "ETS_dummy"] = value_ets

cum_count = count.groupby(["Year"]).sum()
cum_count["Pricing_dummy"] = cum_count.Tax_dummy + cum_count.ETS_dummy

survival = count.groupby(["Year"]).sum()/len(count.loc[count.Year==2016,:])

for col in survival.columns:
    survival[col+"_surv"] = 1 - survival[col]


#Cumulative adoption and survival

markers = ['d','3','>','<']

fig = plt.figure(figsize=(20,12))

fig.patch.set_alpha(0.0)
i = 0

ax1 = fig.add_subplot(121)

for dummy in ["Tax_dummy", "ETS_dummy", "FiT_dummy", "RPS_dummy"]:
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

for dummy in ["Tax_dummy_surv", "ETS_dummy_surv", "FiT_dummy_surv", "RPS_dummy_surv"]:
    ax2.plot(survival.index, survival[dummy], label=labels[i], color=colors[i], marker=markers[i])
    i += 1

ax2.set_ylabel("% of countries without policy", fontsize=24)
ax2.tick_params(axis="x", labelsize=20)
ax2.tick_params(axis="y", labelsize=20)
ax2.set_xlim([1990, 2016])
ax2.set_ylim(0,1.01)

fig.legend(labels=["carbon tax", "emissions trading", "FiT", "RPS"],
           ncol=4,
           loc="lower center", fontsize=24)


plt.savefig(path_root+"/Manuscript/Tex_file/Figures/description/AdoptionII.pdf")



# Diffusion regressor: stringency

titles = ["Tax", "ETS", "FiT", "RPS"]

figname = ["oecd", "noecd"]

m = 0

for countries in [countries_oecd, countries_noecd]:
    
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(18, 14))
    
    i = 0
    j = 0
    l = 0
    
    for feature in ["Tax_imp", "ETS_imp", "FiT_imp", "RPS_imp"]:
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
    
    plt.savefig(path_root+"/Manuscript/Tex_file/Figures/stringency_"+figname[m]+".pdf")
    m +=1 

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


#Preliminary data analysis
#Countries dropped: 
# - Not in trade matrix: Andorra, Liechtenstein, Monaco, Puerto Rico, Sao Tome and Principe
# - Not in cepii matrix: Romania [??], Serbia [??], South Sudan, Timor-Leste, West Bank and Gaza


#List of countries
ctry_list = db.Country.unique()
ctry_len = len(ctry_list)

db.set_index(["Country", "Year"], inplace=True)


#Notes: ODA-weighted variables are only available for countries which receive ODA

pol_names = {"Tax":"tax", "ETS":"ets", "FiT":"fit", "RPS":"rps"}

#Features included in main models
fea_tech = ['pat_ds_imp', 'imports'] #, , 'pat_ds_exp'
fea_other = ['Duration', 'gdp_pc_ppp', 'ff_perc_tot', 'co2_int', 'v2x_libdem']
fea_pol_main = {}

for policy in pol_names:
    fea_pol_main[policy] = [x for x in db.columns.unique() if ((re.match(pol_names[policy], x.lower()) != None) & (re.search(r'\_EU|\_ODA|\_contig|_wto|_exp', x) == None))]

variable_names = {}

for policy in pol_names:
    variable_names[policy] = {policy+'_dummy':policy,
                      policy+'_avg': "Mean global policy stringency", 
                      policy+'_imp':"Leakage risk",
                      policy+'_exp':"Leakage risk - exp",
                      policy+'_comcol':"Policies of countries with common coloniser", 
                      policy+'_comlang_off':"Policies of countries with common language",
                      policy+'_comrelig':"Policies of countries with common religion", 
                      policy+'_contig':"Policies of contiguous countries", 
                      policy+'_fta_hmr':"Policies of PTA partners", 
                      policy+'_fta_wto':"Policies of RTA partners", 
                      policy+'_impexp':"Policies of trade partners",
                      policy+'_EU':"Common market",
                      policy+'_ODA':"Policies of ODA donors",
                      'imports':'Imports (% GDP)',
                      'ff_perc_tot':"% electricity from fossil fuel",
                      'co2_int':"CO2 intensity",
                      'gdp_pc_ppp':"GDP per capita (PPP)",
                      'pat_ds_exp':"Knowledge stock - exports", 
                      'pat_ds_imp':"Knowledge stock (imports)", 
                      'v2x_libdem':"Polyarchy index", 
                      'v2x_polyarchy':"Liberal democracy index"}

#Summary statistics of relational features
sum_stat_pol = pd.DataFrame()

for policy in pol_names:
    
    temp = db[fea_pol_main[policy]].describe(percentiles=[0.5]).transpose()

    temp.reset_index(inplace=True)
    temp.rename(columns={"index":"Variable"}, inplace=True)
    
    temp["Variable"].replace(to_replace=variable_names[policy], inplace=True)

    temp["Policy"] = policy

    sum_stat_pol = pd.concat([sum_stat_pol, temp])

sum_stat_pol = sum_stat_pol.sort_values(by=["Variable", "Policy"])

#Summary statistics of other and tech features
sum_stat_other =  db[fea_tech+fea_other].describe(percentiles=[0.5]).transpose()
sum_stat_other.reset_index(inplace=True)
sum_stat_other.rename(columns={"index":"Variable"}, inplace=True)
sum_stat_other["Policy"] = ""

#Summary statistics table of all features    
sum_stat = pd.concat([sum_stat_pol, sum_stat_other])
sum_stat = sum_stat.round(decimals=3)
sum_stat["count"] = sum_stat["count"].astype(int) 

sum_stat = sum_stat[["Variable", "Policy", "count", "mean", "std", "min", "50%", "max"]]
sum_stat.columns = ["Variable", "Policy", "Count", "Mean", "Median", "Std", "Min", "Max"]

sum_stat["Variable"].replace(to_replace=variable_names["Tax"], inplace=True)

sum_stat_main = sum_stat.loc[~sum_stat.Variable.isin(["Common market", "Free riding - exp"]), :]

sum_stat_main.reset_index(inplace=True)
sum_stat_main.drop("index", axis=1, inplace=True)

name = ""

for i in range(1, len(sum_stat_main)):
    if name != sum_stat_main["Variable"][i]:
        name = sum_stat_main["Variable"][i]
        
    else:
        sum_stat_main["Variable"][i] = ""

#consider including additional columns before writing if you want more information in the table

sum_stat_main = sum_stat_main.reindex([31, 0, 1, 30, 6, 7, 8, 9, 2, 3, 4, 5, 32, 33,
                       10, 11, 12, 13, 14, 
                       15, 16,
                       17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                       34, 35, 36, 37, 38])

sum_stat_main.to_latex(path_root+"/Manuscript/Tex_file/Tables/sum_stat.tex", index=None)


# Duration analysis
#- Performed for each policy
#- Types of model (parametric): Exponential, Weibull, Cox`

#Model (structure) considerations

# neither Lovely and Popp nor Simmons & Elkins do it (at least not for the diffusion regressors)
#- Check whether to use the CoxTimeVaryingFitter
#- We have multiple observations per country (as Lovely & Popp, 2011) and hence should use robust standard errors
#- Be mindful of trending time-varying covariates - this migh be confounded with the baseline hazard and needs careful econometric treatment
#- There is most likely some collinearity between some of the features - drop them or adjust model or include penalizer
#- Consider introducing a penalty term to shrink coefficients and/or make model sparsers
#- It is also possible to provide a different weight to each observation

#Note: 
#- The lifelines API for fitting hazard models requires to specify 
#   1. set of features; 
#   2. the duration variable and 3. the event variable


#The script must be flexible enough to:
#1. execute estimation for all policies
#2. and all models
#3. accommodate inclusion of different features depending on the policy analysed


# create a scaler object
scaler = MinMaxScaler()

# estimators

class ExponentialAFTFitter(ParametricRegressionFitter):

    # this class property is necessary, and should always be a non-empty list of strings.
    _fitted_parameter_names = ['lambda_']

    def _cumulative_hazard(self, params, t, Xs):
        # params is a dictionary that maps unknown parameters to a numpy vector.
        # Xs is a dictionary that maps unknown parameters to a numpy 2d array
        beta = params['lambda_']
        X = Xs['lambda_']
        lambda_ = np.exp(np.dot(X, beta))
        return t / lambda_



penalty = 0.1

cox = CoxPHFitter(penalizer=penalty)
cox_tv = CoxTimeVaryingFitter(penalizer=penalty) 
#see https://lifelines.readthedocs.io/en/latest/Time%20varying%20survival%20regression.html for how to set up the model (and dataframe!)
#for the 'cox time varying' model, the dataframe needs to have a 'start' and 'stop' column, instead of a 'duration' column.


wei = WeibullAFTFitter()
exp = ExponentialAFTFitter(penalizer=penalty) #multivariate exponential model needs to be defined - see lifelines doc

#AFT models are simply survival models where we allow the survival time to depend on a subject's covariates - see lifelines doc.

model_names = ["exp", "wei", "cox", "cox_tv"]
models = {"cox":cox, "wei":wei, "exp":exp, "cox_tv":cox_tv}
model_labels = {"cox":"CoxHPFitter", "wei":"WeibullFitter", "exp":"ExponentialAFTFitter", "cox_tv":"CoxTimeVaryingFitter"}


#Variables of main specifications
#Should be included in main spec: _
# only policy_imp (not exp)
# not Total_price sector
# only one variable for FTA - pick the most relevant, which is fta_hmr (fta_bb is a multilevel variable, not a dummy and fta_wto is for regional agreements)
# only one variable capturing the political regime


fit = {}
summary = {}

for policy in ["FiT"]:#pol_names:

    #select associated dataframe
    df = pd.read_csv(path_root+"/Data/Data output/Dataset/CPdiff_hz_"+pol_names[policy]+".csv")
    
    #keep only features relevant to the policy for which the model is estimated and exclude some variables not part of the main model
    features = fea_pol_main[policy]+fea_tech+fea_other+["Country"]
    
    df = df[features]
    df = df.dropna()
    df_dur = np.array(df["Duration"].astype(int)) #has to be of int type for exponential model
    df_ctry = np.array(df["Country"])
    df.drop(["Country", "Duration"], axis=1, inplace=True)

    #data scaling - https://towardsdatascience.com/data-normalization-with-pandas-and-scikit-learn-7c1cc6ed6475
    #we don't want to scale the "Duration" variable so we take it out and reintroduce it afterwards
    df_norm = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)
    df_norm["Duration"] = df_dur
    df_norm["Country"] = df_ctry

    fit[policy] = {}
    summary[policy] = {}

    for model in ["exp"]:#model_names:
        #select event variable
        pol_dummy = [x for x in df.columns if re.search("dummy", x) != None]
        
        df_norm[pol_dummy[0]] = df_norm[pol_dummy[0]].astype(int)
        
        #fit model
        if model == "exp":
            df_temp = df_norm.drop("Country", axis=1)
            #df_norm["constant"] = 1.0
            # the below variables maps {dataframe columns, formulas} to parameters
            regressors = {
                'lambda_': df_temp.columns.difference(['Duration', pol_dummy[0]])}
            
            fit[policy][model] = models[model].fit(df_temp, "Duration", #only Weibull - and Gompertz - distributions imply a hazard function with duration dependence 
                                                   pol_dummy[0], 
                                                   regressors=regressors,
                                                   robust=True)
            fig = plt.figure(figsize=(16,12))
            plot = fit[policy][model].plot()
            fig.savefig(path_root+"/Manuscript/Tex_file/Figures/box_plot_coefs_"+policy+"_"+model+".pdf")

                
        if model == "cox":
            df_temp = df_norm.drop("Country", axis=1)
            fit[policy][model] = models[model].fit(df_temp, duration_col="Duration", 
                                                   event_col=pol_dummy[0], 
                                                   robust=True)#, cluster_col="Country")
            
            fig = plt.figure(figsize=(16,12))
            plot = fit[policy][model].plot()
            fig.savefig(path_root+"/Manuscript/Tex_file/Figures/box_plot_coefs_"+policy+"_"+model+".pdf")

        if model == "cox_tv":
            df_temp = df_norm
            
            df_temp["start"] = df_temp["Duration"] - 1
            df_temp.rename(columns={"Duration":"stop"}, inplace=True)
            
            fit[policy][model] = models[model].fit(df_temp, id_col="Country", 
                                                   event_col=pol_dummy[0],
                                                   start_col="start", stop_col="stop")#, cluster_col="Country")
            
            fig = plt.figure(figsize=(16,12))
            plot = fit[policy][model].plot()
            fig.savefig(path_root+"/Manuscript/Tex_file/Figures/box_plot_coefs_"+policy+"_"+model+".pdf")

        elif model == "wei":
            df_temp = df_norm.drop("Country", axis=1)
            fit[policy][model] = models[model].fit(df_temp, duration_col="Duration", 
                                                   event_col=pol_dummy[0], 
                                                   robust=True)#, cluster_col="Country")
            fig = plt.figure(figsize=(16,12))
            plot = fit[policy][model].plot()
            fig.savefig(path_root+"/Manuscript/Tex_file/Figures/box_plot_coefs_"+policy+"_"+model+".pdf")
            
        summary[policy][model] = models[model].summary


#model results to tex table
for policy in pol_names:

    all_models = pd.DataFrame()

    for model in ["wei", "cox", "cox_tv"]:#model_names:
        temp = summary[policy][model].reset_index()
        temp = temp[["covariate", "coef", "se(coef)"]]
        
        temp = temp.round(decimals=3)

        temp["se(coef)"] = temp["se(coef)"].apply(lambda x: "("+str(x)+")")
                
        temp = temp.melt(id_vars="covariate")
        temp.sort_values(by=["covariate", "variable"], inplace=True)
        
        if len(all_models) == 0:
            all_models = temp
        else:
            all_models = all_models.merge(temp, on=["covariate", "variable"]) 
    
    #try and reindex dataframe
                                         
    name = ""

    for i in range(1, len(all_models)):
        if name != all_models.iloc[i, 0]:
            name = all_models.iloc[i, 0]
        
        else:
            all_models.iloc[i, 0] = ""

    all_models.drop("variable", axis=1, inplace=True)
    all_models.columns = ["Variable", "Weibull", "Cox", "Exponential"]
    all_models["Variable"].replace(to_replace=variable_names[policy], inplace=True)
    all_models.to_latex(path_root+"/Manuscript/Tex_file/Tables/results_"+policy+".tex", index=None)
        

#Check the goodness of fit of the model - based on Akaike Information Criterion
#(And select best specification)

#access AIC values with .AIC_ or .AIC_partial_
#see https://lifelines.readthedocs.io/en/latest/fitters/regression/WeibullAFTFitter.html

#Partial effects can also be plotted with command 'plot_partial_effects_on_outcome'
cov_part_plots = {"fta":"Tax_fta_bb", "imp":"Tax_imp"}
        
#for joint effects, use the structure below
i=0
j=0

#fig, axs = plt.subplots(2,2, sharey=True, figsize=(16, 12))
fig = plt.figure(figsize=(16,12))

for policy in pol_names:
    for model in ["cox"]:#model_names:
        axs[i, j] = fit[policy][model].plot_partial_effects_on_outcome(covariates=["gdp_pc_ppp"], 
                                                  #values=np.array([[-0.1, 10000], [0.0, 20000], [0.1, 30000], [0.2, 40000], [0.3, 50000]]), 
                                                  values=np.array([[10000], [20000], [30000], [40000], [50000]]), 
                                                  cmap='coolwarm')
        axs[i, j].set_title(policy)        
        plt.savefig(path_root+"/Manuscript/Tex_file/Figures/partial_gdp_"+policy+"_"+model+".pdf")
        
        
    if j == 0:
        j += 1
    elif j > 0:
        i += 1
        j = 0
    
#For the cox model, check the proportional hazard assumption

#ROBUSTNESS CHECKS    
#It's worth investigating whether:
#- drivers of adoption differ acoss groups of countries (e.g. early vs late tecchnology adopters, low vs high income countries)
#- should we consider a lagged structure?, i.e. reflecting the assumption that information and technology take time to diffuse

rob_check_pol = {"cmarket":"_EU", "oda":"_ODA", "contig":"_contig"}
rob_check = {"ksexp":["pat_ds_exp"], "ksglob":["cum_cap_glob"], "techdep":["cap_ws_ds_exp", "cap_ws_ds_imp"]}

rob_check_list = list(rob_check_pol.keys())+list(rob_check.keys())

fea_rob_check = {}

for policy in pol_names:
    fea_rob_check[policy] = {}
    for rc in rob_check_pol:
        fea_rob_check[policy][rc] = [x for x in db.columns.unique() if (re.search(policy+rob_check_pol[rc], x))]
    for rc in rob_check:
        fea_rob_check[policy][rc] = rob_check[rc]

fit = {}
summary = {}

for policy in pol_names:

    #select associated dataframe
    df = pd.read_csv(path_root+"/Data/Data output/Dataset/CPdiff_hz_"+pol_names[policy]+".csv")
    df = df.loc[~df.Country.isin([ctry_drop]), :]

    summary[policy] = {}
    fit[policy] = {}

    for rc in rob_check_list:

        #keep only features relevant to the policy for which the model is estimated and exclude some variables not part of the main model
        features = [policy+'_dummy']+fea_rob_check[policy][rc]+fea_other#+["Country"]
        
        temp_df = df[features]
        temp_df = temp_df.dropna()
        temp_df_dur = temp_df["Duration"].astype(float) #has to be of int type for exponential model
        temp_df.drop("Duration", axis=1, inplace=True)
    
        #data scaling - https://towardsdatascience.com/data-normalization-with-pandas-and-scikit-learn-7c1cc6ed6475
        #we don't want to scale the "Duration" variable so we take it out and reintroduce it afterwards
        df_norm = pd.DataFrame(scaler.fit_transform(temp_df), columns=temp_df.columns)
        df_norm["Duration"] = np.array(temp_df_dur)
    
        for model in ["wei"]:
            #select event variable
            pol_dummy = [x for x in df.columns if re.search("dummy", x) != None]
            
            df_norm[pol_dummy[0]] = df_norm[pol_dummy[0]].astype(float)
            
            #fit model
    
            if model == "wei":
                fit[policy][model] = models[model].fit(df_norm, duration_col="Duration", 
                                                       event_col=pol_dummy[0], 
                                                       robust=True)#, cluster_col="Country")
                fig = plt.figure(figsize=(16,12))
                plot = fit[policy][model].plot()
                fig.savefig(path_root+"/Manuscript/Tex_file/Figures/box_plot_coefs_"+policy+"_"+rc+".pdf")
                plt.close()
                
            summary[policy][rc] = models[model].summary
            

