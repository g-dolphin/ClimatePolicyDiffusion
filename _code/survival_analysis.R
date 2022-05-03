setwd("/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/data/dataset")

library(devtools)
library(plm)
library(dplyr)
library(stargazer)
library(xtable)
library(survival)
library(eha)
library(rms)

#library(Hmisc)
#library(SurvRegCensCov)
#library(flexsurv)
#library(survminer)

#CoxPH model command specification based on appendix from "An R companion to applied regression"

#For the CoxPH with time-varying covariates model, create start and stop variables
#Introducing the 'FiT_avg' variable in the CoxPH regression leads to convergence problems 


# ------------------------------------- SUMMARY STATISTICS -------------------------------------------- #

df_allobs <- read.csv("CPdiff_dataset_v2.csv")
df_allobs_cov <- subset(df_allobs, select = -c(Country, Year, 
                                               tax_imp, ets_imp, fit_imp, rps_imp, pricing_imp, techpol_imp,
                                               tax_exp,ets_exp, fit_exp, rps_exp, pricing_exp, techpol_exp,
                                               tax_eu, ets_eu, fit_eu, rps_eu, pricing_eu, techpol_eu,
                                               patents_ds_exp, capacity_ws_ds_exp, patents_ds_eu, capacity_ws_ds_eu,
                                               #clim_pol_avg, 
                                               capacity_ds_global, patents_ds_global, 
                                               patents_ds, Duration, exports, co2_tot, co2_pc,
                                               pop_tot, v2x_libdem, co2_int, innovator,
                                               gdp_pc_ppp_sc, pop_tot_sc, ff_perc_tot_sc, co2_int_sc, v2x_polyarchy_sc))

df_allobs_cov <- select(df_allobs_cov, !ends_with(c("_fta_wto", "_contig", "_ODA", "_col45")))

stargazer(df_allobs_cov,
          type = "latex",
          covariate.labels = c("Tax", "ETS", "FiT", "RPS", "All pricing policies", "All techno. policies",
                               "Communication - private - Tax", "Communication - private - ETS", "Communication - private - FiT", "Communication - private - RPS", "Communication - private - All pricing", "Communication - private - All techno.",
                               "Cultural sim. - language - Tax", "Cultural sim. - religion - Tax", "Cultural sim. - colon. - Tax", "Communication - PTA - Tax",
                               "Cultural sim. - language - ETS", "Cultural sim. - religion - ETS", "Cultural sim. - colon. - ETS", "Communication - PTA - ETS",
                               "Cultural sim. - language - FiT", "Cultural sim. - religion - FiT", "Cultural sim. - colon. - FiT", "Communication - PTA - FiT",
                               "Cultural sim. - language - RPS", "Cultural sim. - religion - RPS", "Cultural sim. - colon. - RPS", "Communication - PTA - RPS",
                               "Cultural sim. - language - All pricing", "Cultural sim. - religion - All pricing", "Cultural sim. - colon. - All pricing", "Communication - PTA - All pricing",
                               "Cultural sim. - language - All techno.", "Cultural sim. - religion - All techno.", "Cultural sim. - colon. - All techno.", "Communication - PTA - All techno.",
                               "Knowledge stock - local (imports)", "Technology stock - local (imports)",
                               "Knowledge stock x imports", "Technology stock x imports",
                               "Leakage - Tax", "Leakage - ETS", "Leakage - FiT", "Leakage - RPS", "Leakage - All pricing", "Leakage - All techno.",
                               "Free riding - Tax", "Free riding - ETS", "Free riding - FiT", "Free riding - RPS", "Free riding - All pricing", "Free riding - All techno.",
                               "$\\%$ elec. from fossil fuels$\ ^*$", "$\\%$ elec. ff. x CO2 intensity$\ ^*$",
                               "Import shares", "GDP per capita (PPP, 2017 constant USD)$\ ^*$", "Electoral democracy index$\ ^*$"),
          style = "aer",
          font.size = "footnotesize",
          summary.stat = c("n", "min", "median", "max", "sd"),
          title = "Summary statistics",
          notes.align = "l",
          notes.append = FALSE,
          notes = "\\parbox[t]{\\textwidth}{Variables marked with $\ ^*$ are scaled in the regressions so that a one-unit change represents a 10$\\%$ deviation from the mean. This aids with interpretation and follows Lovely and Popp (2011)}",
          out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/Tables/sumstatR_v2.tex")


# ------------------------------------- MAIN MODEL ESTIMATIONS -------------------------------------------- #

dataframes <- list(fit = "CPdiff_hz_fit_v2.csv", rps = "CPdiff_hz_rps_v2.csv", techno = "CPdiff_hz_techpol_v2.csv",
                   ets = "CPdiff_hz_ets_v2.csv", tax = "CPdiff_hz_tax_v2.csv", pricing = "CPdiff_hz_pricing_v2.csv")

weibull <- list()
exponential <- list()
cox <- list()

for (name in names(dataframes)) {
  print(name)

  df = read.csv(dataframes[[name]])
  
  df$start = df$Duration - 1
  df$stop = df$Duration

  outcome <- names(select(df, "policy_"))

  if (name %in% list("fit", "rps", "techno")) {
  
    #WEIBULL
    #Accessing with [[ or $ is similar. However, it differs for [ in that, indexing with [ will return us a data frame but the other two will reduce it into a vector.
    surv_wei <- Surv(df$Duration, df[[outcome[1]]])
    weibull[[name]] <- phreg(surv_wei ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc,  
                            data = df, 
                            dist='weibull')
    
    #The coefficient estimated by the survreg command are AFT coefficients - they need to be transformed into
    #PH model coefficients. https://www.rdocumentation.org/packages/SurvRegCensCov/versions/1.4/topics/ConvertWeibull
    #weibullph[[i]] <- ConvertWeibull(weibull[[i]], conf.level = 0.95)
    
    #EXPONENTIAL
    exponential[[name]] <- phreg(surv_wei ~ policy_share + policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                           data = df, 
                           dist='weibull',
                           shape = 1)
    
    #COX
    surv_cox <- Surv(df$start, df$stop, df[[outcome[1]]])
    cox[[name]] <- coxph(surv_cox ~ policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc,
                      ties = "breslow",
                      data = df, 
                      cluster = Country)
  }
  
  else if (name %in% list("tax", "ets", "pricing")) {
    
    #WEIBULL
    #Accessing with [[ or $ is similar. However, it differs for [ in that, indexing with [ will return us a data frame but the other two will reduce it into a vector.
    surv_wei <- Surv(df$Duration, df[[outcome[1]]])
    weibull[[name]] <- phreg(surv_wei ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc,  
                          data = df, 
                          dist='weibull')
    
    #The coefficient estimated by the survreg command are AFT coefficients - they need to be transformed into
    #PH model coefficients. https://www.rdocumentation.org/packages/SurvRegCensCov/versions/1.4/topics/ConvertWeibull
    #weibullph[[i]] <- ConvertWeibull(weibull[[i]], conf.level = 0.95)
    
    #EXPONENTIAL
    exponential[[name]] <- phreg(surv_wei ~ policy_share + policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                              data = df, 
                              dist='weibull',
                              shape = 1)
    
    #COX
    surv_cox <- Surv(df$start, df$stop, df[[outcome[1]]])
    cox[[name]] <- coxph(surv_cox ~ policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc,
                      ties = "breslow",
                      data = df, 
                      cluster = Country)
    
  }
  
  
  
  #-------------SURVIVAL PLOTS----------------#
  #https://www.datacamp.com/community/tutorials/survival-analysis-R
  ## Weibull plots
#  weiplot <- ggforest(cox[[i]], main = paste("Survival plot ", "test") )
#  ggsave(paste("wei_main_", name), weiplot, )
  
  ## Cox plot
#  coxplot <- ggforest(cox[[i]], data = df, main = paste("Forest plot ", "test") )
#  ggsave(paste("cox_main_", name), coxplot, path = "/Users/GD/Documents/Activities/Work/Research//Users/GD/Documents/Activities/Work/Research/Research_projects/2020/Climate_policy_diffusion/Manuscript/Tex_file/Figures/")
  
}

# Output table - latex
# https://www.jakeruss.com/cheatsheets/stargazer/
# Note: Weibull with cluster option causes trouble for stargazer function...

stargazer(cox["fit"], cox["rps"], cox["ets"], cox["tax"], weibull["fit"], weibull["rps"], weibull["ets"], weibull["tax"], 
          exponential["fit"], exponential["rps"], exponential["ets"], exponential["tax"],
          type = "latex", out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/tables/results_v2.tex",
          style = "aer", intercept.bottom = TRUE, font.size = "tiny", label = "resultsmain",
          #single.row = TRUE,
          no.space = TRUE,
          column.sep.width = "1pt",
          covariate.labels = c("Free riding", "Leakage", 
                               "Cultural sim. - colon.", "Cultural sim. - com. lang", "Cultural sim - religion", 
                               "Communication - official", "Communication - private", 
                               "Local knowledge stock (KS)",  "Local KS x import share",
                               "$\\%$ elec. from fossil fuel", "GDP per capita", "Electoral dem. index", "Constant"), 
          model.numbers          = FALSE,
          model.names            = TRUE,
          dep.var.labels.include = TRUE,
          dep.var.labels   = c("Cox", "Weibull", "Exponential"),
          column.labels = c("FiT", "RPS", "ETS", "Tax", "FiT", "RPS","ETS", "Tax", "FiT", "RPS","ETS", "Tax"),
          title = "Estimation results") #, float.env = "sidewaystable"

stargazer(cox["techno"], cox["pricing"], weibull["techno"], weibull["pricing"], exponential["techno"], exponential["pricing"],
          type = "latex", out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/tables/results_aggpol_v2.tex",
          style = "aer", intercept.bottom = TRUE, font.size = "tiny", label = "resultsagg",
          #single.row = TRUE,
          no.space = TRUE,
          column.sep.width = "1pt",
          covariate.labels = c("Free riding", "Leakage", 
                               "Cultural sim. - colon.", "Cultural sim. - com. lang", "Cultural sim - religion", 
                               "Communication - official", "Communication - private", 
                               "Local knowledge stock (KS)",  "Local KS x import share",
                               "$\\%$ elec. from fossil fuel", "GDP per capita", "Electoral dem. index", "Constant"), 
          model.numbers          = FALSE,
          model.names            = TRUE,
          dep.var.labels.include = TRUE,
          dep.var.labels   = c("Cox", "Weibull", "Exponential"),
          column.labels = c("FiT or RPS", "Any pricing", "FiT or RPS", "Any pricing", "FiT or RPS", "Any pricing"),
          title = "Policy aggregates") #, float.env = "sidewaystable"

# Additionally, one could report
# AIC values for each model using AIC()
# Number of events

# ----------------------------------------- ROBUSTNESS CHECK REGRESSIONS ---------------------------------- #
#Common market (EU: patents_ds_eu)
#Peer pressure (ODA: policy_ODA)
#Disembodied knowledge (EXPORT WEIGHTED PATENTS: patents_ds_exp)
#Local vs global knowledge? (COUNTRY-SPECIFIC GLOBAL KNOWlEDGE STOCK: patents_ds_global, [patents_ds_global_imp])
#Technology development vs technology deployment? (IMPORT WEIGHTED INSTALLED CAPACITY: capacity_ws_ds_imp, capacity_ws_ds_imp_imp)

#eu <- c("Country", "Duration", "policy_dummy", "policy_avg", "policy_imp", "policy_comcol", "policy_comlang_off", "policy_comrelig", "policy_fta_hmr", "policy_impexp", "patents_ds_eu", "ff_perc_tot", "gdp_pc_ppp", "v2x_polyarchy")
#ks_exp <- c("Country", "Duration", "policy_dummy", "policy_avg", "policy_imp", "policy_comcol", "policy_comlang_off", "policy_comrelig", "policy_fta_hmr", "policy_impexp", "patents_ds_exp", "ff_perc_tot", "gdp_pc_ppp", "v2x_polyarchy")
#ks_glob <- c("Country", "Duration", "policy_dummy", "policy_avg", "policy_imp", "policy_comcol", "policy_comlang_off", "policy_comrelig", "policy_fta_hmr", "policy_impexp", "patents_ds_global", "ff_perc_tot", "gdp_pc_ppp", "v2x_polyarchy")
#ts <- c("Country", "Duration", "policy_dummy", "policy_avg", "policy_imp", "policy_comcol", "policy_comlang_off", "policy_comrelig", "policy_fta_hmr", "policy_impexp", "capacity_ws_ds_imp", "ff_perc_tot", "gdp_pc_ppp", "v2x_polyarchy")
#oda <- c("Country", "Duration", "policy_dummy", "policy_avg", "policy_imp", "policy_comcol", "policy_comlang_off", "policy_comrelig", "policy_fta_hmr", "policy_impexp", "policy_ODA", "patents_ds_imp", "patents_ds_imp_imp", "ff_perc_tot", "gdp_pc_ppp", "v2x_polyarchy")

#rc <- list(eu)#, ks_exp, ks_glob, ts, oda)

#NB: survreg(dist = "weibull") works without problem
#The coefficient estimated by the survreg command are AFT coefficients - they need to be transformed into
#PH model coefficients. https://www.rdocumentation.org/packages/SurvRegCensCov/versions/1.4/topics/ConvertWeibull
#weibullph[[i]] <- ConvertWeibull(weibull[[i]], conf.level = 0.95)

cox_rc <- list()

for (name in names(dataframes)) {
  print(name)
  
  df = read.csv(dataframes[[name]])
  
  df$start = df$Duration - 1
  df$stop = df$Duration

  outcome <- names(select(df, "policy_"))

  surv_wei <- Surv(df$Duration, df[[outcome[1]]])
  surv_cox <- Surv(df$start, df$stop, df[[outcome[1]]])

  
  if (name %in% c("fit", "rps", "techno")) {
    #EU 
    try(cox_rc[[paste(name,"_eu", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds_eu + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                              data = df, 
                              ties = 'breslow',
                              cluster = Country))
    
  
    #Disembodied knowledge 
    try(cox_rc[[paste(name,"_dk", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds_exp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))
  
    #Global knowledge stock 
    try(cox_rc[[paste(name,"_gk", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds_global + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))
  
    #Technology supply 
    try(cox_rc[[paste(name,"_ts", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + capacity_ds_global + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))
  
    
    #ODA - only test for FiTs and RPS because few - if any - countries receiving ODA have implemented carbon tax or ETS
    #try(cox_rc[[i+4]] <- coxph(surv_cox ~  policy_imp + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + policy_ODA + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
    #                         data = df, 
    #                         ties = 'breslow',
    #                         cluster = Country))
    
    #Colonial relationship post 1945
    try(cox_rc[[paste(name,"_col45", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + policy_col45 + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                             data = df, 
                             ties = 'breslow',
                             cluster = Country))
    
    #Domestic knowledge stock
    try(cox_rc[[paste(name,"_doms", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))
    
  }
  
  if (name %in% c("tax", "ets", "pricing")) {
    #EU 
    try(cox_rc[[paste(name,"_eu", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds_eu + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                             data = df, 
                             ties = 'breslow',
                             cluster = Country))
    
    #Disembodied knowledge 
    try(cox_rc[[paste(name,"_dk", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds_exp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))
    
    #Global knowledge stock 
    try(cox_rc[[paste(name,"_gk", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds_global + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))
    
    #Technology supply 
    try(cox_rc[[paste(name,"_ts", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + capacity_ds_global + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                               data = df, 
                               ties = 'breslow',
                               cluster = Country))

    #Colonial relationship post 1945
#    try(cox_rc[[paste(name,"_col45")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + policy_col45 + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
#                                                data = df, 
#                                                ties = 'breslow',
#                                                cluster = Country))
    
    #Domestic knowledge stock
    try(cox_rc[[paste(name,"_doms", sep="")]] <- coxph(surv_cox ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + patents_ds + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc, 
                                               data = df, 
                                               ties = 'breslow',
                                               cluster = Country))
    
  }

  }

# Output table - latex
stargazer(cox_rc["fit_eu"], cox_rc["fit_dk"], cox_rc["fit_gk"], cox_rc["fit_ts"], cox_rc["fit_col45"], cox_rc["fit_doms"], 
          cox_rc["rps_eu"], cox_rc["rps_dk"], cox_rc["rps_gk"], cox_rc["rps_ts"], cox_rc["rps_col45"], cox_rc["rps_doms"],
          type = "latex", out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/tables/results_rcI_v2.tex",
          style = "aer", intercept.bottom = TRUE, font.size = "tiny", label = "resultsrcI",
          #single.row = TRUE,
          no.space = TRUE,
          column.sep.width = "1pt",
          covariate.labels = c("Leakage", 
                               "Cultural sim. - colon.", "Cultural sim. - language", "Cultural sim. - relig.", 
                               "Communication - official", "Communication - private",  
                               "Local knowledge stock (KS)", "Local KS x import share",
                               "Knowledge stock EU", "Disembodied knowledge",
                               "Global foreign KS", "Global renewables cap.", "Peer pressure - colonizer", "Domestic KS",
                               "$\\%$ elec. from fossil fuel", "GDP per capita", "Electoral dem. index", "Constant"), 
          model.numbers          = FALSE,
          model.names            = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("SM", "DK", "Glob.", "Tech. Dep.", "Colonizer", "Dom.", "SM", "DK", "Glob.", "Tech. Dep.", "Colonizer", "Dom."),
          title = "Robustness checks - FiT $\\&$ RPS") #, float.env = "sidewaystable"

stargazer(cox_rc["ets_eu"], cox_rc["ets_dk"], cox_rc["ets_gk"], cox_rc["ets_ts"], cox_rc["ets_doms"], 
          cox_rc["tax_eu"], cox_rc["tax_dk"], cox_rc["tax_gk"], cox_rc["tax_ts"], cox_rc["tax_doms"],
          type = "latex", out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/tables/results_rcII_v2.tex",
          style = "aer", intercept.bottom = TRUE, font.size = "tiny", label = "resultsrcII",
          #single.row = TRUE,
          no.space = TRUE,
          column.sep.width = "1pt",
          covariate.labels = c("Leakage", 
                               "Cultural sim. - colon.", "Cultural sim. - language", "Cultural sim. - relig.", 
                               "Communication - official", "Communication - private",  
                               "Local knowledge stock (KS)", "Local KS x import share", "Patent stock EU", "Disembodied KS",
                               "Global foreign KS", "Global renewables cap.", "Domestic KS", #"Peer pressure - colonizer",
                               "$\\%$ elec. from fossil fuel", "GDP per capita", "Electoral dem. index", "Constant"), 
          model.numbers          = FALSE,
          model.names            = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("SM", "DK", "Glob.", "Tech. Dep.", "Dom.", "SM", "DK", "Glob.", "Tech. Dep.", "Dom."),
          title = "Robustness checks - ETS $\\&$ Tax") #, float.env = "sidewaystable"

stargazer(cox_rc["techno_eu"], cox_rc["techno_dk"], cox_rc["techno_gk"], cox_rc["techno_ts"], cox_rc["techno_doms"], 
          cox_rc["pricing_eu"], cox_rc["pricing_dk"], cox_rc["pricing_gk"], cox_rc["pricing_ts"], cox_rc["pricing_doms"],
          type = "latex", out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/tables/results_rcIII_v2.tex",
          style = "aer", intercept.bottom = TRUE, font.size = "tiny", label = "resultsrcII",
          #single.row = TRUE,
          no.space = TRUE,
          column.sep.width = "1pt",
          covariate.labels = c("Leakage", 
                               "Cultural sim. - colon.", "Cultural sim. - language", "Cultural sim. - relig.", 
                               "Communication - official", "Communication - private",  
                               "Local knowledge stock (KS)", "Local KS x import share", "Patent stock EU", "Disembodied KS",
                               "Global foreign KS", "Global renewables cap.", "Domestic KS", #"Peer pressure - colonizer",
                               "$\\%$ elec. from fossil fuel", "GDP per capita", "Electoral dem. index", "Constant"), 
          model.numbers          = FALSE,
          model.names            = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = c("SM", "DK", "Glob.", "Tech. Dep.", "Dom.", "SM", "DK", "Glob.", "Tech. Dep.", "Dom."),
          title = "Robustness checks - All energy policies $\\&$ All pricing policies") #, float.env = "sidewaystable"

# Non-innovating countries

cox_noninn <- list()

for (name in list("fit", "rps", "ets", "tax")) {
  print(name)
  df = read.csv(dataframes[[name]])
  
  df <- subset(df, innovator == 0)
  
  df$start = df$Duration - 1
  df$stop = df$Duration
  
  outcome <- names(select(df, "policy_"))
  
  #Cox
  surv_wei <- Surv(df$Duration, df[[outcome[1]]])

  # Implementing different model specications for different policies
  if (name %in% list("fit", "rps", "ets")) {
    cox_noninn[[paste(name,"_noninn", sep="")]] <- coxph(surv_wei ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc,  
                            data = df, 
                            tie = 'breslow',
                            cluster = Country)#,
  }
  
  else {
    cox_noninn[[paste(name,"_noninn", sep="")]] <- coxph(surv_wei ~  policy_leakage_index + policy_comcol + policy_comlang_off + policy_comrelig + policy_fta_hmr + policy_impexp + patents_ds_imp + patents_ds_imp_imp + ff_perc_tot_sc + gdp_pc_ppp_sc + v2x_polyarchy_sc,  
                                                         data = df, 
                                                         tie = 'breslow',
                                                         cluster = Country)#,
  }

}

stargazer(cox_noninn["fit_noninn"], cox_noninn["rps_noninn"], cox_noninn["ets_noninn"], cox_noninn["tax_noninn"], 
          type = "latex", out = "/Users/gd/OneDrive - rff/Documents/Research/projects/climate_policy_diffusion/manuscript/Tex_file/tables/results_noninn.tex",
          style = "aer", intercept.bottom = TRUE, font.size = "tiny", label = "resultsnoninn",
          #single.row = TRUE,
          no.space = TRUE,
          column.sep.width = "1pt",
          covariate.labels = c("Leakage", 
                               "Cultural sim. - colon.", "Cultural sim. - com. lang", "Cultural sim - religion", 
                               "Communication - official", "Communication - private", 
                               "Local knowledge stock (KS)",  "Local KS x import share",
                               "$\\%$ elec. from fossil fuel", "GDP per capita", "Electoral dem. index", "Constant"), 
          #dep.var.labels.include = FALSE,
          model.numbers          = FALSE,
          model.names            = TRUE,
          dep.var.labels.include = FALSE,
          #dep.var.labels   = c("Cox", "Exponential", "Weibull"),
          column.labels = c("FiT", "RPS", "ETS", "Tax", "FiT", "RPS","ETS", "Tax", "FiT", "RPS","ETS", "Tax"),
          title = "Non-innovating countries") #, float.env = "sidewaystable"

