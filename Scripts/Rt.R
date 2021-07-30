# File:    Rt
# Author:  Janani Mohan
# Date:    2021-06-12
## Reference Code : https://github.com/calldrj/COVID19.Effective.Reproduction.Rate
install.packages("smoother")
install.packages("magrittr")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("lubridate")
install.packages("purrr")
install.packages("plotly")
install.packages("tidyr")
install.packages("zoo")
install.packages("HDInterval")
install.packages("jsonlite")
install.packages("ggplot2")
install.packages("rgdal")
install.packages("rgeos")
install.packages("webshot")
webshot::install_phantomjs()

# Processing functions
# 1. SMOOTH THE DATA IN GAUSSIAN-WINDOW OF ONE-WEEK INTERVAL
library(smoother)   # for smth()
library(magrittr)   # for pipe %>% operator
library(dplyr)      # for mutate(), rename() of columns
library(tidyverse)
library(lubridate)
# Function Smooth()
# * Input: csv dataframe of observations the selected district's date, cases
# * Output: dataframe of observations with the district's cases, smoothed cases, and date
Smooth.Cases <- function(Cases) {
  Cases %>% arrange(date) %>%
    mutate(Cases_Smth=round(smth(Cases, window= 7, tails=TRUE))) %>%
    dplyr::select(date, Cases, Cases_Smth)
}

# 2. VISUALIZE DATA FOR SANITY CHECK
library(plotly)   # for interactive plot ggplotly
Plot.Smth <- function(Smoothed_Cases) {
  plot <- Smoothed_Cases %>% ggplot(aes(x=date, y=Cases)) +
    geom_line(linetype='dotted', color='#429890') + 
    geom_line(aes(y=Cases_Smth), color='#E95D0F') +
    labs(title='Daily Confirmed Cases (Original & Smoothed)', x=NULL, y=NULL) +
    theme(plot.title=element_text(hjust=0.5, color='steelblue'))
} 

# 3. COMPUTE THE EFFECTIVE REPRODUCTION RATE & LOG-LIKELIHOOD 
library(purrr)    # for map() and map
library(tidyr)    # for unnest

RT_MAX <- 10      # the max value of Effective Reproduction Rate Rt
# Generate a set of RT_MAX * 100 + 1 Effective Reproduction Rate value Rt
rt_set <- seq(0, RT_MAX, length=RT_MAX * 1000 + 1)
# Gamma = 1/serial interval
# The serial interval of COVID-19 is defined as the time duration between a primary case-patient (infector) 
# having symptom onset and a secondary case-patient (infectee) having symptom onset. The mean interval was 3.96 days.
# https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
GAMMA = 1/4
# Comp.Likelihood()
# * Input: csv dataframe of observations with the selected district's date, cases, smoothed cases
# * Output: dataframe of observations with the district's cases, smoothed cases, Rt, Rt's log-likelihood
Comp.Log_Likelihood <- function(Acc_Cases) {
  likelihood <- Acc_Cases %>% filter(Cases_Smth > 0) %>%
    # Vectorize rt_set to form Rt column
    mutate(Rt=list(rt_set),
           # Compute lambda starting from the second to the last observation
           Lambda=map(lag(Cases_Smth, 1), ~ .x * exp(GAMMA * (rt_set - 1))),
           # Compute the log likelihood for every observation
           Log_Likelihood=map2(Cases_Smth, Lambda, dpois, log=TRUE)) %>%
    # Remove the first observation
    slice(-1) %>%
    # Remove Lambda column
    select(-Lambda) %>%
    # Flatten the table in columns Rt, Log_Likelihood
    unnest(Log_Likelihood, Rt)
}

# 4. PLOT LIKELIHOOD OF THE EFFECTIVE REPRODUCTION RATE
Plot.Likelihood <- function(likelihoods) {
  likelihoods %>% ggplot(aes(x=Rt, y=Log_Likelihood, group=date)) +
    geom_line(color='#E95D0F', alpha=0.4) +
    labs(title='Daily Likelihood of Rt', subtitle=count) +
    theme(plot.title=element_text(hjust=0.5, color='steelblue'))
}

# 5. COMPUTE THE POSTERIOR OF THE EFFECTIVE REPRODUCTION RATE 
library(zoo)     # for rollapplyr
# Function Comp.Posterior()
# * Input: csv dataframe of observations with the selected state's date, cases, smoothed cases, Rt, Rt's log-likelihood
# * Output: dataframe of observations with the state's cases, smoothed cases, Rt, Rt's posterior
Comp.Posterior <- function(likelihood) {
  likelihood %>% arrange(date) %>%
    group_by(Rt) %>%
    # Compute the posterior for every Rt by a sum of 7-day log-likelihood
    mutate(Posterior=exp(rollapplyr(Log_Likelihood, 7, sum, partial=TRUE))) %>%
    group_by(date) %>%
    # Normalize the posterior 
    mutate(Posterior=Posterior/sum(Posterior, na.rm=TRUE)) %>%
    # Fill missing value of posterior with 0
    mutate(Posterior=ifelse(is.nan(Posterior), 0, Posterior)) %>%
    ungroup() %>%
    # Remove Likelihood column
    dplyr::select(-Log_Likelihood)
}

# 6. PLOT POSTERIOR OF THE EFFECTIVE REPRODUCTION RATE
Plot.Posterior <- function(posteriors) {
  posteriors %>% ggplot(aes(x=Rt, y=Posterior, group=date)) +
    geom_line(color='#E95D0F', alpha=0.2) +
    labs(title='Daily Posterior of Rt', subtitle=count) +
    coord_cartesian(xlim=c(0.2, 5)) +
    theme(plot.title=element_text(hjust=0.5, color='steelblue'))
}

# 7. ESTIMATE THE EFFECTIVE REPRODUCTION RATE
library(HDInterval)
# Function Estimate.Rt()
# * Input: csv dataframe of observations with the selected state's cases, smoothed cases, Rt, Rt's posterior
# * Output: dataframe of observations with the state's Rt, Rt_max, Rt_min
Estimate.Rt <- function(posteriors) {
  posteriors %>% group_by(date) %>%
    summarize(Rts_sampled=list(sample(rt_set, 10000, replace=TRUE, prob=Posterior)),
              Rt_MLL=rt_set[which.max(Posterior)]) %>%
    mutate(Rt_MIN=map_dbl(Rts_sampled, ~ hdi(.x)[1]),
           Rt_MAX=map_dbl(Rts_sampled, ~ hdi(.x)[2])) %>%
    dplyr::select(-Rts_sampled)
}

# 8. PLOT THE THE EFFECTIVE REPRODUCTION RATE'S APPROXIMATION
Plot.Rt <- function(Rt_estimated) {
  plot <- Rt_estimated %>% ggplot(aes(x=date, y=Rt_MLL)) +
    geom_point(color='#429890', alpha=0.5, size=1) +
    geom_line(color='#E95D0F') +
    geom_hline(yintercept=1, linetype='dashed', color='red') +
    geom_ribbon(aes(ymin=Rt_MIN, ymax=Rt_MAX), fill='black', alpha=0.5) +
    labs(title='Estimated Effective Reproduction Rate Rt', x='Time', y='Rt') +
    coord_cartesian(ylim=c(0, 4)) +
    theme(plot.title=element_text(hjust=0.5, color='steelblue'))
}

# 9. ESTIMATE THE EFFECTIVE GROWTH RATE
Estimate.Gt <- function(Cases) {
  growth_rate = Cases %>%
    # first sort by year
    arrange(date) %>%
    mutate(Diff_date = date - lag(date),  # Difference in time (just in case there are gaps)
           Diff_growth = Cases - lag(Cases), # Difference in route between years
           Rate_percent = (Diff_growth / as.double(Diff_date))/Cases * 100) %>%
     # Fill missing value of Rate_percent with 0
     mutate(Rate_percent=ifelse(is.nan(Rate_percent), 0, Rate_percent)) 
}

# 10. PLOT THE THE EFFECTIVE GROWTH RATE'S APPROXIMATION
Plot.Gt <- function(Gt_estimated) {
  plot <- Gt_estimated %>% ggplot(aes(x=date, y=Rate_percent)) +
    geom_line(linetype='dotted', color='#429890') + 
    labs(title='Estimated Growth Rate of Covid 19', x=NULL, y=NULL) +
    theme(plot.title=element_text(hjust=0.5, color='steelblue'))
} 

# 11. PERFORM DATA PROCESSING
Get.Active_Cases <- function(raw_data) {
  active_cases_per_district <- data.frame(raw_data$histories.summary$Thiruvananthapuram$active, 
                                          raw_data$histories.summary$Kollam$active,
                                          raw_data$histories.summary$Pathanamthitta$active,
                                          raw_data$histories.summary$Alappuzha$active,
                                          raw_data$histories.summary$Kottayam$active,
                                          raw_data$histories.summary$Idukki$active,
                                          raw_data$histories.summary$Ernakulam$active,
                                          raw_data$histories.summary$Thrissur$active,
                                          raw_data$histories.summary$Palakkad$active,
                                          raw_data$histories.summary$Malappuram$active,
                                          raw_data$histories.summary$Kozhikode$active,
                                          raw_data$histories.summary$Wayanad$active,
                                          raw_data$histories.summary$Kannur$active,
                                          raw_data$histories.summary$Kasaragod$active,
                                          raw_data$histories.date)
}

Get.Processed.Data <- function() {
  #Get data from Json endpoint
  library(jsonlite)
  raw_data <- as.data.frame(fromJSON("https://keralastats.coronasafe.live/histories.json"))
  active_cases_per_district <- Get.Active_Cases(raw_data)
  colnames(active_cases_per_district) <- c('Thiruvananthapuram',
                                           'Kollam',
                                           'Pathanamthitta',
                                           'Alappuzha',
                                           'Kottayam',
                                           'Idukki',
                                           'Ernakulam',
                                           'Thrissur',
                                           'Palakkad',
                                           'Malappuram',
                                           'Kozhikode',
                                           'Wayanad',
                                           'Kannur',
                                           'Kasaragod',
                                           'date')
  active_cases_per_district$date <- dmy(active_cases_per_district$date)
  # Remove all invalid dates and filter by date range
  start_date <- dmy("01-01-2021")
  end_date <- today()
  active_cases_per_district <- active_cases_per_district %>% 
    filter(date >= start_date & date <= end_date)
  
  #Add lockdown level to data
  active_cases_per_district$lockdown_level <- "level 5"
  active_cases_per_district$lockdown_level <- ifelse(active_cases_per_district$date >=start_date & active_cases_per_district$date <= ymd("2021-05-07") , "Level 4", 
                                                     ifelse(active_cases_per_district$date >= ymd("2021-05-08") & active_cases_per_district$date <= ymd("2021-05-16"), "Level 2", 
                                                            ifelse(active_cases_per_district$date >= ymd("2021-05-17") & active_cases_per_district$date <= ymd("2021-05-21"), "Level 3",
                                                                   ifelse(active_cases_per_district$date >= ymd("2021-05-22") & active_cases_per_district$date <= end_date, "Level 2", "Level 5"))))
  active_cases_per_district <- active_cases_per_district[!duplicated(active_cases_per_district[c('date')]),]
}

library(ggplot2)
library(rgdal)
library(rgeos)
Plot.Kerala.Map <- function() {
  districts_shape = readOGR("./Data", "district")
  class(districts_shape)
  names(districts_shape)
  # print(districts_shape$DISTRICT)
  plot(districts_shape, main = "Administrative Map of Kerala")
}

Plot.Original.Smoothed.Cases.Trend <- function(districts = c('Ernakulam')) {
  active_cases_per_district <- Get.Processed.Data()
  # Plot the original and smoothed cases
  df_cv19 <- list()                         # initialize list of plots for each of districts
  for (i in 1:length(districts)) {
    district <- districts[i]
    df_S <- active_cases_per_district[c(district, "date")] %>% dplyr::select(date, district) %>% 
      rename(Cases=district) %>% 
      Smooth.Cases()
    gplot <- df_S %>% Plot.Smth()
    cl <- rainbow(15)
    plot <- ggplotly(gplot) %>% add_annotations(text=district,
                                                font=list(size=14, color = cl[i]),
                                                xref='paper', yref='paper', x=0, y=0, showarrow=FALSE)
    if (i == 1) {
      plot <- plot %>% add_annotations(text='.... Original',
                                       font=list(size=14, color='#429890'),
                                       xref='paper', yref='paper', x=0.2, y=0, showarrow=FALSE) %>%
        add_annotations(text='--- Smoothed',
                        font=list(size=14, color='#E95D0F'),
                        xref='paper', yref='paper', x=0.4, y=0, showarrow=FALSE)
    }
    df_cv19[[i]] <- plot
  }
  df_cv19 %>% subplot(nrows=length(districts), shareX = TRUE, shareY = TRUE)
}


# Plot the Likelihood of Rt for each district in the list
Plot.Likelihood.Trend <- function(districts = c('Ernakulam')) {
  active_cases_per_district <- Get.Processed.Data()
  df_lh19 <- list()                          # reset list of plots for each of districts
  for (i in 1:length(districts)) {
    district <- districts[i]
    likelihood_by_district <- active_cases_per_district[c(district, "date")] %>% dplyr::select(date, district) %>% 
      rename(Cases=district) %>% 
      Smooth.Cases %>%
      Comp.Log_Likelihood() 
    gplot <- likelihood_by_district %>% Plot.Likelihood()
    
    plot <- ggplotly(gplot) %>% add_annotations(text=district,
                                                font=list(size=14, color='#B51C35'),
                                                xref='paper', yref='paper', x=0, y=1, showarrow=FALSE)
    
    df_lh19[[i]] <- plot
  }
  df_lh19 %>% subplot(nrows=length(districts), shareX=TRUE)
}

# Plot the posteriors of Rt for each district in the list
Plot.Posterior.Trend <- function(districts = c('Ernakulam')) {
  active_cases_per_district <- Get.Processed.Data()
  df_ps19 <- list()                          # reset list of plots for each of districts
  for (i in 1:length(districts)) {
    district <- districts[i]
    posterior_by_district <- active_cases_per_district[c(district, "date")] %>% select(date, district) %>% 
      rename(Cases=district) %>% 
      Smooth.Cases %>%
      Comp.Log_Likelihood() %>%
      Comp.Posterior()
    gplot <- posterior_by_district %>% Plot.Posterior()
    
    plot <- ggplotly(gplot) %>% add_annotations(text=district,
                                                font=list(size=14, color='#B51C35'),
                                                xref='paper', yref='paper', x=0, y=1, showarrow=FALSE)
    
    df_ps19[[i]] <- plot
  }
  df_ps19 %>% subplot(nrows=length(districts), shareX=TRUE)
}

# Compute and plot the max, min, and most-likely values of Rt for each district in the list
Plot.Rt.Trend <- function(districts = c('Ernakulam')) {
  active_cases_per_district <- Get.Processed.Data()
  df_rt19 <- list()                        # reset list of plots for each of district
  df_rtl19 <- list()
  Rt_estimate_list <- list()               # initialize a list of estimated Rt dataframe for each district in the list
  for (i in 1:length(districts)) {
    district <- districts[i]
    df_R <- active_cases_per_district[c(district, "date", "lockdown_level")] %>% dplyr::select(date, district, lockdown_level) %>% 
      rename(Cases=district) %>% 
      Smooth.Cases %>%
      Comp.Log_Likelihood() %>%
      Comp.Posterior() %>%
      Estimate.Rt()
    Rt_estimate_list[[district]] <- df_R
    gplot <- df_R %>% Plot.Rt()
    
    plot <- ggplotly(gplot) %>% add_annotations(text=district,
                                                font=list(size=14, color='#B51C35'),
                                                xref='paper', yref='paper', x=1, y=1, showarrow=FALSE)
    df_rt19[[i]] <- plot
    options(repr.plot.width = 20, repr.plot.height = 8)
    df_R$lockdown_level = active_cases_per_district$lockdown_level[active_cases_per_district$date >= df_R$date[1]] 
    glplot <- df_R  %>%
      ungroup() %>%
      ggplot(aes(x = date, y = Rt_MLL)) +
      geom_col(aes(fill = lockdown_level)) +
      geom_hline(yintercept = 1, linetype = 'dotted') +
      geom_errorbar(aes(ymin = Rt_MIN, ymax = Rt_MAX), width = 0.2) +
      scale_fill_manual(values = c("Level 2" = 'orange', "Level 3" = 'yellow', "Level 4" = 'green')) +
      labs(
        title = 'Most Recent Rt by lockdown',
        x = '', y = ''
      ) +
      theme(axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5))
    options(repr.plot.width = 12, repr.plot.height = 5)
    lplot <- ggplotly(glplot) %>% add_annotations(text=district,
                                        font=list(size=14, color='#B51C35'),
                                        xref='paper', yref='paper', x=1, y=1, showarrow=FALSE)
    df_rtl19[[i]] <- lplot
  }
  df_rt19 %>% subplot(nrows=length(district), shareX=TRUE)
  df_rtl19 %>% subplot(nrows=length(district), shareX=TRUE)
}

# Compute and plot the growth rate
Plot.GrowthRate.Trend <- function(districts = c('Ernakulam')) {
  active_cases_per_district <- Get.Processed.Data()
  df_gt19 <- list()                         # initialize list of plots for each of districts
  for (i in 1:length(districts)) {
    district <- districts[i]
    df_G <- active_cases_per_district[c(district, "date")] %>% select(date, district) %>% 
      rename(Cases=district) %>% 
      Estimate.Gt()
    gplot <- df_G %>% Plot.Gt()
    
    plot <- ggplotly(gplot) %>% add_annotations(text=district,
                                                font=list(size=14, color='#B51C35'),
                                                xref='paper', yref='paper', x=1, y=1, showarrow=FALSE)
    df_gt19[[i]] <- plot
  }
  df_gt19 %>% subplot(nrows=length(districts), shareX=TRUE)
}

# For SEIR Modelling for Ernakulam
# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)

# MODEL INPUTS:

# Vaccine coverage
p <- 0.3

# Total population size
N <- 3430000

# Vector storing the initial number of people in each compartment (at time step 0)
initial_state_values <- c(S = (1-p)*(N-1),   # a proportion 1-p of the total population is susceptible
                          I = 1,             # the epidemic starts with a single infected person
                          R = p*(N-1))       # a proportion p of the total population is vaccinated/immune

# Vector storing the parameters describing the transition rates in units of days^-1
parameters <- c(beta = 0.75,      # the infection rate, which acts on susceptible
                gamma = 0.07)     # the rate of recovery, which acts on those infected

# TIMESTEPS:

# Vector storing the sequence of time steps to solve the model at
times <- seq(from = 0, to = 200, by = 1)   # from 0 to 200 days in daily intervals

# SIR MODEL FUNCTION: 

# The model function takes as input arguments (in the following order): time, state and parameters
sir_model <- function(time, state, parameters) {  
  
  with(as.list(c(state, parameters)), {  # tell R to look for variable names within the state and parameters objects    
    
    # Calculating the total population size N (the sum of the number of people in each compartment)
    N <- S+I+R
    
    # Defining lambda as a function of beta and I:
    lambda <- beta * I/N
    
    # The differential equations
    dS <- -lambda * S               # people move out of (-) the S compartment at a rate lambda (force of infection)
    dI <- lambda * S - gamma * I    # people move into (+) the I compartment from S at a rate lambda, 
    # and move out of (-) the I compartment at a rate gamma (recovery)
    dR <- gamma * I                 # people move into (+) the R compartment from I at a rate gamma
    
    # Return the number of people in the S, I and R compartments at each timestep 
    # (in the same order as the input state variables)
    return(list(c(dS, dI, dR))) 
  })
  
}

# MODEL OUTPUT (solving the differential equations):

# Solving the differential equations using the ode integration algorithm
output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters))

output_long <- melt(as.data.frame(output), id = "time")                  # turn output dataset into long format

# Adding a column for the prevalence proportion to the long-format output
output_long$prevalence <- output_long$value/sum(initial_state_values)

# Plot the prevalence proportion
ggplot(data = output_long,                                               # specify object containing data to plot
       aes(x = time, y = prevalence, colour = variable, group = variable)) +  # assign columns to axes and groups
  geom_line() +                                                          # represent data as lines
  xlab("Time (days)")+                                                   # add label for x axis
  ylab("Prevalence (proportion)") +                                      # add label for y axis
  labs(colour = "Compartment",                                           # add legend title
       title = "Prevalence of infection, susceptibility and recovery over time in Ernakulam")   # add plot title    
