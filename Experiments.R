library(tidyverse)
library(MagmaClustR)

##### Data import #####

## Full data can be requested here: https://gustodatavault.sg/about/request-for-data
## An illustrative sub-sample is used to demonstrate reproducibility 
db = read_csv('Data/BMI_birth_to_10years.csv')

##### Utility functions for prediction, evaluation, plotting #####

MSE = function(obs, pred)
{
  input = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mix_pred = pred %>%
    filter(Input %in% input) %>%
    pull(Mean)
  
  (value - mix_pred)^2 %>%
    mean() %>%
    return()
}

WCIC = function(obs, pred, level)
{
  t = obs %>% pull(Input)
  value = obs %>% pull(Output)
  
  mean = pred %>% arrange(Input) %>% pull(Mean)
  sd = pred %>% arrange(Input) %>% pull(Var) %>% sqrt
  
  CI_inf = mean - qnorm(1 - level/2) * sd
  CI_sup = mean + qnorm(1 - level/2) * sd
  
  100 * ((CI_inf < value) & (value < CI_sup)) %>%
    mean %>%
    return()
}

eval = function(db, mod, db_train = NULL, hyperpost = NULL, remove_random = F,
                ratio_missing = 0.5, start_input_test = 6, name = 'MagmaClust')
{
  if(hyperpost %>% is.null()){
    hyperpost = hyperposterior_clust(
      trained_model = mod,
      data = db_train,
      grid_inputs = unique(db$Input) %>% sort())
  }
  
    floop_j = function(j){
      cat('ID n°',j, '\n \n')

      ID_pred = j
      db_j = db %>% 
        filter(ID %in% ID_pred)

      if(remove_random) {
        size = max(length(db_j$Input) * (1-ratio_missing), 1)
        input_test = sample(db_j$Input, size = size)
      } else {
        input_test = db_j %>% filter(Input < start_input_test) %>% pull(Input)
      }
      
      db_test = db_j %>% 
        filter(Input %in% input_test) %>% 
        dplyr::select(- ID)
      
      test_point = db_j %>% 
        filter(!(Input %in% input_test))
      
      input_pred = test_point$Input

      pred = pred_magmaclust(db_test,
                     trained_model = mod,
                     grid_inputs = input_pred,
                     hyperpost = hyperpost,
                     plot = FALSE)
      
      mixt_pred = pred$mixture_pred %>% 
        mutate(Input = input_pred)
      
      eval = tibble('Method' = name,
                    'ID' = j, 
                    'MSE' = test_point %>% MSE(mixt_pred),
                    'WCIC' = test_point %>% WCIC(mixt_pred, 0.05))
    
      if(remove_random)
      {
        eval = eval %>% mutate('Ratio_missing' = ratio_missing)
      }
      
      eval %>%
        return()
    } 
    db$ID %>%
      unique %>% 
      lapply(floop_j) %>%
      bind_rows %>% 
      return()
}

loop_pred = function(db, mod, hyperpost = NULL, remove_random = F,
                     ratio_missing = 0.5, start_input_test = 6, 
                     grid_input = NULL, name = NULL)
{
  if(hyperpost %>% is.null()){
    hyperposterior_clust(
      trained_model = mod,
      data = db,
      grid_inputs = unique(db$Input) %>%
        union(unique(mod$ini_args$data$Input)) %>%
                             sort())
  }
  
  floop_j = function(j){
    cat('ID n°',j, '\n \n')

    ID_pred = j
    db_j = db %>% 
      filter(ID %in% ID_pred)
    
    if(remove_random) {
      size = max(length(db_j$Input) * (1-ratio_missing), 1)
      input_test = sample(db_j$Input, size = size)
    } else {
      input_test = db_j %>% filter(Input < start_input_test) %>% pull(Input)
    }

    db_test = db_j %>%
      filter(Input %in% input_test)

    test_point = db_j %>%
      filter(!(Input %in% input_test))
    
    input_pred = test_point$Input
    
    if(!(grid_input %>% is.null)){
      input_pred = input_pred %>% 
        c(grid_input) %>% 
        sort()
    }
    
    if(name == 'Splines'){
      if(nrow(db_test) > 4){
        spline_fit = smooth.spline(db_test$Input, db_test$Output)
        pred = predict(spline_fit, input_pred)
      } else {
        pred = list(y = NA)
      }
      
      mixt_pred = tibble('ID' = ID_pred, 'Input' = input_pred,
                         'Mean' = pred$y, 'Var' = 0)
      
    } else{
      
      pred = pred_magmaclust(db_test,
                             trained_model = mod,
                             grid_inputs = input_pred,
                             hyperpost = hyperpost,
                             plot = FALSE)

      mixt_pred = pred$mixture_pred
      
      input_pred = mixt_pred$Input
    }
      res = mixt_pred %>% 
      arrange(Input) %>%
      mutate(Input = round(Input, 3)) %>% 
        left_join(test_point %>% mutate(Input = round(Input, 3)),
                  by = c('ID', 'Input'))
      
      if(!is.null(name)){
        res = res %>% mutate('Method' = name)
      }
      
      return(res)
  } 
  db$ID %>%
    unique %>% 
    lapply(floop_j) %>%
    bind_rows %>% 
    return()
}

summarise_res = function(res, by = Method, digits = 4){
  res %>% 
    dplyr::select(- ID) %>%
    group_by({{ by }}) %>% 
    reframe(across(all_of(c('MSE', 'WCIC')),
                     list('Mean' = mean, 'SD' = sd), na.rm = TRUE)) %>% 
    mutate(across(MSE_Mean:WCIC_SD, \(x) round(x, digits))) %>% 
    mutate('MSE' = paste0(MSE_Mean, ' (', MSE_SD, ')'),
           'WCIC' =  paste0(WCIC_Mean, ' (', WCIC_SD, ')')) %>%
    dplyr::select(c({{ by }}, MSE, WCIC)) %>% 
    return()
}

plot_error = function(pred){
  db_error = pred %>%
    mutate(Error = Mean - Output) %>% 
    mutate(CI_inf = - qnorm(0.975) * sqrt(Var)) %>% 
    mutate(CI_sup = qnorm(0.975) * sqrt(Var)) %>% 
    arrange(Var) %>% 
    rowid_to_column(var = 'Index') %>% 
    mutate(Index = as.double(Index))
  
  gg_error = ggplot(db_error) +
    geom_point(aes(x = Error, y = Index), size = 0.2) + 
    geom_ribbon(aes(
      y = Index, 
      xmin = CI_inf,
      xmax = CI_sup
    ),
    alpha = 0.3,
    fill = "#FA9FB5"
    ) + 
    geom_vline(aes(xintercept = 0), col = "#DB15C1") +  #xlim(c(-15, 15)) +
    theme_classic() + xlab("Error (Mean pred - True value)")  
}

pred_risk = function(db, mod, hyperpost, input, threshold,
                     nb_samples = 1000, name = 'MagmaClust'){
  ## Predict at the input of interest
  pred = pred_magmaclust(
    data = db, 
    trained_model = mod, 
    hyperpost = hyperpost,
    grid_inputs = input, 
    get_full_cov = TRUE,
    plot = FALSE
  )
   
  ## Generate samples of the posterior at this input 
  samples = sample_magmaclust(
    pred_clust = pred,
    nb_samples = nb_samples)
  
  ## Compute the risk of going over the threshold 
  samples %>% 
    dplyr::filter(Input %in% input) %>% 
    summarise(Risk = mean(Output > threshold) * 100) %>% 
    mutate('ID' = unique(db$ID), 'Name' = name, .before = 1) %>% 
    return()
}

loop_pred_risk = function(db, mod, hyperpost, input, start_pred = 2,
                 nb_samples = 1000, name = 'MagmaClust'){
  if(hyperpost %>% is.null()){
    hyperposterior_clust(
      trained_model = mod,
      data = db,
      grid_inputs = unique(db$Input) %>%
        union(unique(mod$ini_args$data$Input)) %>%
        union(input) %>% 
        sort())
  }
  
  floop_i = function(i){
    cat('ID n°',i, '\n \n')
    
    ## Select on individual 
    ID_pred = i
    db_i = db %>% 
      filter(ID %in% ID_pred) %>% 
      filter(Input < start_pred)
    
    ## Change overweight threshold according to the Gender
    if(unique(db_i$Gender) == 'Female'){
      threshold = 22
    } else {
      threshold = 22.8
    }
    
    ## Remove the Gender column
    db_i = db_i %>% dplyr::select(- Gender)
    
    ## Compute the risk for the i-th individual
    risk_i = pred_risk(
      db = db_i,
      mod = mod,
      hyperpost = hyperpost,
      input = input,
      threshold = threshold,
      nb_samples = nb_samples,
      name = name
      )
  }
  db$ID %>%
    unique %>% 
    lapply(floop_i) %>%
    bind_rows() %>% 
    return()
}

##### Training #####

## when using full data, increase this value to 600
nb_training_indivs = 6
nb_total_indivs = n_distinct(db$ID)

## Split the dataset into training/testing sets
db_train = db %>% 
  filter(ID %in% unique(db$ID)[1:nb_training_indivs]) %>% 
  select(- Gender)

db_test = db %>% 
  filter(ID %in% unique(db$ID)[(nb_training_indivs+1):nb_total_indivs]) %>%
  select(- Gender)

## Train MagmaClust with 5 clusters 
## WARNING: retraining with full data would take a few hours 
#mod = train_magmaclust(db_train, nb_cluster = 5)

## To avoid retraining, the pre-trained model can be extracted instead
mod = read_rds('Training/trained_model.rds')

## Pre-compute all mean processes on a fine grid to speed up all experiments
grid_inputs = seq(0, 10.3, 0.1)
all_inputs = unique(db$Input) %>% union(grid_inputs) %>% sort()
h_post = hyperposterior_clust(trained_model = mod, grid_inputs = all_inputs)

##### Illustration #####

## Add Jens-Bailey illustration
raw_db_jb = read_csv('Prediction/Jens-Bailey_illu_3_subsets.csv')

## Select an illustrative example 
set.seed(7)
## In the paper, with full data Figure illu: 601, Figure forecast: 611
ID_input = unique(db$ID)[1] 

## Extract corresponding data from Jenss-Bayley's predictions
db_jb = raw_db_jb %>% mutate(Input = days / 365, Output = missingC) %>% 
  dplyr::select(Input, Output)

## Extract data from the testing individual
sub_db_ID = db %>% filter(ID %in% ID_input) %>% select(- Gender)
input_test = sub_db_ID %>% filter(Input < 6) %>% pull(Input)

## Differentiate observed and testing inputs
sub_db_pred = sub_db_ID %>% filter(Input %in% input_test)
sub_db_test = sub_db_ID %>% filter(!(Input %in% input_test))

## Splines prediction and plotting
spline_fit <- smooth.spline(sub_db_pred$Input, sub_db_pred$Output)
preds <- predict(spline_fit, grid_inputs) %>% as_tibble

ggplot(preds) + 
  geom_line(aes(x = x, y = y)) +
  geom_point(data = sub_db_test, aes(x = Input, y = Output), 
             col = 'red', size = 2) + 
  geom_point(data = sub_db_pred, aes(x = Input, y = Output), 
             col = 'black', size = 2) +
  xlab('Age (years)') + ylab('BMI') +  ylim(c(12,20)) + theme_classic()


## Jenss-Bayley plotting

ggplot(db_jb %>% slice(1:365)) + 
  geom_line(aes(x = Input, y = Output)) +
  geom_point(data = sub_db_test, aes(x = Input, y = Output),
             col = 'red', size = 2) + 
  geom_point(data = sub_db_pred, aes(x = Input, y = Output),
             col = 'black', size = 2) +
  xlab('Age (years)') + ylab('BMI') + ylim(c(12,20)) + theme_classic()


## MagmaClust prediction and plotting
pred = pred_magmaclust(sub_db_pred,
                       trained_model = mod,
                       grid_inputs = grid_inputs,
                       hyperpost = h_post,
                       get_hyperpost = T,
                       get_full_cov = T,
                       plot = F)

plot_magmaclust(pred, alpha_samples = 0.3, samples = T, size_data = 5) +
  geom_point(data = sub_db_test,
             aes(x = Input, y = Output), 
             col = 'red', 
             size = 2) + 
  geom_point(data = sub_db_pred, 
             aes(x = Input, y = Output), 
             col = 'black', 
             size = 2) + 
  xlab('Age (years)') + xlim(c(0,10.2)) + ylim(c(13,23)) +
  ylab('BMI') + labs(title=NULL) #+ scale_color_brewer(palette="Set3")

#### Illustration of forecasting 

forcast_to_10 = read_csv('Prediction/forecast_from2-6_to_10_with_Output.csv')
JB_to_10 = read_csv('Prediction/JB_forecast.csv') %>% mutate(Input = Days/365)

## Plot forecasting illustration
ggplot(JB_to_10) + geom_line(aes(x = Input, y = thresh_0to4)) +
  geom_point(data = sub_db_test,
             aes(x = Input, y = Output), 
             col = 'red', 
             size = 2) + 
  geom_point(data = sub_db_pred, 
             aes(x = Input, y = Output), 
             col = 'black', 
             size = 2) + 
  xlab('Age (years)')  + ylim(c(13,23)) +
  ylab('BMI') + theme_classic()
 
## Plot errors of MagmaClust predictions
forcast_to_10 %>% 
  filter(Method == '4to10_MagmaClust_5clusters') %>% 
  plot_error() + xlim(c(-15, 12))

##### Illustration risk overweight #####

## Select example and sub-sample of data

## illu men illu: "010-21686" overweight "010-21980",
## illu women: bad pred "020-66136" good pred "010-20748"
ID_input = unique(db$ID)[1]
sub_db_ID = db %>% filter(ID %in% ID_input) %>% select(- Gender)
input_test = sub_db_ID %>% filter(Input < 6) %>% pull(Input)

## Differentiate observed and testing inputs
sub_db_pred = sub_db_ID %>% filter(Input %in% input_test)
sub_db_test = sub_db_ID %>% filter(!(Input %in% input_test))

## MagmaClust predictions
pred = pred_magmaclust(sub_db_pred,
                       trained_model = mod,
                       grid_inputs = grid_inputs,
                       hyperpost = h_post,
                       get_hyperpost = T,
                       get_full_cov = T,
                       plot = F)

## Draw samples from the posterior
set.seed(4) #illu men: 4, illu women: 4, bad fit : 12
samples = sample_magmaclust(pred, nb_samples = 100)

## Extract IDs of samples above the overweight threshold
ID_ow = samples %>% 
  filter(Output > 22) %>%
  distinct(Sample) %>% 
  pull()

## Retrieve the full trajectories of those crossing the threshold
samples_ow = samples %>% filter(Sample %in% ID_ow)

plot_magmaclust(pred = pred, samples = T) +
  geom_line(data = samples_ow,
            aes(x = Input, y = Output, group = Sample),
            colour = "black", alpha = 1) +
  geom_hline(data = samples_ow,
              aes(yintercept = 22),
              col = "black", linetype = 'dashed', size = 0.8, alpha = 0.8) +
  geom_point(data = sub_db_pred,
             aes(x = Input, y = Output), col = 'black', size = 2) +
  geom_point(data = sub_db_test,
             aes(x = Input, y = Output),
             col = 'red', size = 2) +
  xlab('Age (years)') + ylab('BMI') + xlim(c(0,10.4)) + ylim(c(13,24))

##### Evaluation all individuals #####

## Extract IDs of testing individuals
ID_test = db %>% 
  filter(ID %in% unique(db$ID)[nb_training_indivs:nb_total_indivs]) %>%
  count(ID) %>%
  filter(n > 5) %>%
  pull(ID)

## Extract testing data
test_db = db %>% filter(ID %in% ID_test) %>% select(- Gender)

## Evaluate quality of MagmaClust predictions on testing individuals
res = eval(db = test_db, mod = mod, hyperpost = h_post, remove_random = F,
           start_input_test = 2, name = 'MagmaClust_5clusters_2to10')

## Averaged summary
summary_res = summarise_res(res, by = Method)

##### Missing data reconstruction #####

## Loop the predictions for different ratios of missing data 
floop = function(i){
  eval(db = db_test, mod = mod, remove_random = T,
       ratio_missing = i, hyperpost = h_post,
       name = 'Splines') %>% 
    mutate('Ratio_missing' = i) %>% 
    return()
}
res_missing_ratios = c(0.1, 0.25, 0.5, 0.75, 0.9) %>% 
  lapply(floop) %>% 
  bind_rows()

## Averaged summary
summary_res_missing_ratios = summarise_res(res_missing_ratios, by = Ratio_missing)

##### Prediction and plot of errors sorted by increasing uncertainty #####

## Forecasting on testing points for all testing individuals
pred_0to10years6 = loop_pred(db_test, mod, start_input_test = 6,
                            hyperpost = h_post,
                            name = '6to10_MagmaClust_5clusters')

pred_0to10years6 %>%
  group_by(ID) %>%
  reframe(
    Method = Method,
    ID = ID,
    MSE = (Output - Mean)^2,
    WCIC = 100*mean((Output < Mean+1.96*sqrt(Var)) & (Output > Mean-1.96*sqrt(Var)))
  ) %>%
  summarise_res()

## Missing data reconstruction for all testing individuals
missing_0to10years = loop_pred(db = db_test, 
                               mod = mod, remove_random = T,
                               ratio_missing = 0.5, hyperpost = h_post,
                               name = 'MagmaClust_5clusters_missing_50%')

missing_0to10years %>%
  group_by(ID) %>%
  reframe(
    Method = Method,
    ID = ID,
    MSE = (Output - Mean)^2,
    WCIC = 100*mean((Output < Mean+1.96*sqrt(Var)) & (Output > Mean-1.96*sqrt(Var)))
  ) %>%
  summarise_res()

##### Splines evaluation #####

## Extract testing individuals
ID_test = db %>% 
  filter(ID %in% unique(db$ID)[nb_training_indivs:nb_total_indivs]) %>%
  count(ID) %>%
  filter(n > 8) %>%
  pull(ID)

test_db = db %>% filter(ID %in% ID_test) %>% select(- Gender)

## Forecasting on testing points for all testing individuals
spline_pred = loop_pred(db = test_db, mod = '', hyperpost = '',
                       remove_random = F, start_input_test = 6,
                       name = 'Splines')

spline_pred %>%
  group_by(ID) %>%
  reframe(
    Method = Method,
    ID = ID,
    MSE = (Output - Mean)^2,
    WCIC = 100*mean((Output < Mean+1.96*sqrt(Var)) & (Output > Mean-1.96*sqrt(Var)))
    ) %>%
  summarise_res()

## Missing data reconstruction on testing points for all testing individuals
spline_miss = loop_pred(db = test_db, mod = '', hyperpost = '',
                       remove_random = T, ratio_missing = 0.5,
                       name = 'Splines')

spline_miss %>%
  group_by(ID) %>%
  reframe(
    Method = Method,
    ID = ID,
    MSE = (Output - Mean)^2,
    WCIC = 100*mean((Output < Mean+1.96*sqrt(Var)) & (Output > Mean-1.96*sqrt(Var)))
  ) %>%
  summarise_res()

##### Plot mean processes #####

## Tailored function to plot mean processes
plot_mean_proc = function(mod, add_data = T){
  pred = mod$hyperpost$pred
  
  db = data_allocate_cluster(mod)
  
  tib = c()
  for(i in names(pred)){
    tib = tib %>% bind_rows(
      pred[[i]] %>% mutate('Cluster' = i)
    )
  }
  
  gg = ggplot(tib) + 
    geom_line(aes(x = Input, y = Mean, col = Cluster), linetype = 'dashed') + 
    theme_classic() + xlim(c(0, 10.4))
  
  if(add_data){
    gg = gg +
    geom_point(
      data=db,
      aes(x = Input, y = Output,
          col = Cluster), 
      alpha = 0.2, size = 1) +
    scale_color_brewer(palette = "Set1")
    }
 
  return(gg)
}

plot_mean_proc(mod, F)

##### Compute risk overweight at 10 years ######

## Compute risk of overweight at 10 years from starting age 2, 4, 6 and 8 years.
res_risk = c()
for(i in c(2, 4, 6, 8)){
  name = paste0(i, 'to10_years')
  
  res_risk = res_risk %>% 
    bind_rows(
      loop_pred_risk(
        db = db,
        mod = mod,
        hyperpost = h_post,
        input = 10,
        start_pred = i,
        nb_samples = 10000,
        name = name
      )
    )
}

## Graph of overweight risk evaluation per forecasting period
risk = res_risk %>% left_join(db %>% 
  filter(Input > 9.5) %>%
  filter(Input < 10.5) %>%
  group_by(ID, Gender) %>% 
  summarise(Output = max(Output)) %>% 
  mutate('Status'  = if_else(
    Gender == 'Female',
    if_else(Output > 22, 'Overweight', 'Not overweight'),
    if_else(Output > 22.8, 'Overweight', 'Not overweight'))),
  by = 'ID') %>% 
  drop_na()

## Plot summary figure of overweight risks coloured by true status at 10 years
ggplot(risk) +
  geom_jitter(aes(x = Name, y = Risk, col = Status),
              position = position_jitterdodge(), alpha = 0.8) +
  # geom_boxplot(aes(x = Name, y = Risk, col = Status),
  #              outlier.shape = NA, alpha = 0.5) +
  theme_classic() + ylab('Overweight probability') +
  xlab('Forecasting period') + labs(col = "Status at 10 years")


## Number of children per group
risk %>% group_by(Status, Name) %>% summarise(N = n())

## Compute True/False positives/negatives for confusion matrices
risk %>% mutate(Eval = if_else(
  Risk>5, 
  if_else(Status == 'Overweight', 'TP', 'FP'),
  if_else(Status == 'Overweight', 'FN', 'TN')
  )) %>%
  group_by(Name) %>% 
  count(Eval) %>%  
  pivot_wider(names_from = Eval, values_from = n, values_fill = 0) %>% 
  mutate(Sensitivity = TP / (TP + FN), 
         Specificity = TN / (TN + FP), 
         Accuracy = (TP + TN) / (TP + FN + TN + FP))
 
         
