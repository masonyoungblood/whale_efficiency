#load libraries
library(brms)
library(lme4)
library(broom.mixed)
library(mice)
library(ggplot2)
library(ggstar)

#create production constraint model based on james et al. (2021)
#data requires two named columns: duration (of elements or intervals) and sequence (ID not type)
prod_null <- function(data){
  #get unique songs and durations
  unique_seqs <- unique(data$sequence)
  seq_durs <- as.numeric(sapply(unique_seqs, function(x){sum(data$duration[which(data$sequence == x)])}))
  
  #get shuffled versions of data
  sampled_rows <- sample(nrow(data))
  sampled_durs <- data$duration[sampled_rows]
  
  #double everything so the algorithm can loop back around if it needs to
  sampled_durs <- c(sampled_durs, sampled_durs)
  
  #set starting points and sampling window for cumulative duration calculation
  ind <- 1
  m <- 1
  window <- max(table(data$sequence))*100 #by default, 10 times the maximum signal length in units
  
  #create empty sequences list to fill
  sequences <- list()
  
  #iterate through unique songs and create sequences
  while(m <= length(unique_seqs)){
    #get cumulative durations from ind to the end of the sampling window
    cum_durs <- cumsum(sampled_durs[ind:(ind + window)])
    
    #set the targets as the value at which the cumulative durations surpass the desired duration
    rel_position <- min(which(cum_durs > seq_durs[m])) #relative position in the cumulative duration window
    abs_position <- ind + rel_position - 1 #absolute position in the sampled vectors
    
    #if syllable duration surpasses song duration at first position, do exception handling
    if(rel_position == 1){
      #add the single-unit sequence and iterate ind
      sequences[[m]] <- data.frame(duration = sampled_durs[ind], sequence = unique_seqs[m], length = rel_position)
      ind <- ind + 1
    } else{
      #if the difference between the cumulative duration at the target and the desired duration is less that half of the length of target
      if(cum_durs[rel_position] - seq_durs[m] < sampled_durs[abs_position]/2){
        #add the sequence and iterate ind
        sequences[[m]] <- data.frame(duration = sampled_durs[ind:abs_position], sequence = unique_seqs[m], length = rel_position)
        ind <- abs_position + 1
      } else{
        #otherwise, add abbreviated sequence and iterate ind
        sequences[[m]] <- data.frame(duration = sampled_durs[ind:(abs_position - 1)], sequence = unique_seqs[m], length = rel_position - 1)
        ind <- abs_position
      }
    }
    
    #iterate m
    m <- m + 1
  }
  
  #return combined sequences
  return(do.call(rbind, sequences))
}

#function that creates simulated datasets and fits the actual and simulated models
menz_fit <- function(data, n_sim_data = 10, use_mean = FALSE, rm_singles = TRUE, cores = 7){
  #store boolean of whether data are from multiple studies
  multiple_studies <- !is.data.frame(data)
  
  #if singletons are to be removed, then remove them
  #note that this happens upstream of the null model, so that it is constructed with the same sequences that are modelled in the real data
  if(rm_singles){
    #if data is a list of datasets (from different studies on the same species)
    if(multiple_studies){
      #loop through each dataset
      data <- lapply(1:length(data), function(x){
        #identify sequences with length of 1
        length_table <- data.frame(table(data[[x]]$sequence))
        to_rm <- which(length_table$Freq == 1)
        
        #if there are any, then remove them
        if(length(to_rm) > 0){
          return(data[[x]][-which(data[[x]]$sequence %in% length_table$Var1[to_rm]), ])
        } else{
          #otherwise, just return the original object
          return(data[[x]])
        }
      })
    } else{
      #otherwise, identify sequences with length of 1
      length_table <- data.frame(table(data$sequence))
      to_rm <- which(length_table$Freq == 1)
      
      #if there are any, then remove them
      if(length(to_rm) > 0){
        data <- data[-which(data$sequence %in% length_table$Var1[to_rm]), ]
      }
    }
  }
  
  #create simulated datasets
  if(n_sim_data > 0){
    if(multiple_studies){
      #if data is a list of datasets from different studies
      sim_data <- parallel::mclapply(1:n_sim_data, function(x){
        #compute the production null separately and concatenate with study as a new variable
        temp <- do.call(rbind, lapply(1:length(data), function(y){
          return(cbind(prod_null(data[[y]]), study = y))
        }))
        
        #and combine sequence and study into a string so sequences are unique
        temp$sequence <- paste0(temp$study, "-", temp$sequence)
        
        #if mean duration of each song is used then simplify
        if(use_mean){temp <- aggregate(duration ~ ., temp, FUN = mean)}
        
        #and return
        return(temp)
      }, mc.cores = cores)
    } else{
      #otherwise
      sim_data <- parallel::mclapply(1:n_sim_data, function(x){
        #compute the null on the full data
        temp <- prod_null(data)
        
        #if mean duration of each song is used then simplify
        if(use_mean){temp <- aggregate(duration ~ ., temp, FUN = mean)}
        
        #and return
        return(temp)
      }, mc.cores = cores)
    }
  }
  
  #if data is a list of datasets (from different studies on the same species)
  if(multiple_studies){
    #add study as a variable 
    data <- do.call(rbind, parallel::mclapply(1:length(data), function(x){
      #add sequence length to the original data
      data[[x]]$length <- sapply(1:nrow(data[[x]]), function(y){length(which(data[[x]]$sequence == data[[x]]$sequence[y]))})
      
      #add position to the original data
      data[[x]]$position <- unlist(lapply(unique(data[[x]]$sequence), function(y){
        #compute position using (n - 1)/(l - 1), as per james et al. (2021)
        seq_length <- length(which(data[[x]]$sequence == y))
        return((c(1:seq_length) - 1)/(seq_length - 1))
      }))
      
      #and add study as a variable
      return(cbind(data[[x]], study = x))
    }, mc.cores = cores))
    
    #and combine sequence and study into a string so sequences are unique
    data$sequence <- paste0(data$study, "-", data$sequence)
    
    #if mean duration of each song is used then simplify
    if(use_mean){data <- aggregate(duration ~ ., data, FUN = mean)}
  } else{
    #otherwise, add sequence length
    data$length <- sapply(1:nrow(data), function(x){length(which(data$sequence == data$sequence[x]))})
    
    #and add position
    data$position <- unlist(lapply(unique(data$sequence), function(x){
      #compute position using (n - 1)/(l - 1), as per james et al. (2021)
      seq_length <- length(which(data$sequence == x))
      return((c(1:seq_length) - 1)/(seq_length - 1))
    }))
    
    #if mean duration of each song is used then simplify
    if(use_mean){data <- aggregate(duration ~ ., data, FUN = mean)}
  }
  
  #run the full and reduced models, as well as a model of position (if means are not used)
  if(use_mean){
    if(multiple_studies){
      full_model <- lmer(log(duration) ~ log(length) + length + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      reduced_model <- lmer(log(duration) ~ log(length) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      reduced_scaled_model <- lmer(scale(log(duration)) ~ scale(log(length)) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    }
    if(!multiple_studies){
      full_model <- lm(log(duration) ~ log(length) + length, data = data)
      reduced_model <- lm(log(duration) ~ log(length), data = data)
      reduced_scaled_model <- lm(scale(log(duration)) ~ scale(log(length)), data = data)
    }
  } else{
    if(multiple_studies){
      full_model <- lmer(log(duration) ~ log(length) + length + (1|sequence) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      reduced_model <- lmer(log(duration) ~ log(length) + (1|sequence) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      reduced_scaled_model <- lmer(scale(log(duration)) ~ scale(log(length)) + (1|sequence) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      position_scaled_model <- lmer(scale(log(duration)) ~ scale(log(length)) + scale(position) + (1|sequence) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    }
    if(!multiple_studies){
      full_model <- lmer(log(duration) ~ log(length) + length + (1|sequence), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      reduced_model <- lmer(log(duration) ~ log(length) + (1|sequence), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      reduced_scaled_model <- lmer(scale(log(duration)) ~ scale(log(length)) + (1|sequence), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
      position_scaled_model <- lmer(scale(log(duration)) ~ scale(log(length)) + scale(position) + (1|sequence), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
    }
  }
  
  #fit a single model across the simulated datasets
  if(n_sim_data > 0){
    if(use_mean){
      if(multiple_studies){
        sim_full_model <- pool(lapply(1:n_sim_data, function(x){lmer(log(duration) ~ log(length) + length + (1|study), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
        sim_reduced_model <- pool(lapply(1:n_sim_data, function(x){lmer(log(duration) ~ log(length) + (1|study), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
        sim_reduced_scaled_model <- pool(lapply(1:n_sim_data, function(x){lmer(scale(log(duration)) ~ scale(log(length)) + (1|study), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
      }
      if(!multiple_studies){
        sim_full_model <- pool(lapply(1:n_sim_data, function(x){lm(log(duration) ~ log(length) + length, data = sim_data[[x]])}))
        sim_reduced_model <- pool(lapply(1:n_sim_data, function(x){lm(log(duration) ~ log(length), data = sim_data[[x]])}))
        sim_reduced_scaled_model <- pool(lapply(1:n_sim_data, function(x){lm(scale(log(duration)) ~ scale(log(length)), data = sim_data[[x]])}))
      }
    } else{
      if(multiple_studies){
        sim_full_model <- pool(lapply(1:n_sim_data, function(x){lmer(log(duration) ~ log(length) + length + (1|study/sequence), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
        sim_reduced_model <- pool(lapply(1:n_sim_data, function(x){lmer(log(duration) ~ log(length) + (1|study/sequence), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
        sim_reduced_scaled_model <- pool(lapply(1:n_sim_data, function(x){lmer(scale(log(duration)) ~ scale(log(length)) + (1|study/sequence), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
      }
      if(!multiple_studies){
        sim_full_model <- pool(lapply(1:n_sim_data, function(x){lmer(log(duration) ~ log(length) + length + (1|sequence), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
        sim_reduced_model <- pool(lapply(1:n_sim_data, function(x){lmer(log(duration) ~ log(length) + (1|sequence), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
        sim_reduced_scaled_model <- pool(lapply(1:n_sim_data, function(x){lmer(scale(log(duration)) ~ scale(log(length)) + (1|sequence), data = sim_data[[x]], REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}))
      }
    }
  }
  
  #return models
  if(use_mean){
    if(n_sim_data > 0){return(list(actual = list(full = full_model, reduced = reduced_model, reduced_scaled = reduced_scaled_model), prod_null = list(full = sim_full_model, reduced = sim_reduced_model, reduced_scaled = sim_reduced_scaled_model)))}
    if(n_sim_data == 0){return(list(actual = list(full = full_model, reduced = reduced_model, reduced_scaled = reduced_scaled_model)))}
  } else{
    if(n_sim_data > 0){return(list(actual = list(full = full_model, reduced = reduced_model, reduced_scaled = reduced_scaled_model, position_scaled = position_scaled_model), prod_null = list(full = sim_full_model, reduced = sim_reduced_model, reduced_scaled = sim_reduced_scaled_model)))}
    if(n_sim_data == 0){return(list(actual = list(full = full_model, reduced = reduced_model, reduced_scaled = reduced_scaled_model, position_scaled = position_scaled_model)))}
  }
}

#function for creating menzerath plot
menz_plot <- function(data, model, tokens = FALSE, intervals = FALSE, rm_singles = TRUE, original_axes = TRUE, unit = "s", color = "black", alpha = 1, effects_panel = 0.2, prop_dist_to_plot = 1, effects_axis = c(-0.8, 0.4), effects_breaks = 6, y_breaks = 2){
  #get intercept and effects for best fit line
  intercept <- summary(model$actual$reduced_scaled)$coef[1, 1]
  effect <- summary(model$actual$reduced_scaled)$coef[2, 1]
  
  #construct ylabel from function options
  ylab <- paste0(ifelse(!tokens, "Mean ", ""), ifelse(!intervals, "ED (", "IOI ("), unit, ")")
  
  #if multiple datasets then collapse into single
  if(!is.data.frame(data)){
    data <- data.frame(duration = unlist(lapply(1:length(data), function(x){data[[x]]$duration})),
                       sequence = unlist(lapply(1:length(data), function(x){paste0(x, "-", data[[x]]$sequence)})))
  }
  
  #if plot is averaged across sequences, aggregate by average and length
  if(!tokens){
    data <- data.frame(duration = aggregate(. ~ sequence, data, FUN = mean)[, 2], length = aggregate(. ~ sequence, data, FUN = length)[, 2])
  } else{
    #otherwise just add length as a variable
    temp <- aggregate(. ~ sequence, data, FUN = length)
    data$length <- temp[, 2][match(data$sequence, temp$sequence)]
  }
  
  #if there are any non-sequences (length less than 2) remove them
  if(rm_singles & length(which(data$length < 2)) > 0){data <- data[-which(data$length < 2), ]}
  
  #if the original axes of the data are to be plotted
  if(original_axes){
    #get original scaling from data
    og_scaled_x <- scale(log(data$length))
    og_scaled_y <- scale(log(data$duration))
    
    #function that predicts durations from lengths, based on intercepts, effects, and means and sds from scaling
    pred_durations <- function(intercept, effect, length, mean_log_length, sd_log_length, mean_log_dur, sd_log_dur){
      return(exp((intercept + effect*((log(length) - mean_log_length)/sd_log_length))*sd_log_dur + mean_log_dur))
    }
    
    #compute best fit line for plotting
    best_fit_data <- data.frame(x = seq(from = min(data$length), to = max(data$length), length.out = 1000))
    best_fit_data$y <- pred_durations(intercept, effect, length = best_fit_data$x,
                                      mean_log_length = attr(og_scaled_x, "scaled:center"), sd_log_length = attr(og_scaled_x, "scaled:scale"),
                                      mean_log_dur = attr(og_scaled_y, "scaled:center"), sd_log_dur = attr(og_scaled_y, "scaled:scale"))
  } else{
    #otherwise
    #overwrite original data with scaled and logged data
    data$length <- scale(log(data$length))
    data$duration <- scale(log(data$duration))
    
    #compute best fit line for plotting
    best_fit_data <- data.frame(x = seq(from = min(data$length), to = max(data$length), length.out = 1000))
    best_fit_data$y <- (best_fit_data$x*effect) + intercept
  }
  
  #construct plot
  plot <- ggplot(data, aes(x = length, y = duration)) +
    geom_point(colour = scales::alpha(color, alpha), size = 1, stroke = 0) +
    geom_line(data = best_fit_data, aes(x = x, y = y), colour = "black", linewidth = 0.5) +
    theme_linedraw(base_size = 6, base_family = "Avenir") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = y_breaks)) + 
    xlim(min(data$length), as.numeric(quantile(data$length, prop_dist_to_plot))) + 
    #ylim(min(data$duration), as.numeric(quantile(data$duration, prop_dist_to_plot))) + 
    xlab(paste0("# of ", ifelse(intervals, "Intervals", "Elements"), " in Sequence")) + ylab(ylab) + 
    theme(axis.text.y = element_text(angle = 90))
  
  #get bounds of confidence interval
  bounds <- confint(model$actual$reduced_scaled, parm = "scale(log(length))", method = "Wald")
  
  #construct effects plot
  effects_plot <- ggplot() + 
    geom_point(aes(x = 1, y = summary(model$actual$reduced_scaled)$coef[2, 1]), size = 1.5, stroke = 0, colour = scales::alpha(color, alpha)) + 
    geom_linerange(aes(x = 1, ymin = bounds[1], ymax = bounds[2]), colour = scales::alpha(color, alpha)) + 
    theme_linedraw(base_size = 6, base_family = "Avenir") + 
    geom_hline(yintercept = 0, lty = "dotted") + 
    scale_y_continuous(name = "Slope", position = "right", limits = effects_axis, labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = effects_breaks)) + 
    scale_x_continuous(name = "...", breaks = c(1), labels = c("A"), limits = c(0.5, 1.5)) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(colour = "white"),
          axis.text.x = element_text(colour = "white"))
  
  #return plot
  return(cowplot::plot_grid(plot, effects_plot, rel_widths = c(1, effects_panel)))
}

#function for comparing the strength of menzerath's law in mysticetes and odonticetes, and between duration and interval data
#data should be a nested list, where top levels are species, and second levels are individual studies (if applicable)
#group should be a vector of 0/1, 0 for mysticetes and 1 for odontocetes
#type should be a vector of 0/1, 0 for duration and 1 for interval
menz_compare <- function(data, group, type, rm_singles = TRUE, cores = 7){
  #combine datasets
  data <- do.call(rbind, parallel::mclapply(1:length(data), function(x){
    #store temporary version of data to modify
    temp <- data[[x]]
    
    #if there are multiple studies
    if(!is.data.frame(temp)){
      #iterate through them and add measures
      temp <- do.call(rbind, lapply(1:length(temp), function(y){
        #add sequence length to the original data
        temp[[y]]$length <- sapply(1:nrow(temp[[y]]), function(z){length(which(temp[[y]]$sequence == temp[[y]]$sequence[z]))})
        
        #and add study as a variable
        return(cbind(temp[[y]], study = y, species = x, group = group[x], type = type[x]))
      }))
      
      #if there are any non-sequences (length less than 2) remove them
      if(rm_singles & length(which(temp$length < 2)) > 0){temp <- temp[-which(temp$length < 2), ]}
      
      #and concatenate what needs to be concatenated
      temp$study <- paste0(temp$species, "-", temp$study)
      temp$sequence <- paste0(temp$study, "-", temp$sequence)
      
      #and add position
      temp$position <- scale(unlist(lapply(unique(temp$sequence), function(x){
        #compute position using (n - 1)/(l - 1), as per james et al. (2021)
        seq_length <- length(which(temp$sequence == x))
        return((c(1:seq_length) - 1)/(seq_length - 1))
      })))
      
      #log transform and scale duration and length
      temp$duration <- scale(log(temp$duration))
      temp$length <- scale(log(temp$length))
    } else{
      #add sequence length to the original data
      temp$length <- sapply(1:nrow(temp), function(y){length(which(temp$sequence == temp$sequence[y]))})
      
      #if there are any non-sequences (length less than 2) remove them
      if(rm_singles & length(which(temp$length < 2)) > 0){temp <- temp[-which(temp$length < 2), ]}
      
      #add study, group, and type
      temp$study <- 1
      temp$species <- x
      temp$group <- group[x]
      temp$type <- type[x]
      
      #and concatenate what needs to be concatenated
      temp$study <- paste0(temp$species, "-", temp$study)
      temp$sequence <- paste0(temp$study, "-", temp$sequence)
      
      #and add scaled position
      temp$position <- scale(unlist(lapply(unique(temp$sequence), function(x){
        #compute position using (n - 1)/(l - 1), as per james et al. (2021)
        seq_length <- length(which(temp$sequence == x))
        return((c(1:seq_length) - 1)/(seq_length - 1))
      })))
      
      #log transform and scale duration and length
      temp$duration <- scale(log(temp$duration))
      temp$length <- scale(log(temp$length))
    }
    
    #and return
    return(temp)
  }, mc.cores = cores))
  
  #run base model and model with position
  base_model <- lmer(duration ~ length + length:group + length:type + (1|sequence) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
  position_model <- model <- lmer(duration ~ length + position + length:group + length:type + position:group + position:type + (1|sequence) + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
  
  #return model
  return(list(base = base_model, position = position_model))
}

#read human phonemic data from the doreco dataset
read_phonemes <- function(path, interpolate = FALSE){
  #read in and format csv file
  data <- read.csv(path, header = FALSE)
  if(length(which(is.na(data[1, ]))) > 0){
    data <- data[, 1:(min(which(is.na(data[1, ]))) - 1)]
  }
  colnames(data) <- data[1, ]
  data <- data[-1, ]
  data$duration <- as.numeric(data$end) - as.numeric(data$start)
  
  #remove exceptional speech events, such as pauses, singing, and disfluencies: https://doreco.huma-num.fr/HowTo
  data <- data[-which(data$ph %in% c("<<fp>>", "<<fs>>", "<<pr>>", "<<fm>>", "<<sg>>", "<<bc>>", "<<id>>", "<<on>>", "<<wip>>", "<<ui>>", "<p:>")), ]
  
  #if interpolation of raw durations with median durations of phonemes is desired (to assess robustness of interpolation of humpback data)
  if(interpolate){
    #get unique phonemes
    unique_phonemes <- data.frame(phoneme = unique(data$ph))
    
    #for each unique phoneme
    unique_phonemes$med_dur <- sapply(unique_phonemes$phoneme, function(x){
      #get the word IDs corresponding to that word
      unique_ids <- unique(data$ph_ID[which(data$ph == x)])
      
      #return the median duration across all tokens of that word
      return(median(sapply(unique_ids, function(y){sum(data$duration[which(data$ph_ID == y)])})))
    })
    
    #replace real durations with median durations of types
    data$duration <- unlist(parallel::mclapply(data$ph, function(x){unique_phonemes$med_dur[which(unique_phonemes$phoneme == x)]}, mc.cores = 7))
  }
  
  #phonemes within words
  word_level <- data[, match(c("duration", "wd_ID"), colnames(data))]
  colnames(word_level) <- c("duration", "sequence")
  word_level$sequence <- as.numeric(factor(word_level$sequence))
  
  #words within sentences
  inds <- rle(data$tx)
  inds_end <- cumsum(inds$lengths)
  inds_start <- c(1, lag(inds_end)[-1] + 1)
  sentence_level <- do.call(rbind, lapply(1:length(inds_end), function(x){
    temp <- data[inds_start[x]:inds_end[x], ]
    data.frame(duration = aggregate(duration ~ wd_ID, temp, sum)$duration, sequence = x)
  }))
  
  #return both
  return(list(words = word_level, sentences = sentence_level))
}

#extract effects from models
extract_freq_effects <- function(models){
  #get vector of L95%, M, and H95% for prod null (reduced version)
  prod_null_est <- summary(models$prod_null$reduced_scaled, conf.int = TRUE)
  prod_null_est <- as.numeric(prod_null_est[which(prod_null_est$term == "scale(log(length))"), c(7, 2, 8)])
  
  #get vector of L95%, M, and H95% for length
  length_est <- summary(models$actual$reduced_scaled)$coef
  length_est <- length_est[which(rownames(length_est) == "scale(log(length))"), 1]
  length_est <- as.numeric(c(length_est, c(confint(models$actual$reduced_scaled, parm = "scale(log(length))", method = "Wald"))))
  length_est <- length_est[c(2, 1, 3)]
  
  #get vector of L95%, M, and H95% for position
  position_est <- summary(models$actual$position_scaled)$coef
  position_est <- position_est[which(rownames(position_est) == "scale(position)"), 1]
  position_est <- as.numeric(c(position_est, c(confint(models$actual$position_scaled, parm = "scale(position)", method = "Wald"))))
  position_est <- position_est[c(2, 1, 3)]
  
  #return both
  return(list(length = length_est, position = position_est, prod_null = prod_null_est))
}

#create label for phylogeny plots
label_maker <- function(data, rm_singles = TRUE, intervals = FALSE){
  #store boolean of whether data are from multiple studies
  multiple_studies <- !is.data.frame(data)
  
  #set number of studies to zero by default
  n_studies <- 1
  
  #if multiple studies are in the data
  if(multiple_studies){
    #store actual number of studies
    n_studies <- length(data)
    
    #collapse data
    data <- data.frame(duration = unlist(lapply(1:length(data), function(x){data[[x]]$duration})),
                       sequence = unlist(lapply(1:length(data), function(x){paste0(x, "-", data[[x]]$sequence)})))
  }
  
  #otherwise just add length as a variable
  temp <- aggregate(. ~ sequence, data, FUN = length)
  data$length <- temp[, 2][match(data$sequence, temp$sequence)]
  
  #if there are any non-sequences (length less than 2) remove them
  if(rm_singles & length(which(data$length < 2)) > 0){data <- data[-which(data$length < 2), ]}
  
  #create text label
  return(paste0(formatC(nrow(data), big.mark = ","), " ", ifelse(intervals, "Intervals", "Elements"), "\n", 
                formatC(length(unique(data$sequence)), big.mark = ","), " ", "Sequences", "\n",
                n_studies, " ", ifelse(multiple_studies, "Studies", "Study")))
}

#create function for extracting model effects and confidence intervals
table_output <- function(model, w_position = FALSE, zla = FALSE){
  #prod_null_all <- as.numeric(summary(model$prod_null$reduced_scaled, conf.int = TRUE, method = "Wald")[2, c(2, 7, 8)])
  
  if(!zla){
    if(w_position){
      length_effect <- summary(model$actual$position_scaled)$coefficients[2, 1]
      length_confints <- as.numeric(confint(model$actual$position_scaled, parm = "scale(log(length))", method = "Wald"))
      position_effect <- summary(model$actual$position_scaled)$coefficients[3, 1]
      position_confints <- as.numeric(confint(model$actual$position_scaled, parm = "scale(position)", method = "Wald"))
      output <- c(length_effect, length_confints, position_effect, position_confints)
    } else{
      length_effect <- summary(model$actual$reduced_scaled)$coefficients[2, 1]
      length_confints <- as.numeric(confint(model$actual$reduced_scaled, parm = "scale(log(length))", method = "Wald"))
      output <- c(length_effect, length_confints)
    }
  }
  
  if(zla){
    effect <- summary(model)$coefficients[2, 1]
    confints <- confint(model, parm = "count", method = "Wald")
    output <- c(effect, confints)
  }
  
  return(output)
}

#fit zla models
zla_fit <- function(data){
  #if data are from multiple studies, combine them (scaling within studies to enable comparison despite different sample sizes)
  if(!is.data.frame(data)){
    data <- do.call(rbind, lapply(1:length(data), function(x){
      data[[x]]$count <- scale(sapply(1:nrow(data[[x]]), function(y){length(which(data[[x]]$type == data[[x]]$type[y]))}))
      data[[x]]$study <- x
      data[[x]]
    }))
    lme4::lmer(scale(log(duration)) ~ count + (1|type) + (1|study), data = data, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))
  } else{
    #otherwise run normally
    data$count <- scale(sapply(1:nrow(data), function(x){length(which(data$type == data$type[x]))}))
    lme4::lmer(scale(log(duration)) ~ count + (1|type), data = data, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))
  }
}

#function for creating zla plot
zla_plot <- function(data, model, original_axes = TRUE, color = "black", aggregated = TRUE, alpha = 1, effects_panel = 0.2, effects_breaks = 8, ylims = NULL, effects_axis = c(-1, 1.2), y_breaks = 2){
  #structure data for plotting
  if(!is.data.frame(data)){
    data <- do.call(rbind, lapply(1:length(data), function(x){
      data[[x]]$raw_count <- sapply(1:nrow(data[[x]]), function(y){length(which(data[[x]]$type == data[[x]]$type[y]))})
      data[[x]]$count <- scale(sapply(1:nrow(data[[x]]), function(y){length(which(data[[x]]$type == data[[x]]$type[y]))}))
      data[[x]]$study <- x
      data[[x]]
    }))
  } else{
    data$raw_count <- sapply(1:nrow(data), function(x){length(which(data$type == data$type[x]))})
    data$count <- scale(sapply(1:nrow(data), function(x){length(which(data$type == data$type[x]))}))
  }
  
  #get intercept and effects for best fit line
  intercept <- summary(model)$coef[1, 1]
  effect <- summary(model)$coef[2, 1]
  
  #if the original axes of the data are to be plotted
  if(original_axes){
    #get original scaling from data
    og_scaled_y <- scale(log(data$duration))
    
    #function that predicts durations from lengths, based on intercepts, effects, and means and sds from scaling
    pred_durations <- function(intercept, effect, count, mean_log_dur, sd_log_dur){
      return(exp((intercept + effect*count)*sd_log_dur + mean_log_dur))
    }
    
    #compute best fit line for plotting
    best_fit_data <- data.frame(x = seq(from = min(data$count), to = max(data$count), length.out = 1000))
    best_fit_data$y <- pred_durations(intercept, effect, count = best_fit_data$x,
                                      mean_log_dur = attr(og_scaled_y, "scaled:center"), sd_log_dur = attr(og_scaled_y, "scaled:scale"))
    best_fit_data$x <- seq(from = min(data$raw_count), to = max(data$raw_count), length.out = 1000)
  } else{
    #otherwise
    #overwrite original data with scaled and logged data
    data$duration <- scale(log(data$duration))
    
    #compute best fit line for plotting
    best_fit_data <- data.frame(x = seq(from = min(data$count), to = max(data$count), length.out = 1000))
    best_fit_data$y <- (best_fit_data$x*effect) + intercept
  }
  
  #construct plot, aggregating durations into types to match menzerath plot
  #with option to specify custom xlim
  if(aggregated){
    if(original_axes){
      data <- aggregate(duration ~ type + raw_count, data, mean)
      colnames(data) <- c("type", "count", "duration")
    } else{
      data <- aggregate(duration ~ type + count, data, mean)
      colnames(data) <- c("type", "count", "duration")
    }
  }
  if(is.null(ylims)){
    plot <- ggplot(data, aes(x = count, y = duration)) +
      geom_point(colour = scales::alpha(color, alpha), size = 1, stroke = 0) +
      geom_line(data = best_fit_data, aes(x = x, y = y), colour = "black", linewidth = 0.5) +
      theme_linedraw(base_size = 6, base_family = "Avenir") +
      xlab("Count") + ylab("Mean ED (s)") + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = y_breaks)) + 
      theme(axis.text.y = element_text(angle = 90))
  } else{
    plot <- ggplot(data, aes(x = count, y = duration)) +
      geom_point(colour = scales::alpha(color, alpha), size = 1, stroke = 0) +
      geom_line(data = best_fit_data, aes(x = x, y = y), colour = "black", linewidth = 0.5) +
      theme_linedraw(base_size = 6, base_family = "Avenir") +
      xlab("Count") + ylab("Mean ED (s)") + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = y_breaks), limits = c(ylims[1], ylims[2])) + 
      theme(axis.text.y = element_text(angle = 90))
  }
  
  #get bounds of confidence interval
  bounds <- confint(model, parm = "count", method = "Wald")
  
  #construct effects plot
  effects_plot <- ggplot() + 
    geom_point(aes(x = 1, y = summary(model)$coef[2, 1]), size = 1.5, stroke = 0, colour = scales::alpha(color, alpha)) + 
    geom_linerange(aes(x = 1, ymin = bounds[1], ymax = bounds[2]), colour = scales::alpha(color, alpha)) + 
    theme_linedraw(base_size = 6, base_family = "Avenir") + 
    geom_hline(yintercept = 0, lty = "dotted") + 
    scale_y_continuous(name = "Slope", position = "right", limits = effects_axis, labels = scales::number_format(accuracy = 0.1), breaks = scales::pretty_breaks(n = effects_breaks)) + 
    scale_x_continuous(name = "...", breaks = c(1), labels = c("A"), limits = c(0.5, 1.5)) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(colour = "white"),
          axis.text.x = element_text(colour = "white"))
  
  #return plot
  return(cowplot::plot_grid(plot, effects_plot, rel_widths = c(1, effects_panel)))
}

#create label for phylogeny plots
zla_label_maker <- function(data){
  #store boolean of whether data are from multiple studies
  multiple_studies <- !is.data.frame(data)
  
  #set number of studies to zero by default
  n_studies <- 1
  
  #if multiple studies are in the data
  if(multiple_studies){
    #store actual number of studies
    n_studies <- length(data)
    
    #collapse data
    data <- data.frame(type = unlist(lapply(1:length(data), function(x){paste0(x, "-", data[[x]]$type)})),
                       duration = unlist(lapply(1:length(data), function(x){data[[x]]$duration})))
  }
  
  #create text label
  return(paste0(formatC(nrow(data), big.mark = ","), " Elements", "\n", 
                formatC(length(unique(data$type)), big.mark = ","), " ", "Types", "\n",
                n_studies, " ", ifelse(multiple_studies, "Studies", "Study")))
}

#read human phonemic data from the doreco dataset for zla
zla_read_phonemes <- function(path){
  #read in and format csv file
  data <- read.csv(path, header = FALSE)
  if(length(which(is.na(data[1, ]))) > 0){
    data <- data[, 1:(min(which(is.na(data[1, ]))) - 1)]
  }
  colnames(data) <- data[1, ]
  data <- data[-1, ]
  data$duration <- as.numeric(data$end) - as.numeric(data$start)
  
  #remove exceptional speech events, such as pauses, singing, and disfluencies: https://doreco.huma-num.fr/HowTo
  data <- data[-which(data$ph %in% c("<<fp>>", "<<fs>>", "<<pr>>", "<<fm>>", "<<sg>>", "<<bc>>", "<<id>>", "<<on>>", "<<wip>>", "<<ui>>", "<p:>")), ]
  
  #phonemes
  phonemes_zla <- data.frame(type = data$ph, duration = data$duration)
  
  #words
  words_zla <- aggregate(duration ~ wd_ID, data, sum)
  words_zla$wd <- data$wd[match(words_zla$wd_ID, data$wd_ID)]
  words_zla <- data.frame(type = words_zla$wd, duration = words_zla$duration)
  
  #return both
  return(list(phonemes = phonemes_zla, words = words_zla))
}

#function for assessing zipf's law of abbreviation across species
#data should be a nested list, where top levels are species, and second levels are individual studies (if applicable)
#group should be a vector of 0/1, 0 for mysticetes and 1 for odontocetes
zipf_compare <- function(data, group, include_group = TRUE, cores = 7){
  #combine datasets
  data <- do.call(rbind, parallel::mclapply(1:length(data), function(x){
    #store temporary version of data to modify
    temp <- data[[x]]
    
    #if data are from multiple studies, combine them (scaling within studies to enable comparison despite different sample sizes)
    if(!is.data.frame(temp)){
      temp <- do.call(rbind, lapply(1:length(temp), function(y){
        temp[[y]]$count <- scale(sapply(1:nrow(temp[[y]]), function(z){length(which(temp[[y]]$type == temp[[y]]$type[z]))}))
        temp[[y]]$duration <- scale(log(temp[[y]]$duration))
        temp[[y]]$species <- x
        temp[[y]]$study <- y
        if(include_group){temp[[y]]$group <- group[x]}
        return(temp[[y]])
      }))
      
      #concatenate what needs to be concatenated
      temp$study <- paste0(temp$species, "-", temp$study)
    } else{
      #otherwise run normally
      temp$count <- scale(sapply(1:nrow(temp), function(x){length(which(temp$type == temp$type[x]))}))
      temp$duration <- scale(log(temp$duration))
      temp$species <- x
      temp$study <- 1
      if(include_group){temp$group <- group[x]}
      
      #concatenate what needs to be concatenated
      temp$study <- paste0(temp$species, "-", temp$study)
    }
    
    #and return
    return(temp) 
  }, mc.cores = cores))
  
  #run model
  if(include_group){model <- lmer(duration ~ count + count:group + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}
  if(!include_group){model <- lmer(duration ~ count + (1|study), data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))}
  
  #return model
  return(model)
}
