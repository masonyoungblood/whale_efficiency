#set working directory
setwd("~/Documents/Work/Spring 2024/Whale ZLA/whale_song_efficiency/docs")

# SPERM WHALES ------------------------------------------------------------

sperm_gero_2016 <- readxl::read_xlsx("data/sperm_gero_2016.xlsx", sheet = 2)
sperm_gero_2016 <- do.call(rbind, parallel::mclapply(1:nrow(sperm_gero_2016), function(x){
  temp <- sperm_gero_2016[x, c(2:10)]
  temp <- temp[which(temp > 0)]
  return(data.frame(duration = as.numeric(temp), sequence = x))
}, mc.cores = 7))

sperm_hersh_2022 <- read.csv("data/sperm_hersh_2022.csv")
sperm_hersh_2022 <- do.call(rbind, parallel::mclapply(1:nrow(sperm_hersh_2022), function(x){
  temp <- sperm_hersh_2022[x, c(13:42)]
  temp <- temp[which(temp > 0)]
  return(data.frame(duration = as.numeric(temp), sequence = x))
}, mc.cores = 7))

sperm_vachon_2022 <- suppressWarnings(readxl::read_xlsx("data/sperm_vachon_2022.xlsx", sheet = 2))
sperm_vachon_2022 <- do.call(rbind, lapply(c(1:nrow(sperm_vachon_2022))[-which(rowSums(sperm_vachon_2022[, c(15:54)], na.rm = TRUE) == 0)], function(x){
  temp <- sperm_vachon_2022[x, c(15:54)]
  temp <- temp[which(temp > 0 & !is.na(temp))]
  return(data.frame(duration = as.numeric(temp), sequence = x))
}))

sperm_data <- list(gero_2016 = sperm_gero_2016, hersh_2022 = sperm_hersh_2022, vachon_2022 = sperm_vachon_2022)
save(sperm_data, file = "data/processed/sperm_data.RData")

# FIN WHALES --------------------------------------------------------------

fin_romagosa_2024 <- read.csv("data/fin_romagosa_2024/INI_all_corrected_by_romagosa.csv")
fin_romagosa_2024 <- data.frame(duration = fin_romagosa_2024$INI, sequence = as.numeric(factor(fin_romagosa_2024$date)))

fin_wood_2022 <- do.call(rbind, parallel::mclapply(2:8, function(x){
  temp <- readxl::read_xlsx("data/fin_wood_2022.xlsx", sheet = x)
  song_starts <- which(!is.na(temp$Pattern))
  song_stops <- c(song_starts[2:length(song_starts)]-2, nrow(temp))
  song_id <- rep(NA, nrow(temp))
  for(y in 1:length(song_starts)){song_id[song_starts[y]:song_stops[y]] <- y}
  temp$song_id <- song_id
  return(do.call(rbind, lapply(unique(temp$song_id)[-1], function(y){
    temp_temp <- as.numeric(format(temp$`IPI (s)`[which(temp$song_id == y)], "%OS"))
    return(data.frame(duration = temp_temp, sequence = paste0(x, "-", y)))
  })))
}, mc.cores = 7))
fin_wood_2022$sequence <- as.numeric(factor(fin_wood_2022$sequence))

fin_best_2022 <- read.csv("data/fin_best_2022.csv")
fin_best_2022 <- do.call(rbind, parallel::mclapply(unique(fin_best_2022$songid), function(x){
  temp <- fin_best_2022$IPI[which(fin_best_2022$songid == x)]
  temp <- temp[which(!is.na(temp))]
  return(data.frame(duration = as.numeric(temp), sequence = x))
}, mc.cores = 7))

fin_data <- list(wood_2022 = fin_wood_2022, best_2022 = fin_best_2022, romagosa_2024 = fin_romagosa_2024)
save(fin_data, file = "data/processed/fin_data.RData")

# KILLER WHALES -----------------------------------------------------------

killer_sharpe_2017 <- readxl::read_xlsx("data/killer_sharpe_2017.xlsx")
killer_sharpe_2017 <- do.call(rbind, lapply(1:nrow(killer_sharpe_2017), function(x){
  temp <- na.omit(as.numeric(killer_sharpe_2017[x, -1]))
  data.frame(duration = temp, sequence = rep(x, length(temp)))
}))

killer_data <- killer_sharpe_2017
save(killer_data, file = "data/processed/killer_data.RData")

killer_selbmann_2023 <- read.csv("data/killer_selbmann_2023.csv")
song_gaps <- which(killer_selbmann_2023$pause > 1.72) #pause gap based on analysis here: https://www.nature.com/articles/s41598-023-48349-1#Sec4
killer_selbmann_2023$sequence <- NA
killer_selbmann_2023$sequence[1:song_gaps[1]] <- 1
for(x in 1:(length(song_gaps)-1)){
  killer_selbmann_2023$sequence[(song_gaps[x]+1):song_gaps[x+1]] <- x+1
}
killer_selbmann_2023$sequence[(max(song_gaps)+1):nrow(killer_selbmann_2023)] <- max(killer_selbmann_2023$sequence, na.rm = TRUE)+1
killer_selbmann_2023$duration <- killer_selbmann_2023$end_time-killer_selbmann_2023$start_time
killer_selbmann_2023 <- killer_selbmann_2023[, c(7, 6)]

killer_sequence_data <- killer_selbmann_2023
save(killer_sequence_data, file = "data/processed/killer_sequence_data.RData")

# BLUE WHALES -------------------------------------------------------------

#explored before: https://www.animalbehaviorandcognition.org/uploads/journals/53/2%20Jolliffe_et_al_ABC_10(3).pdf
#this parsing is complicated and follows this step by step: https://royalsocietypublishing.org/doi/full/10.1098/rsos.180241#d1e2938
blue_lewis_2018 <- do.call(rbind, lapply(2:43, function(x){
  temp <- readxl::read_xlsx("data/blue_lewis_2018.xlsx", sheet = x)
  temp <- temp[which(temp[, 7] == "tagged whale"), ]
  temp <- data.frame(type = temp$`Call Type`, start = temp$`Call start time`, end = temp$`Call end time`)
  return(temp)
}))
blue_lewis_2018 <- blue_lewis_2018[-which(is.na(blue_lewis_2018$end)), ]
blue_lewis_2018 <- blue_lewis_2018[order(blue_lewis_2018$start), ]
blue_lewis_2018$pause <- c(0, as.numeric(blue_lewis_2018$start[2:nrow(blue_lewis_2018)]-blue_lewis_2018$start[1:(nrow(blue_lewis_2018)-1)]))
blue_lewis_2018$type[which(blue_lewis_2018$type == "A NE Pacific")] <- "A"
blue_lewis_2018$type[which(blue_lewis_2018$type == "B NE Pacific")] <- "B"

#remove all D calls, since they never appear in phrases by definition
blue_lewis_2018 <- blue_lewis_2018[-which(blue_lewis_2018$type == "D"), ]

#groups of A and B with pauses less than 70 are repetitive phrases: https://royalsocietypublishing.org/doi/full/10.1098/rsos.180241#d1e2938
#the original paper they cite used 120 seconds, but found that 70 seconds was the avg, and they appear to have used that too...
phrases <- c(1)
for(x in 2:nrow(blue_lewis_2018)){
  if(blue_lewis_2018$type[x] %in% c("A", "B") & blue_lewis_2018$type[x-1] %in% c("A", "B") & blue_lewis_2018$pause[x] < (70)){
  #if(blue_lewis_2018$type[x] %in% c("A", "B") & blue_lewis_2018$type[x-1] %in% c("A", "B") & blue_lewis_2018$pause[x] < (70 + 29)){
  #if(blue_lewis_2018$type[x] %in% c("A", "B") & blue_lewis_2018$type[x-1] %in% c("A", "B") & blue_lewis_2018$pause[x] < (120)){
    phrases <- c(phrases, phrases[x-1])
  } else{
    phrases <- c(phrases, phrases[x-1]+1)
  }
}
blue_lewis_2018$sequence <- phrases
blue_lewis_2018$length <-  sapply(1:nrow(blue_lewis_2018), function(x){length(which(blue_lewis_2018$phrase == blue_lewis_2018$phrase[x]))})
blue_lewis_2018$duration <- as.numeric(blue_lewis_2018$end - blue_lewis_2018$start)
blue_lewis_2018 <- blue_lewis_2018[, c(7, 5)]

blue_data <- blue_lewis_2018
save(blue_data, file = "data/processed/blue_data.RData")

# HUMPBACK WHALES ---------------------------------------------------------

phrases <- readLines("data/humpback_owen_2019/phrases.txt")
units <- read.csv("data/humpback_owen_2019/units.txt", sep = "\t")
units <- aggregate(Dur ~ Sound, data = units, FUN = median)
units$Sound <- gsub("\\.", "-", units$Sound)

humpback_owen_2019 <- do.call(rbind, lapply(1:length(phrases), function(x){
  temp <- strsplit(phrases[x], ", ")[[1]][-c(1:4)]
  match(temp, units$Sound)
  return(data.frame(duration = units$Dur[match(temp, units$Sound)], sequence = x))
}))
humpback_owen_2019$duration[which(is.na(humpback_owen_2019$duration))] <- median(humpback_owen_2019$duration, na.rm = TRUE)

humpback_schall_2021 <- list.files("data/humpback_schall_2021", "*.txt")
humpback_schall_2021 <- do.call(rbind, lapply(1:length(humpback_schall_2021), function(x){
  temp <- read.csv(paste0("data/humpback_schall_2021/", humpback_schall_2021[x]), sep = "\t")
  data.frame(duration = temp$End.Time..s.-temp$Begin.Time..s., sequence = x)
}))

humpback_schall_2022 <- list.files("data/humpback_schall_2022", "*.txt")
humpback_schall_2022 <- humpback_schall_2022[-grep("sel.", humpback_schall_2022)]
humpback_schall_2022 <- do.call(rbind, lapply(1:length(humpback_schall_2022), function(x){
  temp <- read.csv(paste0("data/humpback_schall_2022/", humpback_schall_2022[x]), sep = "\t")
  data.frame(duration = temp$End.Time..s.-temp$Begin.Time..s., sequence = x)
}))

humpback_data <- humpback_owen_2019
save(humpback_data, file = "data/processed/humpback_data.RData")

humpback_phrase_data <- list(schall_2021 = humpback_schall_2021, schall_2022 = humpback_schall_2022)
save(humpback_phrase_data, file = "data/processed/humpback_phrase_data.RData")

# MINKE WHALES ------------------------------------------------------------

#sequences of calls, which are different from the bioduck or whatever that I requested from other labs
minke_martin_2022 <- readxl::read_xlsx("data/minke_martin_2022.xlsx", sheet = 3)

#parse into intervals
minke_martin_2022 <- do.call(rbind, lapply(1:length(unique(minke_martin_2022$TrackNum)), function(x){
  temp <- minke_martin_2022[which(minke_martin_2022$TrackNum == unique(minke_martin_2022$TrackNum)[x]), ]
  intervals <- diff(temp$JulianTimeOfEmission)*86400 #in seconds
  below <- which(intervals <= (0.63+0.36+0.36+0.36)*60) #three standard deviations above the mean of state 2 in seconds
  if(length(below) > 0){
    
    #https://www.tutorialspoint.com/how-to-split-a-vector-into-smaller-vectors-of-consecutive-values-in-r
    seq_inds <- split(below, cumsum(c(1, diff(below) != 1)))
    
    #get intervals
    temp <- do.call(rbind, lapply(1:length(seq_inds), function(y){
      data.frame(duration = intervals[seq_inds[[y]]], sequence = paste0(x, "-", y))
    }))
    return(temp)
  }
}))
minke_martin_2022$sequence <- as.numeric(factor(minke_martin_2022$sequence))

#parse into elements
minke_martin_2022_elements <- readxl::read_xlsx("data/minke_martin_2022.xlsx", sheet = 3)
minke_martin_2022_elements <- do.call(rbind, lapply(1:length(unique(minke_martin_2022_elements$TrackNum)), function(x){
  temp <- minke_martin_2022_elements[which(minke_martin_2022_elements$TrackNum == unique(minke_martin_2022_elements$TrackNum)[x]), ]
  intervals <- diff(temp$JulianTimeOfEmission)*86400 #in seconds
  below <- which(intervals <= (0.63+0.36+0.36+0.36)*60) #three standard deviations above the mean of state 2 in seconds
  if(length(below) > 0){
    
    #https://www.tutorialspoint.com/how-to-split-a-vector-into-smaller-vectors-of-consecutive-values-in-r
    seq_inds <- split(below, cumsum(c(1, diff(below) != 1)))
    
    #get intervals
    temp <- do.call(rbind, lapply(1:length(seq_inds), function(y){
      data.frame(duration = intervals[seq_inds[[y]]], sequence = paste0(x, "-", y))
    }))
    return(temp)
  }
}))
minke_martin_2022$sequence <- as.numeric(factor(minke_martin_2022$sequence))

minke_data <- minke_martin_2022
save(minke_data, file = "data/processed/minke_data.RData")

# BOWHEAD WHALES ----------------------------------------------------------

bowhead_erbs_2021 <- read.csv("data/bowhead_erbs_2021.csv")
bowhead_erbs_2021 <- data.frame(duration = bowhead_erbs_2021$X..Delta.Time..s.., sequence = as.numeric(factor(bowhead_erbs_2021$stype)))

bowhead_data <- bowhead_erbs_2021
save(bowhead_data, file = "data/processed/bowhead_data.RData")

# RIGHT WHALES ------------------------------------------------------------

#manually went through and converted all variations of pattern duration column name to "Pattern Duration"
#added "PLACEHOLDER" in "Pattern Duration" to mark start of patterns in GS2-TP (old B) analyzed excels/M4 2013-14 Pattern B analyzed.xlsx
#WILL HAVE TO DO THIS MANUALLY FOR ALL OF THEM AND CHECK... too unreliable as is
#leaving out PU results in normal looking curve
#possible issue is that the "number gunshots" column for PU is not a reliable indicator of the start and end of songs

right_crance_2019 <- lapply(1:4, function(b){
  if(b == 1){directory <- "data/right_crance_2019/GS1-PF (old A) analyzed excels"}
  if(b == 2){directory <- "data/right_crance_2019/GS4-DG (old  D) analyzed excels"}
  if(b == 3){directory <- "data/right_crance_2019/GS3-PU (old C) analyzed excels"}
  if(b == 4){directory <- "data/right_crance_2019/GS2-TP (old B) analyzed excels"}
  
  files <- paste0(directory, "/", list.files(directory))
  if(b == 4){files <- files[-grep("end only", files)]}
  
  #for each file
  output <- lapply(1:length(files), function(x){
    #read in the file
    temp <- readxl::read_xlsx(files[x])
    
    #remove extra unnecessary terminal columns in file 13
    if(files[x] == "data/right_crance_2019/GS2-TP (old B) analyzed excels/M4 2013-14 Pattern B analyzed.xlsx"){temp <- temp[, -which(colnames(temp) %in% c("# terminals", "# terminal gunshots"))]}
    
    #remove final summary bit and extract ici columns
    temp <- temp[-which(is.na(temp[, 1])), ]
    
    # #extract main phrases dependending on the song type
    # if(b == 1){ici_cols <- which(colnames(temp) == "ICI gunshots")}
    # if(b == 2){ici_cols <- which(colnames(temp) %in% c("ICI doubles", "ICI rest"))}
    # if(b == 3){ici_cols <- which(colnames(temp) == "ICI gunshots")}
    # if(b == 4){ici_cols <- which(colnames(temp) == "ICI rest main")}
    
    ici_cols <- grep("ici", colnames(temp), ignore.case = TRUE)
    ici_cols <- ici_cols[which(ici_cols < grep("ipi", colnames(temp), ignore.case = TRUE))]
    
    #combined the ici values
    ici <- as.numeric(temp[, ici_cols[1]][[1]])
    if(length(ici_cols) > 1){
      for(i in 2:length(ici_cols)){
        next_ici <- as.numeric(temp[, ici_cols[i]][[1]])
        ici[which(!is.na(next_ici))] <- next_ici[which(!is.na(next_ici))]
        rm(next_ici)
      }
    }
    
    #get the start and end of each sequence (based on my manually identified bounds)
    seq_starts <- which(!is.na(temp[, grep("mason_bounds", colnames(temp))][[1]]))
    seq_ends <- c(c(seq_starts-1)[-1], nrow(temp))
    
    return(do.call(rbind, lapply(1:length(seq_starts), function(y){
      data.frame(duration = c(na.omit(ici[seq_starts[y]:seq_ends[y]])), sequence = paste0(gsub(".*/", "", files[x]), " ", y))
    })))
    
    # dur_cols <- c(grep("duration", colnames(temp), ignore.case = TRUE), grep("delta time", colnames(temp), ignore.case = TRUE))
    # dur_cols <- dur_cols[-which(dur_cols == grep("pattern", colnames(temp), ignore.case = TRUE))]
    # 
    # #combine the duration values
    # dur <- as.numeric(temp[, dur_cols[1]][[1]])
    # if(length(dur_cols) > 1){for(i in 2:length(dur_cols)){
    #   next_dur <- as.numeric(temp[, dur_cols[i]][[1]])
    #   dur[which(!is.na(next_dur))] <- next_dur[which(!is.na(next_dur))]
    #   rm(next_dur)
    # }}
    # 
    # #get the start and end of each sequence (based on the number of gunshots column)
    # seq_starts <- which(!is.na(temp[, grep("mason_bounds", colnames(temp))][[1]]))
    # seq_ends <- c(c(seq_starts-1)[-1], nrow(temp))
    # 
    # return(do.call(rbind, lapply(1:length(seq_starts), function(y){
    #                       data.frame(duration = c(na.omit(dur[seq_starts[y]:seq_ends[y]])), sequence = paste0(gsub(".*/", "", files[x]), " ", y))
    #                   })))
  })
  
  output <- do.call(rbind, output)
  output$sequence <- as.numeric(factor(output$sequence))
  return(output)
})

right_crance_2019[[1]]$sequence <- paste0("1-", right_crance_2019[[1]]$sequence)
right_crance_2019[[2]]$sequence <- paste0("2-", right_crance_2019[[2]]$sequence)
right_crance_2019[[3]]$sequence <- paste0("3-", right_crance_2019[[3]]$sequence)
right_crance_2019[[4]]$sequence <- paste0("4-", right_crance_2019[[4]]$sequence)
right_crance_2019 <- do.call(rbind, right_crance_2019)

right_data <- right_crance_2019
save(right_data, file = "data/processed/right_data.RData")

# NARROW-RIDGED FINLESS PORPOISES -----------------------------------------

#pulse packet parsing (lower level, makes more sense)
narrow_terada_2022 <- readxl::read_xlsx("data/narrow_terada_2022/narrow_terada_2022.xlsx")
to_iter <- which(!is.na(narrow_terada_2022$`Number of pulse`) & narrow_terada_2022$`Number of pulse` != 0)
narrow_terada_2022 <- do.call(rbind, lapply(1:length(to_iter), function(x){
  temp_inds <- (to_iter[x]-(narrow_terada_2022$`Number of pulse`[to_iter[x]]-1)):to_iter[x]
  data.frame(duration = narrow_terada_2022$IPI[temp_inds], sequence = x)
}))
narrow_terada_2022 <- narrow_terada_2022[-which(is.na(narrow_terada_2022$duration)), ]

narrow_data <- narrow_terada_2022
save(narrow_data, file = "data/processed/narrow_data.RData")

# BOTTLENOSE DOLPHINS -----------------------------------------------------

# bottlenose_moore_2020 <- read.csv("data/bottlenose_moore_2020.csv")
# bottlenose_moore_2020 <- bottlenose_moore_2020[which(bottlenose_moore_2020$Label == "IPI"), ] #only include solo males
# bottlenose_moore_2020 <- data.frame(duration = bottlenose_moore_2020$Duration..s., sequence = bottlenose_moore_2020$train)

# bottlenose_data <- bottlenose_moore_2020
# save(bottlenose_data, file = "data/processed/bottlenose_data.RData")

bottlenose_stepanov_2023 <- readxl::read_xlsx("data/bottlenose_stepanov_2023.xlsx")
bottlenose_data <- data.frame(duration = bottlenose_stepanov_2023$DURATION, sequence = bottlenose_stepanov_2023$SEQNO)
bottlenose_data$duration[which(bottlenose_data$duration == 0)] <- median(bottlenose_data$duration)
save(bottlenose_data, file = "data/processed/bottlenose_data.RData")

# RISSO'S DOLPHINS --------------------------------------------------------

rissos_arranz_2016 <- readxl::read_xlsx("data/rissos_arranz_2016.xlsx", sheet = 2)
pulsed_inds <- which(rissos_arranz_2016$`Inter-click interval (s)` < 0.016) #threshold for pulsed sounds in https://journals.biologists.com/jeb/article/219/18/2898/15466/Discrimination-of-fast-click-series-produced-by
split_inds <- split(pulsed_inds, cumsum(c(1, diff(pulsed_inds) != 1))) #https://stackoverflow.com/questions/24837401/find-consecutive-values-in-vector-in-r
rissos_arranz_2016 <- do.call(rbind, lapply(1:length(split_inds), function(x){
  data.frame(duration = rissos_arranz_2016$`Inter-click interval (s)`[split_inds[[x]]], sequence = x)
}))

#replace small number of zeros with mean
rissos_arranz_2016$duration[which(rissos_arranz_2016$duration == 0)] <- mean(rissos_arranz_2016$duration[-which(rissos_arranz_2016$duration == 0)])

rissos_data <- rissos_arranz_2016
save(rissos_data, file = "data/processed/rissos_data.RData")

# SEI WHALES --------------------------------------------------------------

sei_macklin_2024 <- readxl::read_xlsx("data/sei_macklin_2024.xlsx")

#iterate through and get sequence ids
sequence <- c()
counter <- 1
for(i in 1:nrow(sei_macklin_2024)){
  if(sei_macklin_2024$NoDownsweep[i] == "x"){
    sequence <- c(sequence, NA)
  } else{
    if(as.numeric(sei_macklin_2024$NoDownsweep[i]) == 1){
      counter <- counter + 1
      sequence <- c(sequence, counter)
    }
    if(as.numeric(sei_macklin_2024$NoDownsweep[i]) > 1){
      sequence <- c(sequence, counter)
    }
  }
}

unique_sequences <- na.omit(unique(sequence))
sei_macklin_2024 <- do.call(rbind, lapply(unique_sequences, function(x){
  data.frame(duration = na.omit(sei_macklin_2024$Duration95[which(sequence == x)]), sequence = x)
}))
# sei_macklin_2024 <- do.call(rbind, lapply(unique_sequences, function(x){
#   duration <- na.omit(sei_macklin_2024$Intracall.Spacing[which(sequence == x)])
#   if(length(duration) > 1){
#     data.frame(duration = duration, sequence = x)
#   }
# }))

sei_data <- sei_macklin_2024
save(sei_data, file = "data/processed/sei_data.RData")

# HEAVISIDE'S DOLPHINS ----------------------------------------------------

files <- list.files("data/heavisides_martin_2018/NBHF and BB bps")
files <- files[grep(".mat", files)]

heavisides_data <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/heavisides_martin_2018/NBHF and BB bps/", files[x]))$clicks[, 1]), sequence = x)
}))
save(heavisides_data, file = "data/processed/heavisides_data.RData")

# COMMERSON'S DOLPHINS ----------------------------------------------------

files <- list.files("data/commersons_martin_2021/NBHF bps")
files <- files[grep(".mat", files)]

commersons_nbhf <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/commersons_martin_2021/NBHF bps/", files[x]))$clicks[, 1]), sequence = x)
}))

files <- list.files("data/commersons_martin_2021/BB bps")
files <- files[grep(".mat", files)]

commersons_bb <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/commersons_martin_2021/BB bps/", files[x]))$clicks[, 1]), sequence = max(commersons_nbhf$sequence) + x)
}))

commersons_data <- rbind(commersons_nbhf, commersons_bb)
save(commersons_data, file = "data/processed/commersons_data.RData")

# HECTOR'S DOLPHINS -------------------------------------------------------

files <- list.files("data/hectors_nielsen_2024/NBHF bps")
files <- files[grep(".mat", files)]

hectors_nbhf <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/hectors_nielsen_2024/NBHF bps/", files[x]))$clicks[, 1]), sequence = x)
}))

files <- list.files("data/hectors_nielsen_2024/BB bps")
files <- files[grep(".mat", files)]

hectors_bb <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/hectors_nielsen_2024/BB bps/", files[x]))$clicks[, 1]), sequence = max(hectors_nbhf$sequence) + x)
}))

hectors_data <- rbind(hectors_nbhf, hectors_bb)
save(hectors_data, file = "data/processed/hectors_data.RData")

# PEALE'S DOLPHINS --------------------------------------------------------

files <- list.files("data/peales_martin_2024/NBHF bps")
files <- files[grep(".mat", files)]

peales_nbhf <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/peales_martin_2024/NBHF bps/", files[x]))$clicks[, 1]), sequence = x)
}))

files <- list.files("data/peales_martin_2024/BB bps")
files <- files[grep(".mat", files)]

peales_bb <- do.call(rbind, lapply(1:length(files), function(x){
  data.frame(duration = diff(R.matlab::readMat(paste0("data/peales_martin_2024/BB bps/", files[x]))$clicks[, 1]), sequence = max(peales_nbhf$sequence) + x)
}))

peales_data <- rbind(peales_nbhf, peales_bb)
save(peales_data, file = "data/processed/peales_data.RData")
