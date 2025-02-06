#written by Dan, Deus and Chikondi
#===============================================================================

#===============================================================================
#dataset for serotype-specific carriage
#===============================================================================

#load trackcare dataset and transform to long format
t0 <-
  rio::import(here::here("data", "trackcare_kids_serotyping.csv")) %>%
  dplyr::mutate(date = as.Date(collection_date, '%d/%m/%Y')) %>%
  dplyr::select(-starts_with('binary')) %>%
  dplyr::select(sampleid, date, starts_with('st'), pneumo_pos,piab_final) %>%
  melt(., id.vars=c('sampleid', 'date','pneumo_pos','piab_final'))  %>%
  dplyr::rename(st=variable) %>%
  dplyr::mutate(st_pos=if_else(value<45,1,0)) %>%
  dplyr::filter(!is.na(pneumo_pos)) %>%
  dplyr::arrange(st, pneumo_pos) %>%
  dplyr::mutate(st=toupper(st),
                st=gsub('ST','',st),
                st=if_else(st %in% c('23FAB'), '23F', st))

#compute the prevalence of each serotype
t1 <- 
  t0 %>%
  dplyr::group_by(st, pneumo_pos) %>%
  dplyr::summarize(N_samps = n(), N_pos = sum(st_pos), pct = mean(st_pos)) %>%
  dplyr::ungroup() %>%
  dcast( ., st ~ pneumo_pos, value.var = 'pct')

#download file and store in data folder
download.file('https://github.com/DanWeinberger/ABCs_pneumococcal_data/raw/refs/heads/master/Data/ABCs_st_1998_2021.rds', here::here("data", "abc_temp.RDS"))

#load ABCs dataset from data folder and compute total IPD for kids <5y
a1.kids <-
  rio::import(here::here("data", "abc_temp.RDS")) %>%
  dplyr::rename(agec = "Age.Group..years.",
                year = Year,
                st = IPD.Serotype,
                N_IPD = Frequency.Count) %>%
  dplyr::filter(agec %in% c("Age <2","Age 2-4") & year %in% c(2018, 2019, 2020, 2021)) %>%
  dplyr::mutate(st = if_else(st == '16', '16F', st),
                st = if_else(st %in% c('15B','15C', '15B/C'), '15BC', st)) %>%
  dplyr::group_by(st) %>%
  dplyr::summarize(N_IPD_kid = sum(N_IPD)) %>%
  dplyr::ungroup()

#load ABCs dataset from data folder and compute total IPD for older individuals >=5y
a1 <-
  rio::import(here::here("data", "abc_temp.RDS")) %>%
  dplyr::rename(agec = "Age.Group..years.",
                year = Year,
                st = IPD.Serotype,
                N_IPD = Frequency.Count) %>%
  dplyr::filter(!(agec %in% c("Age <2","Age 2-4")) & year %in% c(2018,2019,2020,2021)) %>%
  dplyr::mutate( st = if_else(st=='16','16F', st),
                 st = if_else(st %in% c('15B','15C', '15B/C'), '15BC', st)) %>%
  dplyr::group_by(st) %>%
  dplyr::summarize(N_IPD_older = sum(N_IPD)) %>%
  dplyr::ungroup() %>%
  dplyr::full_join(a1.kids, by ='st') %>%
  dplyr::mutate(N_IPD_older = if_else(is.na(N_IPD_older), 0, N_IPD_older),
                N_IPD_kid = if_else(is.na(N_IPD_kid), 0, N_IPD_kid))

#plot correlation between kids and older individuals
a1 %>%
  ggplot(aes(x = log(N_IPD_older+1), y = log(N_IPD_kid+1))) +
  geom_point() +
  theme_classic()

#augment the kids IPD data with adults where no observations
log.pred.kid.ipd <- predict(lm(log(N_IPD_kid+1) ~ log(N_IPD_older+1), data = a1))
a1 <- a1 %>% dplyr::mutate(N_IPD_kid = if_else(N_IPD_kid == 0, exp(log.pred.kid.ipd), N_IPD_kid))

#invasiveness data (from Navajo paper)
inv1 <- 
  read.csv('https://raw.githubusercontent.com/weinbergerlab/Invasiveness_Navajo/main/Results/mcmc_invasive_single_stage.csv') %>%
  dplyr::select('st','log.inv.age1') %>%
  dplyr::mutate(st = if_else(st %in% c('15B','15C', '15B/C'), '15BC', st))

#export invasiveness dataset to datra folder
rio::export(inv1, here::here("data", "invasiveness.csv"))

#reconcile inferred prevalence from ABCs IPD data with trackcare data and compute probability of false positive
b1 <- 
  base::merge(a1, inv1, by.x ='st', by.y ='st', all = T) %>%
  dplyr::mutate(log.inv.age1 = if_else(is.na(log.inv.age1), mean(log.inv.age1, na.rm = T), log.inv.age1), #if missing, assign the mean invassiveness
                N_IPD_kid = if_else(is.na(N_IPD_kid), 0.5, N_IPD_kid),
                invasiveness = exp(log.inv.age1),
                carr_est = N_IPD_kid/invasiveness/10000, #estimate prevalence (C)
                carr_prop = carr_est/sum(carr_est),  #pre
                
                #change serotypes to match those in the pcr panel
                st = if_else(st %in% c('17F','17A'), '17', st),
                st = if_else(st %in% c('7F','7A'), '7FA', st),
                st = if_else(st %in% c('9V','9A'), '9VA', st),
                st = if_else(st %in% c('11A','11D','11E'), '11ADE', st),
                st = if_else(st %in% c('9L','9N'), '9LN', st),
                st = if_else(st %in% c('22F','22A'), '22FA', st),
                st = if_else(st %in% c('6A','6B','6C','6D'), '6ABCD', st),
                #st = if_else(st %in% c('23F','23A','23B'), '23FAB', st),
                #st = if_else(st %in% c('23F','23A','23B'), '23FAB', st),
                st = if_else(st %in% c('28F','28A'), '28FA', st),
                st = if_else(st %in% c('33F','33A','37'), '33FA_37', st)) %>%
  
  dplyr::group_by(st) %>%
  dplyr::summarize(pre_test_prob = sum(carr_prop)) %>%
  dplyr::full_join(t1, by = 'st') %>%
  dplyr::ungroup() %>%
  dplyr::rename(false_pos = '0') %>%
  dplyr::mutate(false_pos = if_else(is.na(false_pos), 0, false_pos)) %>%
  dplyr::select(st, pre_test_prob, false_pos) %>%
  dplyr::mutate(sensitivity = 1,
                probability_pos = (sensitivity*pre_test_prob)/(sensitivity*pre_test_prob + (1-pre_test_prob)*false_pos)) 

#compute adjusted prevalence of the piab positive prevalence
t2 <- 
  t1 %>%
  dplyr::rename(raw_prev_piab_pos ='1') %>%
  dplyr::left_join(b1, by='st') %>%
  dplyr::mutate(adjusted_prevalence = raw_prev_piab_pos *probability_pos) %>%
  dplyr::filter(!is.na(adjusted_prevalence)) %>%
  mutate(st = if_else(st =='17', 'st17',
                      if_else(st =='7FA', 'st7fa',
                              if_else(st =='15BC', 'st15bc',
                                      if_else(st =='10A', 'st10a',
                                              if_else(st =='20', 'st20',
                                                      if_else(st =='9VA', 'st9va',
                                                              if_else(st =='23B', 'st23b',
                                                                      if_else(st =='11ADE', 'st11ade',
                                                                              if_else(st =='19A', 'st19a',
                                                                                      if_else(st =='9LN', 'st9ln',
                                                                                              if_else(st =='16F', 'st16f',
                                                                                                      if_else(st =='21', 'st21',
                                                                                                              if_else(st =='22FA', 'st22fa',
                                                                                                                      if_else(st =='3', 'st3',
                                                                                                                              if_else(st =='23F', 'st23f',
                                                                                                                                      if_else(st =='23A', 'st23a',
                                                                                                                                              if_else(st =='19F', 'st19f',
                                                                                                                                                      if_else(st =='28FA', 'st28fa',
                                                                                                                                                              if_else(st =='33FA_37', 'st33fa_37', NA_character_)))))))))))))))))))) %>%
  dplyr::filter(!is.na(st))

#plot piab prevalence against adjusted prevalence
t2 %>%
  ggplot(aes(x = raw_prev_piab_pos, y = adjusted_prevalence, label = st, color=st))+
  geom_text() +
  theme(legend.position = "none")

#Markov model on serotype-specific estimates

#the scaling of invasiveness doesn't really matter, especially if trying to use it in a different population because one would have to use a scaling factor to match the incidence in that other population anyway.
#load dataset
trackcare <- 
  rio::import(here("data", "trackcare_kids_serotyping.csv")) %>%
  dplyr::select(-V1) %>%
  dplyr::mutate(cdate = dmy(collection_date)) %>%
  dplyr::filter(!is.na(piab_final)) %>% #remove any mission final spn status
  dplyr::select(cdate, pid, season, piab_final, pneumo_pos, binary_st17:binary_st33fa_37)

#rename binary columns to their serotype
for (col in 6:ncol(trackcare)){
  colnames(trackcare)[col] <-  sub(".*binary_", "", colnames(trackcare)[col])
}

#rename columns and their values
trackcare <-
  trackcare %>%
  dplyr::rename(seas = season, piab = piab_final, spn_o = pneumo_pos) %>%
  dplyr::mutate(seas = as.factor(if_else(seas == "season1", "first", "second")))

#wrangle follow-up variables and generate new spn variable based on individual serotype status than overall status 
trackcare <-
  trackcare %>%
  dplyr::arrange(pid, cdate) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(dys = as.integer(cumsum(as.integer(cdate - lag(cdate, default = first(cdate)))))) %>%
  dplyr::ungroup() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(spn_n = sum(c(st17,st14,st7fa,st15bc,st10a,st20,st2,st9va,st23b,st11ade,st19a,st9ln,st16f,st21,st22fa,st34,st3,st6abcd,st23fab,st23a,st19f,st23f,st28fa,st33fa_37), na.rm = TRUE),
                spn_n = as.integer(if_else(spn_n >= 1, 2, 1)),
                spn_o = as.integer(if_else(spn_o ==1, 2, 1))) %>% 
  dplyr::mutate(obstx = as.integer(if_else(dys == 0, TRUE, FALSE)),
                obsty = as.integer(if_else(spn_o==1, TRUE, FALSE))) %>%
  dplyr::ungroup() %>%
  dplyr::select(cdate, dys, pid, seas, piab, spn_o, spn_n, obstx, obsty, everything())


#===============================================================================
#markov models for serotype-specific carriage adjusting for false positivity for all samples
#===============================================================================

#PERIOD 1 and 2 SAMPLES

#change dataset from wider to longer
trackcare_s <-
  trackcare %>%
  dplyr::select(everything(), -cdate, -seas, -piab, -spn_o, -spn_n, -obstx, -obsty) %>%
  dplyr::select(dys, pid, st17, st7fa, st15bc, st10a, st20, st9va, st23b, st11ade, st19a, st9ln, st16f, st21, st22fa , st3 , st23fab, st23a, st19f, st28fa, st33fa_37) %>%
  dplyr::rename("st23f" = "st23fab")

var <- c('st17','st7fa','st15bc','st10a','st20','st9va','st23b','st11ade','st19a','st9ln','st16f','st21','st22fa','st3','st23f','st23a','st19f','st28fa','st33fa_37')
trackcare_s[,var] <- sapply(trackcare_s[,var],function(x) ifelse(x == 1, 2L, 1L))

#Markov model including new carriage status, spn_n
#show transition frequency
trackcare_s <- arrange(trackcare_s, pid, dys)
statetable.msm(st28fa, pid, data = trackcare_s) #as an example of number of transitions 

#model list
modelfit <- list()
trackcare_x <-list()
t3 <- t2 %>% dplyr::select(st, probability_pos) %>% dplyr::rename("prob_pos" = "probability_pos")

for (st_name in var) {
  
  trackcare_x[[st_name]] <- 
    trackcare_s %>% dplyr::select(dys, pid, st_name) %>%
    mutate(st = if_else(st_name %in% c(t3 %>% rownames(.)[1]), st_name, NA_character_)) %>%
    dplyr::left_join(t3) %>%
    dplyr::select(everything(), -st)
  
  colnames(trackcare_x[[st_name]]) <- c('dys', 'pid', 'nstate', 'prob_pos')
  
  #initiate transition intensity matrix Q
  Qmatrix <- rbind(c(0.02, 0.02), c(0.02, 0.02))
  rownames(Qmatrix) <- c("clear", "carry")
  colnames(Qmatrix) <- c("clear", "carry")
  Qmatrix
  
  #run the Markov model
  modelfit[[st_name]] <- msm::msm(nstate ~ dys, 
                                  subject = pid, 
                                  data = trackcare_x[[st_name]],
                                  qmatrix = Qmatrix,
                                  opt.method = "bobyqa", control=list(maxfun=10000000))
}

modelrate <- list()
for (st_name in var) {
  modelrate[[st_name]] = qmatrix.msm(modelfit[[st_name]], ci = "normal", cl = 0.95)
}

#create a data frame of acquisitions rate and duration of carriage
modelDF1 <- data.frame(st = as.character(), acq = as.numeric(), acqL = as.numeric(), acqU = as.numeric(),  acq2 = as.numeric(), acqL2 = as.numeric(), acqU2 = as.numeric())
modelDF2 <- data.frame(st = as.character(), cle = as.numeric(), cleL = as.numeric(), cleU = as.numeric())

t4 <- t3 %>% base::split(list(.$st))

#acquisition rate
for (st_name in var) {
  modelDF1[st_name,] <- c(st_name, 
                          modelrate[[st_name]]$estimates[1,2]*t4[[st_name]][,2], 
                          modelrate[[st_name]]$L[1,2]*t4[[st_name]][,2],
                          modelrate[[st_name]]$U[1,2]*t4[[st_name]][,2],
                          modelrate[[st_name]]$estimates[1,2], 
                          modelrate[[st_name]]$L[1,2],
                          modelrate[[st_name]]$U[1,2])
}

#clearance rate
for (st_name in var) {
  modelDF2[st_name,] <- c(st_name, 
                          1/modelrate[[st_name]]$estimates[2,1],
                          1/modelrate[[st_name]]$U[2,1],
                          1/modelrate[[st_name]]$L[2,1])
}

#combine rates and probabilities
modelDF <-
  bind_cols(modelDF1, modelDF2 %>% select(everything(), -st)) %>%
  dplyr::filter(acq !=0) %>%
  dplyr::mutate(acq = round(as.numeric(acq), digits = 5),
                acqL = round(as.numeric(acqL), digits = 5),
                acqU = round(as.numeric(acqU), digits = 5),
                acq2 = round(as.numeric(acq2), digits = 5),
                acqL2 = round(as.numeric(acqL2), digits = 5),
                acqU2 = round(as.numeric(acqU2), digits = 5),
                cle = round(as.numeric(cle), digits = 5),
                cleL = round(as.numeric(cleL), digits = 5),
                cleU = round(as.numeric(cleU), digits = 5))

#rio::export(modelDF, file = here("output", "acqdur_overall.csv"))

#===============================================================================
#markov models for serotype-specific carriage adjusting for false positivity for spring 2021 samples
#===============================================================================

#PERIOD 1 SAMPLES

#change dataset from wider to longer
trackcare_s1 <-
  trackcare %>%
  dplyr::filter(seas == 'first') %>%
  dplyr::select(everything(), -cdate, -seas, -piab, -spn_o, -spn_n, -obstx, -obsty) %>%
  dplyr::select(dys, pid, st17, st7fa, st15bc, st10a, st20, st9va, st23b, st11ade, st19a, st9ln, st16f, st21, st22fa , st3 , st23fab, st23a, st19f, st28fa, st33fa_37) %>%
  dplyr::rename("st23f" = "st23fab")

var <- c('st17','st7fa','st15bc','st10a','st20','st9va','st23b','st11ade','st19a','st9ln','st16f','st21','st22fa','st3','st23f','st23a','st19f','st28fa','st33fa_37')
trackcare_s1[,var] <- sapply(trackcare_s1[,var],function(x) ifelse(x == 1, 2L, 1L))

#Markov model including new carriage status, spn_n
#show transition frequency
trackcare_s1 <- arrange(trackcare_s1, pid, dys)
statetable.msm(st28fa, pid, data = trackcare_s1) #as an example of number of transitions 

#model list
modelfit <- list()
trackcare_x <-list()
t3 <- t2 %>% dplyr::select(st, probability_pos) %>% dplyr::rename("prob_pos" = "probability_pos")

for (st_name in var) {
  
  trackcare_x[[st_name]] <- 
    trackcare_s1 %>% dplyr::select(dys, pid, st_name) %>%
    mutate(st = if_else(st_name %in% c(t3 %>% rownames(.)[1]), st_name, NA_character_)) %>%
    dplyr::left_join(t3) %>%
    dplyr::select(everything(), -st)
  
  colnames(trackcare_x[[st_name]]) <- c('dys', 'pid', 'nstate', 'prob_pos')
  
  #initiate transition intensity matrix Q
  Qmatrix <- rbind(c(0.02, 0.02), c(0.02, 0.02))
  rownames(Qmatrix) <- c("clear", "carry")
  colnames(Qmatrix) <- c("clear", "carry")
  Qmatrix
  
  #run the Markov model
  modelfit[[st_name]] <- msm::msm(nstate ~ dys, 
                                  subject = pid, 
                                  data = trackcare_x[[st_name]],
                                  qmatrix = Qmatrix,
                                  opt.method = "bobyqa", control=list(maxfun=10000000))
}

modelrate <- list()
for (st_name in var) {
  modelrate[[st_name]] = qmatrix.msm(modelfit[[st_name]], ci = "normal", cl = 0.95)
}

#create a data frame of acquisitions rate and duration of carriage
modelDF1a <- data.frame(st = as.character(), acq = as.numeric(), acqL = as.numeric(), acqU = as.numeric())
modelDF1b <- data.frame(st = as.character(), acq2 = as.numeric(), acqL2 = as.numeric(), acqU2 = as.numeric())
modelDF2 <- data.frame(st = as.character(), cle = as.numeric(), cleL = as.numeric(), cleU = as.numeric())

t4 <- t3 %>% base::split(list(.$st))

#acquisition rate
for (st_name in var) {
  modelDF1a[st_name,] <- c(st_name, 
                           modelrate[[st_name]]$estimates[1,2]*t4[[st_name]][,2], 
                           modelrate[[st_name]]$L[1,2]*t4[[st_name]][,2],
                           modelrate[[st_name]]$U[1,2]*t4[[st_name]][,2])
}

for (st_name in var) {
  modelDF1b[st_name,] <- c(st_name, 
                           modelrate[[st_name]]$estimates[1,2], 
                           modelrate[[st_name]]$L[1,2],
                           modelrate[[st_name]]$U[1,2])
}

#clearance rate
for (st_name in var) {
  modelDF2[st_name,] <- c(st_name, 
                          1/modelrate[[st_name]]$estimates[2,1],
                          1/modelrate[[st_name]]$U[2,1],
                          1/modelrate[[st_name]]$L[2,1])
}

#combine rates and probabilities
modelDF_per1 <-
  bind_cols(modelDF1a, modelDF1b %>% select(everything(), -st), modelDF2 %>% select(everything(), -st)) %>%
  dplyr::filter(acq !=0) %>%
  dplyr::mutate(acq = round(as.numeric(acq), digits = 5),
                acqL = round(as.numeric(acqL), digits = 5),
                acqU = round(as.numeric(acqU), digits = 5),
                acq2 = round(as.numeric(acq2), digits = 5),
                acqL2 = round(as.numeric(acqL2), digits = 5),
                acqU2 = round(as.numeric(acqU2), digits = 5),
                cle = round(as.numeric(cle), digits = 5),
                cleL = round(as.numeric(cleL), digits = 5),
                cleU = round(as.numeric(cleU), digits = 5))

#rio::export(modelDF_per1, file = here("output", "acqdur_firstseas.csv"))

#===============================================================================
#markov models for serotype-specific carriage adjusting for false positivity for winter/spring 2021/22 samples
#===============================================================================

#change dataset from wider to longer
trackcare_s2 <-
  trackcare %>%
  dplyr::filter(seas == 'second') %>%
  dplyr::select(everything(), -cdate, -seas, -piab, -spn_o, -spn_n, -obstx, -obsty) %>%
  dplyr::select(dys, pid, st17, st7fa, st15bc, st10a, st20, st9va, st23b, st11ade, st19a, st9ln, st16f, st21, st22fa , st3 , st23fab, st23a, st19f, st28fa, st33fa_37) %>%
  dplyr::rename("st23f" = "st23fab")

var <- c('st17','st7fa','st15bc','st10a','st20','st9va','st23b','st11ade','st19a','st9ln','st16f','st21','st22fa','st3','st23f','st23a','st19f','st28fa','st33fa_37')
trackcare_s2[,var] <- sapply(trackcare_s2[,var],function(x) ifelse(x == 1, 2L, 1L))

#Markov model including new carriage status, spn_n
#show transition frequency
trackcare_s2 <- arrange(trackcare_s2, pid, dys)
statetable.msm(st28fa, pid, data = trackcare_s2) #as an example of number of transitions 

#model list
modelfit <- list()
trackcare_x <-list()
t3 <- t2 %>% dplyr::select(st, probability_pos) %>% dplyr::rename("prob_pos" = "probability_pos")

for (st_name in var) {
  
  trackcare_x[[st_name]] <- 
    trackcare_s2 %>% dplyr::select(dys, pid, st_name) %>%
    mutate(st = if_else(st_name %in% c(t3 %>% rownames(.)[1]), st_name, NA_character_)) %>%
    dplyr::left_join(t3) %>%
    dplyr::select(everything(), -st)
  
  colnames(trackcare_x[[st_name]]) <- c('dys', 'pid', 'nstate', 'prob_pos')
  
  #initiate transition intensity matrix Q
  Qmatrix <- rbind(c(0.02, 0.02), c(0.02, 0.02))
  rownames(Qmatrix) <- c("clear", "carry")
  colnames(Qmatrix) <- c("clear", "carry")
  Qmatrix
  
  #run the Markov model
  modelfit[[st_name]] <- msm::msm(nstate ~ dys, 
                                  subject = pid, 
                                  data = trackcare_x[[st_name]],
                                  qmatrix = Qmatrix,
                                  opt.method = "bobyqa", control=list(maxfun=10000000))
}

modelrate <- list()
for (st_name in var) {
  modelrate[[st_name]] = qmatrix.msm(modelfit[[st_name]], ci = "normal", cl = 0.95)
}

#create a data frame of acquisitions rate and duration of carriage
modelDF1a <- data.frame(st = as.character(), acq = as.numeric(), acqL = as.numeric(), acqU = as.numeric())
modelDF1b <- data.frame(st = as.character(), acq2 = as.numeric(), acqL2 = as.numeric(), acqU2 = as.numeric())
modelDF2 <- data.frame(st = as.character(), cle = as.numeric(), cleL = as.numeric(), cleU = as.numeric())

t4 <- t3 %>% base::split(list(.$st))

#acquisition rate
for (st_name in var) {
  modelDF1a[st_name,] <- c(st_name, 
                           modelrate[[st_name]]$estimates[1,2]*t4[[st_name]][,2], 
                           modelrate[[st_name]]$L[1,2]*t4[[st_name]][,2],
                           modelrate[[st_name]]$U[1,2]*t4[[st_name]][,2])
}

for (st_name in var) {
  modelDF1b[st_name,] <- c(st_name, 
                           modelrate[[st_name]]$estimates[1,2], 
                           modelrate[[st_name]]$L[1,2],
                           modelrate[[st_name]]$U[1,2])
}

#clearance rate
for (st_name in var) {
  modelDF2[st_name,] <- c(st_name, 
                          1/modelrate[[st_name]]$estimates[2,1],
                          1/modelrate[[st_name]]$U[2,1],
                          1/modelrate[[st_name]]$L[2,1])
}

#combine rates and probabilities
modelDF_per2 <-
  bind_cols(modelDF1a, modelDF1b %>% select(everything(), -st), modelDF2 %>% select(everything(), -st)) %>%
  dplyr::filter(acq !=0) %>%
  dplyr::mutate(acq = round(as.numeric(acq), digits = 5),
                acqL = round(as.numeric(acqL), digits = 5),
                acqU = round(as.numeric(acqU), digits = 5),
                acq2 = round(as.numeric(acq2), digits = 5),
                acqL2 = round(as.numeric(acqL2), digits = 5),
                acqU2 = round(as.numeric(acqU2), digits = 5),
                cle = round(as.numeric(cle), digits = 5),
                cleL = round(as.numeric(cleL), digits = 5),
                cleU = round(as.numeric(cleU), digits = 5))

#rio::export(modelDF_per2, file = here("output", "acqdur_secondseas.csv"))

#===============================================================================
#===============================================================================

X <-
  bind_rows(
    t0 %>%
      dplyr::group_by(st, pneumo_pos) %>%
      dplyr::summarize(N_samps = n(), N_pos = sum(st_pos), pct = mean(st_pos), pctx=sum(st_pos)/N_samps) %>%
      dplyr::ungroup() %>%
      dplyr::filter(pneumo_pos==1, pct !=0) %>%
      dplyr::select(st, pct) %>%
      dplyr::mutate(seas = 'Overall') %>%
      dplyr::left_join(t2 %>% dplyr::select(st, adjusted_prevalence) %>% dplyr::mutate(st = stringr::str_sub(toupper(st), 3, nchar(st)))),
    
    t0 %>%
      dplyr::filter(date >=date('2021-02-04') & date <=date('2021-06-03')) %>%
      dplyr::group_by(st, pneumo_pos) %>%
      dplyr::summarize(N_samps = n(), N_pos = sum(st_pos), pct = mean(st_pos)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(pneumo_pos==1, pct !=0) %>%
      dplyr::select(st, pct) %>%
      dplyr::mutate(seas = 'Spring 2021') %>%
      dplyr::left_join(t2 %>% dplyr::select(st, adjusted_prevalence) %>% dplyr::mutate(st = stringr::str_sub(toupper(st), 3, nchar(st)))),
    
    t0 %>%
      dplyr::filter(date >=date('2021-11-10') & date <=date('2022-06-22')) %>%
      dplyr::group_by(st, pneumo_pos) %>%
      dplyr::summarize(N_samps = n(), N_pos = sum(st_pos), pct = mean(st_pos)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(pneumo_pos==1, pct !=0) %>%
      dplyr::select(st, pct) %>%
      dplyr::mutate(seas = 'Winter/Spring 2021/22') %>%
      dplyr::left_join(t2 %>% dplyr::select(st, adjusted_prevalence) %>% dplyr::mutate(st = stringr::str_sub(toupper(st), 3, nchar(st))))) %>%
  
  dplyr::mutate(seas = factor(seas, levels = c('Spring 2021','Winter/Spring 2021/22','Overall'))) %>%
  dplyr::filter(st != '6abcd', st != '34') %>%
  tidyr::complete(st, seas, fill = list(adjusted_prevalence = 0)) %>%
  dplyr::mutate(st = if_else(st=='6ABCD', '6A/B/C/D',
                             if_else(st=='7FA', '7F/A',
                                     if_else(st=='9LN', '9L/N',
                                             if_else(st=='9VA', '9V/A',
                                                     if_else(st=='11ADE', '11A/D/E',
                                                             if_else(st=='15BC', '15B/C',
                                                                     if_else(st=='22FA', '22F/A',
                                                                             if_else(st=='28FA', '28F/A',
                                                                                     if_else(st=='33FA_37', '33F/A/37', st)))))))))) %>%
  
  ggplot() + 
  geom_point(aes(x = adjusted_prevalence, y = reorder(st, adjusted_prevalence, decreasing = F), color = seas), stat = 'identity', shape = 4, size = 2.5, stroke = 2.5, position = position_dodge(width = 0.3)) + 
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(A)", x = "Adjusted prevalence", y = "Pneumococcal serotype") +
  scale_x_continuous(limit = c(0, 0.12), breaks = seq(0, 0.12, 0.02), labels = scales::percent_format(accuracy = 1)) + 
  guides(color = guide_legend(title = "")) +
  theme(legend.position = 'none', axis.text=element_text(size=16), legend.text = element_text(size=10)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#combined plot for serotype acquisition of the seasons above
seascol = c('Overall'='#619CFF', 'Spring 2021'='#F8766D', 'Winter/Spring 2021/22'='#00BA38')
Y <-
  bind_rows(
    data_frame(modelDF_per1) %>% dplyr::mutate(seas = 'Spring 2021') %>% dplyr::filter(st !='st19a'), 
    data_frame(modelDF_per2) %>% dplyr::mutate(seas = 'Winter/Spring 2021/22') %>% dplyr::filter(st !='st23a', st !='st23b'),
    data_frame(modelDF) %>% dplyr::mutate(seas = 'Overall')) %>%
  dplyr::filter(!is.na(acqL)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(st !='st9ln') %>%
  dplyr::mutate(st = stringr::str_sub(st, 3, nchar(st)),
                st = factor(st, levels = c('33fa_37','19a','15bc','23b','20','11ade','19f','10a','23f','23a','21','3','17','22fa','16f','9va','28fa','7fa')),
                seas = factor(seas, levels = c('Overall','Spring 2021','Winter/Spring 2021/22'))) %>%
  dplyr::mutate(st = toupper(st)) %>%
  dplyr::mutate(st = factor(if_else(st=='6ABCD', '6A/B/C/D',
                             if_else(st=='7FA', '7F/A',
                                     if_else(st=='9LN', '9L/N',
                                             if_else(st=='9VA', '9V/A',
                                                     if_else(st=='11ADE', '11A/D/E',
                                                             if_else(st=='15BC', '15B/C',
                                                                     if_else(st=='22FA', '22F/A',
                                                                             if_else(st=='28FA', '28F/A',
                                                                                     if_else(st=='33FA_37', '33F/A/37', st))))))))), levels = c('33F/A/37','19A','15B/C','23B','20','11A/D/E','19F','10A','23F','23A','21','3','17','22F/A','16F','9V/A','28F/A','7F/A'))) %>%
  ggplot() +
  geom_point(aes(x=st, y=log(acq), size = cle, color = seas),  shape=1, stroke = 2.5, position = position_jitter(width = 0,seed=1988), stat = "identity") +
  geom_errorbar(aes(x=st, y=log(acq),  ymin = log(acqL), ymax = log(acqU), color = seas), width = 0, size = 1.2, position = position_jitter(width = 0, seed=1988), stat = "identity") +
  
  geom_point(aes(x=st, y=log(acq2), size = cle), color = 'gray70', stroke = 2, shape=1, position = position_jitter(width = 0.2, seed=1988), stat = "identity") +
  geom_errorbar(aes(x=st, y=log(acq2),  ymin = log(acqL2), ymax = log(acqU2)), color = 'gray70', width = 0, size = 0.8, position = position_jitter(width = 0.2, seed=1988), stat = "identity") +
  scale_color_manual(values = seascol) +
  
  scale_y_continuous(breaks = c(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11), limits = c(-11, -2.5)) + 
  #coord_cartesian(ylim = c(-11, -2.5)) +
  theme_bw(base_size = 14, base_family = 'American Typewriter') + 
  labs(title = "(B)", x = "Pneumococcal serotype", y = "Daily carriage acquisition log_rate") + 
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(seas ~.) +
  theme(legend.text = element_text(size = 12)) +
  guides(color = guide_legend(title="Study period"), size = guide_legend(title="Carriage\nduration")) +
  theme(plot.margin=grid::unit(c(0.5,0,0,0), "mm")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(strip.text.y = element_text(size = 20), strip.background=element_rect(fill="gray90"))

#combined plots
ggsave(here::here("output", "ST_estimates.png"),
       plot = ((X | Y | plot_layout(ncol = 2, width = c(1,3)))),
       width = 20, height = 13, unit = "in", dpi = 300)

