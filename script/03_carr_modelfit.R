#written by Dan, Deus and Chikondi
#===============================================================================

#===============================================================================
#data manipulation for analysis
#===============================================================================

#load carriage data
tc_data <- 
  rio::import(here::here("data","kidsData.rds")) %>%
  dplyr::select(sample_id, collection_date, pid, daycare_location, household_size, race.y, sex.y, age_years, piab_bi, lyta_bi, piab2, ethnicity, season.x) %>%
  
  #data munging
  dplyr::rename("daycare"="daycare_location", "cdate"="collection_date", "hhsize"="household_size", "ctvalue"="piab2", "race"="race.y", "sex"="sex.y", "agey"="age_years", "ethn"="ethnicity", "seas"="season.x") %>%
  dplyr::mutate(cdate = date(cdate),
                seas = if_else(seas == "season1",1L,2L),
                pid = stringr::str_c(pid,seas)) %>%
  dplyr::arrange(pid, cdate) %>%
  group_by(pid) %>%
  dplyr::mutate(dys = cumsum(as.integer(cdate - lag(cdate, default = first(cdate))))) %>%
  ungroup() %>%
  
  dplyr::group_by(daycare) %>%
  dplyr::mutate(dycare = n_distinct(pid)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(seasx = if_else(year(cdate) == 2022 & date(cdate)>='2022-03-20', 1, seas)) %>% #sensitivity analysis
  
  dplyr::mutate(yearc = factor(lubridate::year(cdate)),
                dycare = as.factor(if_else(dycare <=10, "lowdens", 
                                           if_else(dycare >10, "highdens", NA_character_))),
                hhsize = as.factor(if_else(hhsize <=2, "2",
                                           if_else(hhsize >=3, "3+", NA_character_))),
                race = as.factor(if_else(race=="", NA_character_,
                                         if_else(race=="White", "White", "nonWhite"))),
                sex = as.factor(if_else(sex=="Male", "Male",
                                        if_else(sex=="Female", "Female", NA_character_))),
                agey = as.factor(if_else(agey <=2, "0-2y", 
                                         if_else(agey >2, "3+y", NA_character_))),
                piab_bi = as.integer(if_else(is.na(piab_bi), 9L, 
                                             if_else(piab_bi==1L, 2L, 1L))),
                lyta_bi = as.integer(if_else(is.na(lyta_bi), 9L, 
                                             if_else(lyta_bi==1L, 2L, 1L))),
                spnden = as.factor(if_else(piab_bi==1, "None",
                                           if_else(piab_bi==2 & ctvalue >=30, "Low", 
                                                   if_else(piab_bi==2 & ctvalue <30, "High", NA_character_)))),
                ethn = as.factor(if_else(ethn=="Hispanic or Latino", "Hispanic",
                                         if_else(ethn=="NOT Hispanic or Latino", "nonHispanic", NA_character_))),
                seas = as.factor(if_else(seas==1, "Period1",
                                         if_else(seas==2, "Period2", NA_character_))),
                seasx = as.factor(if_else(seasx==1, "Spring",
                                          if_else(seasx==2, "Winter", NA_character_)))) %>%
  dplyr::select(sample_id, cdate, yearc, pid, dys, piab_bi, lyta_bi, dycare, agey, sex, spnden, ethn, race, hhsize, seas, seasx) %>%
  dplyr::filter(piab_bi !=9)

#===============================================================================

tc_dataS <-
  rio::import(here::here("data", "trackcare_kids_serotyping.csv")) %>%
  dplyr::mutate(date = as.Date(collection_date, '%d/%m/%Y')) %>%
  dplyr::select(-starts_with('binary')) %>%
  dplyr::select(sampleid, date, starts_with('st'), pneumo_pos, piab_final) %>%
  tidyr::pivot_longer(cols = c(starts_with('st')),
                      names_to = c("st"),
                      values_to = "Ct_st") %>%
  dplyr::mutate(st=if_else(st %in% c('st23fab'), 'st23f', st)) %>%
  dplyr::rename('sample_id'='sampleid', 'dateS'='date', 'Ct_final'='piab_final') %>%
  dplyr::mutate(pneumo_pos = if_else(pneumo_pos==1, 2L,1L)) %>%

#join with dataset of 'probability of true Spn positive'   
  dplyr::left_join(t3) %>%
  dplyr::filter(st != 'st6abcd', st != 'st34', st != 'st2', st != 'st14') %>%
  dplyr::filter(st != 'st7fa', st != 'st28fa', st != 'st9va') %>% #st != 'st33fa_37'
  #dplyr::mutate(Ct_stx = Ct_st+(Ct_st*(1-prob_pos))) %>% # Ct_st*exp(1-prob_pos))

#join with covariate dataset
  dplyr::left_join(tc_data) %>%
  dplyr::group_by(sample_id, dys) %>%
  
#if st neg then neg, but if st pos on piab then look up the Ct values
  dplyr::mutate(dom_st = if_else(piab_bi==1, 'Uncol',
                                 if_else(piab_bi==2,  st[which.min(Ct_st)], NA_character_))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(dom_st)) %>%
  dplyr::distinct(sample_id, dom_st, .keep_all = TRUE) %>%
  dplyr::arrange(pid, dateS) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(dys = cumsum(as.integer(dateS - lag(dateS, default = first(dateS))))) %>%
  dplyr::ungroup() %>%
  
  dplyr::arrange(pid,dys) %>%
  dplyr::group_by(pid) %>%
  dplyr::mutate(visit_index = row_number(),
                lag_dom_st = lag(dom_st, 1),
                lag_dom_st_adjusted = if_else(visit_index == 1, 'first', lag_dom_st),
                person_index= cur_group_id(),
                state = NA_integer_) %>%
  dplyr::ungroup()

#function to switch between suceptible, infection state 1 and infection state2
col_switch_fun <- function(X) {
  ds1.sub <- tc_dataS %>% dplyr::filter(person_index == X)
  ds1.sub$state <- NA_integer_
  for (k in 1:nrow(ds1.sub)) {
    if (k == 1) {
      ds1.sub$state[k] <- if_else(ds1.sub$dom_st[k] == 'Uncol', 1, 2)
    } else {
      ds1.sub$state[k] <- if_else(ds1.sub$dom_st[k] == 'Uncol', 1,
                                  if_else(ds1.sub$dom_st[k] != 'Uncol' & ds1.sub$lag_dom_st[k] == 'Uncol',2,
                                          if_else(ds1.sub$dom_st[k] == ds1.sub$lag_dom_st[k] & ds1.sub$state[k-1]==2,2,
                                                  if_else(ds1.sub$dom_st[k] == ds1.sub$lag_dom_st[k] & ds1.sub$state[k-1]==3,3,
                                                          if_else(ds1.sub$dom_st[k] != ds1.sub$lag_dom_st[k] & ds1.sub$state[k-1]==2,3,
                                                                  if_else(ds1.sub$dom_st[k] != ds1.sub$lag_dom_st[k] & ds1.sub$state[k-1]==3,2,NA_integer_))))))
                                 
    }
  }
  return(ds1.sub)
}

tc_dataSX <- lapply(1:max(tc_dataS$person_index), function(X) {
  if (nrow(tc_dataS %>% dplyr::filter(person_index == X)) > 0) {
    col_switch_fun(X)
  } else {
    NULL
  }
  }) %>% 
  dplyr::bind_rows() %>%
  dplyr::select(everything(), -cdate, -lyta_bi, -visit_index, -lag_dom_st, -lag_dom_st_adjusted, -person_index) %>%
  dplyr::mutate(sswitch = factor(if_else((lag(state)==2 & state==3) | (lag(state)==3 & state==2), 'on', 'off'))) #phi parameter in the model diagram as covariate

#===============================================================================
#markov models for whole carriage and risk factors (main)
#===============================================================================

# show transition frequency
tc_fup <- arrange(tc_dataSX, pid, dys)
statetable.msm(state, pid, data = tc_fup)

#initiate transition intensity matrix Q
spn_Qmatrix <- rbind(c(0.05, 0.05, 0.05),
                     c(0.05, 0.05, 0.05),
                     c(0.05, 0.05, 0.05))

rownames(spn_Qmatrix) <- c('S', 'I1', 'I2')
colnames(spn_Qmatrix) <- c('S', 'I1', 'I2')
spn_Qmatrix

#run the Markov model, both seasons
spn_modelfit0 <- msm(state ~ dys, subject = pid, data = tc_fup,
                     qmatrix = spn_Qmatrix,
                     qconstraint = c(1,1,2,1,2,1),
                     covariates = list("1-2" = ~ agey + sex + ethn + seas + seasx + hhsize + sswitch,
                                       "2-1" = ~ agey + sex),
                     opt.method = "bobyqa", control = list(maxfun = 10000000))
Overall <- hazard.msm(spn_modelfit0, hazard.scale = 1, cl = 0.95)

#run the Markov model, season 1
spn_modelfit1 <- msm(state ~ dys, subject = pid, data = tc_fup %>% dplyr::filter(seas=="Period1"),
                     qmatrix = spn_Qmatrix,
                     qconstraint = c(1,1,2,1,2,1),
                     covariates = list("1-2" = ~ agey + sex + ethn + hhsize + sswitch, 
                                       "2-1" = ~ agey + sex),
                     opt.method = "bobyqa", control = list(maxfun = 10000000))
Period1 <- hazard.msm(spn_modelfit1, hazard.scale = 1, cl = 0.95)

#run the Markov model, season 2
spn_modelfit2 <- msm(state ~ dys, subject = pid, data = tc_fup %>% dplyr::filter(seas=="Period2"),
                     qmatrix = spn_Qmatrix,
                     qconstraint = c(1,1,2,1,2,1),
                     covariates = list("1-2" = ~ agey + sex + ethn + seasx + hhsize + sswitch, 
                                       "2-1" = ~ agey + sex ),
                     opt.method = "bobyqa", control = list(maxfun = 10000000))
Period2 <- hazard.msm(spn_modelfit2, hazard.scale = 1, cl = 0.95)

#===============================================================================
#ACQUSITION
#===============================================================================

#OVERALL
#create an empty acquisition table
OverallDS <- tibble(
  index = c(1:5),
  dstype = "Overall",
  cova = c(rep("agey", 1), rep("sex", 1), rep("ethn", 1), rep("seas", 1), rep("hhsize", 1)),
  label = c(rep("Age group \n (≤2y) vs 3-6y", 1), rep("Sex \n (female) vs male", 1), rep("Ethnicity \n (Hispanic) vs Non-Hispanic ", 1), rep("Study period \n (Spring 2021) vs Winter/Spring 2021/22", 1), rep("Household size \n (2 members) vs 3+ members", 1)),
  HR = c(rep(0, 5)),
  HRlci = c(rep(0, 5)),
  HRuci = c(rep(0, 5)))

#insert HR from acquisition values in the Overall fitted model
l = 1; m = 1; n = 1
for(i in c("agey3+y", "sexMale", "ethnnonHispanic", "seasPeriod2", "hhsize3+")){
  for(j in 1){
    OverallDS[l, 5] = Overall[[i]][1,j]
    l = l+1
  }
  for(j in 2){
    OverallDS[m, 6] = Overall[[i]][1,j]
    m = m+1
  }
  for(j in 3){
    OverallDS[n, 7] = Overall[[i]][1,j]
    n = n+1
  }
}

#PERIOD 1
#create an empty acquisition table
Period1DS <- tibble(
  index = c(1:4),
  dstype = "Spring 2021",
  cova = c(rep("agey", 1), rep("sex", 1), rep("ethn", 1), rep("hhsize", 1)),
  label = c(rep("Age group \n (≤2y) vs 3-6y", 1), rep("Sex \n (female) vs male", 1), rep("Ethnicity \n (Hispanic) vs Non-Hispanic ", 1), rep("Household size \n (2 members) vs 3+ members", 1)),
  HR = c(rep(0, 4)),
  HRlci = c(rep(0, 4)),
  HRuci = c(rep(0, 4)))

#insert HR from acquisition values in the Overall fitted model
l = 1; m = 1; n = 1
for(i in c("agey3+y", "sexMale", "ethnnonHispanic", "hhsize3+")){
  for(j in 1){
    Period1DS[l, 5] = Period1[[i]][1,j]
    l = l+1
  }
  for(j in 2){
    Period1DS[m, 6] = Period1[[i]][1,j]
    m = m+1
  }
  for(j in 3){
    Period1DS[n, 7] = Period1[[i]][1,j]
    n = n+1
  }
}

#PERIOD 2
#create an empty acquisition table
Period2DS <- tibble(
  index = c(1:4),
  dstype = "Winter/Spring 2021/22",
  cova = c(rep("agey", 1), rep("sex", 1), rep("ethn", 1), rep("hhsize", 1)),
  label = c(rep("Age group \n (≤2y) vs 3-6y", 1), rep("Sex \n (female) vs male", 1), rep("Ethnicity \n (Hispanic) vs Non-Hispanic ", 1), rep("Household size \n (2 members) vs 3+ members", 1)),
  HR = c(rep(0, 4)),
  HRlci = c(rep(0, 4)),
  HRuci = c(rep(0, 4)))

#insert HR from acquisition values in the Overall fitted model
l = 1; m = 1; n = 1
for(i in c("agey3+y", "sexMale", "ethnnonHispanic", "hhsize3+")){
  for(j in 1){
    Period2DS[l, 5] = Period2[[i]][1,j]
    l = l+1
  }
  for(j in 2){
    Period2DS[m, 6] = Period2[[i]][1,j]
    m = m+1
  }
  for(j in 3){
    Period2DS[n, 7] = Period2[[i]][1,j]
    n = n+1
  }
}

#===============================================================================
#CLEARANCE
#===============================================================================

#OVERALL
#create an empty clearance table
OverallDSx <- tibble(
  index = c(1:2),
  dstype = "Overall",
  cova = c(rep("agey", 1), rep("sex", 1)),
  label = c(rep("Age group \n (≤2y) vs 3-6y", 1), rep("Sex \n (female) vs male", 1)),
  HR = c(rep(0, 2)),
  HRlci = c(rep(0, 2)),
  HRuci = c(rep(0, 2)))

#insert HR from acquisition values in the Overall fitted model
l = 1; m = 1; n = 1
for(i in c("agey3+y", "sexMale")){
  for(j in 1){
    OverallDSx[l, 5] = Overall[[i]][3,j]
    l = l+1
  }
  for(j in 2){
    OverallDSx[m, 6] = Overall[[i]][3,j]
    m = m+1
  }
  for(j in 3){
    OverallDSx[n, 7] = Overall[[i]][3,j]
    n = n+1
  }
}

#PERIOD 1
#create an empty clearance table
Period1DSx <- tibble(
  index = c(1:2),
  dstype = "Spring 2021",
  cova = c(rep("agey", 1), rep("sex", 1)),
  label = c(rep("Age group \n (≤2y) vs 3-6y", 1), rep("Sex \n (female) vs male", 1)),
  HR = c(rep(0, 2)),
  HRlci = c(rep(0, 2)),
  HRuci = c(rep(0, 2)))

#insert HR from acquisition values in the Overall fitted model
l = 1; m = 1; n = 1
for(i in c("agey3+y", "sexMale")){
  for(j in 1){
    Period1DSx[l, 5] = Period1[[i]][3,j]
    l = l+1
  }
  for(j in 2){
    Period1DSx[m, 6] = Period1[[i]][3,j]
    m = m+1
  }
  for(j in 3){
    Period1DSx[n, 7] = Period1[[i]][3,j]
    n = n+1
  }
}

#PERIOD 2
#create an empty clearance table
Period2DSx <- tibble(
  index = c(1:2),
  dstype = "Winter/Spring 2021/22",
  cova = c(rep("agey", 1), rep("sex", 1)),
  label = c(rep("Age group \n (≤2y) vs 3-6y", 1), rep("Sex \n (female) vs male", 1)),
  HR = c(rep(0, 2)),
  HRlci = c(rep(0, 2)),
  HRuci = c(rep(0, 2)))

#insert HR from acquisition values in the Overall fitted model
l = 1; m = 1; n = 1
for(i in c("agey3+y", "sexMale")){
  for(j in 1){
    Period2DSx[l, 5] = Period2[[i]][3,j]
    l = l+1
  }
  for(j in 2){
    Period2DSx[m, 6] = Period2[[i]][3,j]
    m = m+1
  }
  for(j in 3){
    Period2DSx[n, 7] = Period2[[i]][3,j]
    n = n+1
  }
}

#===============================================================================

#plot the acquisition and clearance hazard ratios
A <- 
  bind_rows(OverallDS, Period1DS, Period2DS) %>% dplyr::mutate(carrtype = "Carriage acquisition") %>%
  dplyr::mutate(dstype = factor(dstype, levels = c('Spring 2021','Winter/Spring 2021/22','Overall'))) %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = dstype)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0, size = 1, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 1, size = 3,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7,8)) + 
  coord_cartesian(xlim = c(c(0, 6))) +
  labs(title = "(A) Acquisition", x = "Hazard Ratio (HR)", y = "") + 
  facet_grid(.~dstype, scales = 'free_y') +
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) + 
  theme(strip.text.x = element_text(size = 16), strip.background = element_rect(fill="gray90")) +
  theme(legend.text = element_text(size = 12), legend.position = 'none') +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

B <- 
  bind_rows(OverallDSx, Period1DSx, Period2DSx) %>% dplyr::mutate(carrtype = "Carriage duration") %>%
  dplyr::mutate(dstype = factor(dstype, levels = c('Spring 2021','Winter/Spring 2021/22','Overall'))) %>%
  ggplot(aes(y = reorder(label, index), x = HR, color = dstype)) +
  geom_errorbarh(aes(xmin = HRlci, xmax = HRuci), height = 0, size = 1, position = position_dodge2(width = 0.5)) +
  geom_point(shape = 1, size = 3,  stroke = 2, position = position_dodge2(width = 0.5)) +  
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 0.6, alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7,8)) + 
  coord_cartesian(xlim = c(c(0, 6))) +
  labs(title = "(B) Clearance", x = "Hazard Ratio (HR)", y = "") + 
  facet_grid(. ~ dstype, scales = 'free_y') +
  theme_bw(base_size = 14, base_family = 'American Typewriter') +
  guides(color = guide_legend(title = "")) +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) + 
  theme(strip.text.x = element_text(size = 16), strip.background = element_rect(fill="gray90")) +
  theme(legend.text = element_text(size = 12), legend.position = 'bottom') +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#combined plots
ggsave(here::here("output", "carr_estimates.png"),
       plot = ((A + B) + plot_layout(nrow = 2, ncol = 1, heights = c(3,1))),
       width = 16, height = 12, unit = "in", dpi = 300)

