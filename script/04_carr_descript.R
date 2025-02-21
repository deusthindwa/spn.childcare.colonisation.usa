#written by Dan, Deus and Chikondi
#===============================================================================

#===============================================================================
#aggregated prevalence
#===============================================================================

#aggregated prevalence of piaB carriage

bind_cols(
  tc_data %>%
    dplyr::mutate(piab_bi = if_else(piab_bi==1, 0, if_else(piab_bi==2, 1, NA_integer_))) %>% 
    dplyr::group_by(seas) %>%
    tally() %>%
    dplyr::ungroup() %>%
    dplyr::rename('N'='n'),
  
  tc_data %>%
    dplyr::mutate(piab_bi = if_else(piab_bi==1, 0, if_else(piab_bi==2, 1, NA_integer_))) %>% 
    dplyr::group_by(seas) %>%
    dplyr::tally(piab_bi) %>%
    dplyr::ungroup() %>%
    dplyr::select(n)) %>%
  
  dplyr::rowwise() %>% 
  dplyr::mutate(p = n/N,
                p_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], 
                p_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  dplyr::ungroup()


#aggregated prevalence lytA of carriage
bind_cols(
  tc_data %>%
    dplyr::mutate(lyta_bi = if_else(lyta_bi==1, 0, if_else(lyta_bi==2, 1, NA_integer_))) %>% 
    dplyr::group_by(seas) %>%
    tally() %>%
    dplyr::ungroup() %>%
    dplyr::rename('N'='n'),
  
  tc_data %>%
    dplyr::mutate(lyta_bi = if_else(lyta_bi==1, 0, if_else(lyta_bi==2, 1, NA_integer_))) %>% 
    dplyr::group_by(seas) %>%
    dplyr::tally(lyta_bi) %>%
    dplyr::ungroup() %>%
    dplyr::select(n)) %>%
  
  dplyr::rowwise() %>% 
  dplyr::mutate(p = n/N,
                p_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], 
                p_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
  dplyr::ungroup()

#===============================================================================
#characteristics of piaB carriage prevalence (figure 1)
#===============================================================================

#manipulating carriage dataset for plotting
tc_datap <-
  tc_data %>%
  dplyr::select(cdate, piab_bi, lyta_bi, seas, spnden) %>%
  dplyr::mutate(seas = if_else(seas=="Period1", "Spring 2021", "Winter/Spring 2021/22")) %>%
  dplyr::mutate(piab_bi = if_else(piab_bi==1, 0, if_else(piab_bi==2, 1, NA_integer_)), 
                lyta_bi = if_else(lyta_bi==1, 0, if_else(lyta_bi==2, 1, NA_integer_)),
                yrmo = str_sub(cdate, 1,7)) %>%
  dplyr::select(yrmo, piab_bi, lyta_bi, seas, spnden)

#manipulating stringency dataset for plotting
stringency <- 
  rio::import(here::here("data","stringency.xlsx")) %>%
  dplyr::mutate(yrmo = str_sub(cdate, 1,7), strin = strin/100) %>%
  dplyr::group_by(yrmo, period) %>%
  dplyr::summarise(strin = mean(strin)) %>%
  dplyr::ungroup()

tc_data %>%
  dplyr::select(cdate, piab_bi, lyta_bi, seas, spnden) %>%
  dplyr::mutate(seas = if_else(seas=="Period1", "Spring 2021", "Winter/Spring 2021/22")) %>%
  dplyr::mutate(piab_bi = if_else(piab_bi==1, 0, if_else(piab_bi==2, 1, NA_integer_)), 
                lyta_bi = if_else(lyta_bi==1, 0, if_else(lyta_bi==2, 1, NA_integer_)),
                yrmo = str_sub(cdate, 1,7)) %>%
  dplyr::select(yrmo, piab_bi, lyta_bi, seas, spnden)


#load schematic model diagrams
fig1a <- rasterGrob(readPNG(here("data","fig1a.png")), interpolate=TRUE)
A <- 
  ggplot() + 
  annotation_custom(fig1a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(A)", x = "Overall carriage model structure", y = "") +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))
  

fig1b <- rasterGrob(readPNG(here("data","fig1b.png")), interpolate=TRUE)
B <- 
  ggplot() + 
  annotation_custom(fig1b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(B)", x = "Individual serotype model structure", y = "") +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

C <-
  bind_rows(
    tc_datap %>% dplyr::group_by(seas, spnden) %>% dplyr::tally(piab_bi) %>% dplyr::filter(spnden !='None') %>% 
      dplyr::left_join(tc_datap %>% dplyr::group_by(seas) %>% dplyr::tally() %>% dplyr::ungroup() %>% rename('N'='n') %>% dplyr::select(seas, N)) %>% 
      dplyr::mutate(p = n/N, meth = "piaB") %>%
      rowwise() %>% 
      dplyr::mutate(p_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], 
                    p_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
      dplyr::ungroup(),
    
    tc_datap %>% dplyr::group_by(seas, spnden) %>% dplyr::tally(lyta_bi) %>% dplyr::filter(spnden !='None') %>% 
      dplyr::left_join(tc_datap %>% dplyr::group_by(seas) %>% dplyr::tally() %>% dplyr::ungroup() %>% rename('N'='n') %>% dplyr::select(seas, N)) %>% 
      dplyr::mutate(p = n/N, meth = "lytA") %>%
      rowwise() %>% 
      dplyr::mutate(p_lci = exactci(n, N, 0.95) [["conf.int"]][[1]], 
                    p_uci = exactci(n, N, 0.95) [["conf.int"]][[2]]) %>%
      dplyr::ungroup()) %>%
  dplyr::mutate(gp = str_c(seas,meth)) %>%
  
  ggplot() +
  geom_point(aes(x = spnden, y = p, color = seas, shape = meth, group=gp), size = 3, stroke = 1.5, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(x = spnden, y = p, ymin = p_lci, ymax = p_uci, color = seas, group = gp), width = 0, size = 1, position = position_dodge(width = 0.6)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_shape_manual(values = c(4,1)) +
  scale_color_manual(values = c('#F8766D', '#00BA38')) +
  scale_y_continuous(breaks=c(0.10, 0.20, 0.30, 0.40, 0.50), limits = c(0.05, 0.45), labels = scales::percent_format(accuracy = 1)) + 
  labs(title = "(C)", x = "Carriage density", y = "Unadjusted prevalence") + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  guides(shape = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 10), legend.position = "none", legend.title = element_text(size = 10)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#plot the data
D <-
  bind_rows(
    bind_cols(tc_datap %>% dplyr::group_by(yrmo, seas) %>% dplyr::tally(piab_bi) %>% dplyr::ungroup(),
              tc_datap %>% dplyr::group_by(yrmo, seas) %>% dplyr::tally() %>% dplyr::ungroup() %>% rename('N'='n') %>% dplyr::select(N)) %>%
      dplyr::mutate(p = n/N, meth='piaB'),
    
    bind_cols(tc_datap %>% dplyr::group_by(yrmo, seas) %>% dplyr::tally(lyta_bi),
              tc_datap %>% dplyr::group_by(yrmo, seas) %>% dplyr::tally() %>% dplyr::ungroup() %>% rename('N'='n') %>% dplyr::select(N)) %>%
      dplyr::mutate(p = n/N, meth='lytA')) %>%
  
  rowwise() %>%
  dplyr::mutate(p_lci = exactci(n, N, 0.95) [["conf.int"]][[1]],
                p_uci = exactci(n, N, 0.95) [["conf.int"]][[2]],
                cutline = str_c(meth,seas)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(stringency) %>%
  
  ggplot() +
  geom_point(aes(x = yrmo, y = p, color = seas, shape=meth, group=cutline), size =3, stroke = 2, , position = position_dodge(width = 0.2)) +
  geom_line(aes(x = yrmo, y = p, color = seas, group = cutline), size = 1, position = position_dodge(width = 0.2)) +
  geom_line(aes(x = yrmo, y = strin, group=1), linetype = "dashed", size = 1) +
  geom_errorbar(aes(x = yrmo, y = p, ymin = p_lci, ymax = p_uci, color = seas, group=meth), width = 0, size = 1, position = position_dodge(width = 0.2)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_shape_manual(values = c(4,1)) +
  scale_color_manual(values = c('#F8766D', '#00BA38')) +
  scale_y_continuous(limit = c(0, 0.95), breaks = seq(0, 0.95, 0.10), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "(D)", x = "Sampling month", y = "Unadjusted prevalence") +
  theme(plot.title = element_text(size = 20), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5), axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size = 12), legend.title = element_text(size = 12)) +
  guides(shape = guide_legend(title = "qPCR gene:"), color = guide_legend(title = "Study period:")) +
  theme(legend.position = "bottom") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#combined plots
ggsave(here("output", "carr_description.png"),
       plot = ((A + B + C) / D + plot_layout(nrow = 2, heights = c(1,1))),
       width = 13, height = 10, unit="in", dpi = 300)
