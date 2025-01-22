#written by Dan, Deus and Chikondi
#===============================================================================

#===============================================================================
#viral circulation plots
#===============================================================================

#load the viral data
viral_data <- 
  rio::import(here::here("data","viral_data.xlsx")) %>%
  dplyr::select(cdate, co_v, covid_19f, flu_a, flu_b, hmpv, paraflu1_2_3_4, rhino, rsv) %>%
  dplyr::mutate(yrmo = str_sub(cdate, 1,7),
                seas = if_else(yrmo %in% c('2021-02','2021-03','2021-04','2021-05','2021-06'), 'Spring 2021',
                               if_else(yrmo %in% c('2021-11','2021-12','2022-01','2022-02','2022-03','2022-04','2022-05','2022-06'), 'Winter/Spring 2021/22', NA_character_)),
                co_v = as.integer(co_v), 
                covid_19f = as.integer(covid_19f), 
                flu_a = as.integer(flu_a), 
                flu_b = as.integer(flu_b), 
                hmpv = as.integer(hmpv), 
                paraflu1_2_3_4 = as.integer(paraflu1_2_3_4), 
                rhino = as.integer(rhino), 
                rsv = as.integer(rsv))

cols <- c('Spring 2021'='#F8766D','Winter/Spring 2021/22'='#00BA38')

A <-
  bind_rows(
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(co_v) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Seasonal Coronaviruses'),
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(covid_19f) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'SARs-CoV-2'),
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(flu_a) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Influenza A'),
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(flu_b) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Influenza B'),
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(hmpv) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Human Metapneumoviruses'),
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(paraflu1_2_3_4) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Parainfluenza 1,2,3,4'),
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(rhino) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Rhinoviruses'),  
    viral_data %>% dplyr::group_by(yrmo,seas) %>% dplyr::tally(rsv) %>% dplyr::ungroup() %>% dplyr::mutate(p=n/sum(n, na.rm = TRUE), pathogen = 'Respiratory Syncytial Virus')) %>%
  
  ggplot() +
  geom_line(aes(x = yrmo, y=n, group=1), size = 1.2) +
  geom_rect(aes(xmin = ('2021-01'), xmax = ('2021-06'), ymin = 0, ymax = Inf, color = 'Spring 2021', alpha= 0.004), alpha = 0.004, fill='#F8766D') +
  geom_rect(aes(xmin = ('2021-11'), xmax = ('2022-06'), ymin = 0, ymax = Inf, color = 'Winter/Spring 2021/22'), alpha = 0.004, fill='#00BA38') +
  scale_colour_manual(name = "Study period:", values = cols) + 
  facet_wrap(.~pathogen, nrow = 4, scales = 'free_y') +
  theme_bw(base_size = 16,  base_family = 'American Typewriter') + 
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5)) +
  labs(title = "", x = "Date of reporting", y = "Number of cases") +
  theme(strip.text.x = element_text(size = 16), strip.text.y = element_text(size = 16), strip.background = element_rect(fill = "gray80")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), legend.position = 'bottom') +
  theme(axis.text=element_text(size=20), legend.text = element_text(size=16))

ggsave(here("output", "viral.png"),
       plot = A,
       width = 24, height = 12, unit = "in", dpi = 300)
