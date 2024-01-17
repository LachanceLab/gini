library(tidyverse)
library(data.table)

loc_simphenos <- "../generated_data/simulations/simphenos_tbl.txt"

simphenos_tbl <- as_tibble(fread(loc_simphenos))

ggplot(simphenos_tbl, aes(x = gini_true, y=gini_obs_p8_r01)) +
  geom_abline(slope=1, color='black') +
  geom_point(aes(shape = as.factor(h2), color=as.factor(Mc))) +
  scale_x_continuous(limits = c(0,1), expand=c(0,0.05)) +
  scale_y_continuous(limits = c(0,1), expand=c(0,0.05)) +
  theme_light() +
  labs(x = 'True Gini',
       y = 'Observed Gini (clumping r2 = 0.2, p < 1E-5)',
       color = 'Polygenicity',
       shape = "Heritability")
ggsave('../generated_figures/gini.true_obs.p5_r2.png',
       width = 180, height=150, units='mm')


simphenos_tbl.plot <- simphenos_tbl %>%
  select(-contains('gini_obs_nopad')) %>%
  pivot_longer(cols = starts_with('gini_obs_'),
               names_prefix = 'gini_obs_',
               names_to = 'setting',
               values_to = 'gini_obs') %>%
  mutate(setting.text = ifelse(setting=='p5_r2','r2 < 0.2, p < 1x10^-5',
                               'r2 < 0.01, p < 5x10^-8'))

ggplot(simphenos_tbl.plot, aes(x = gini_true, y=gini_obs)) +
  geom_abline(slope=1, color='black') +
  geom_point(aes(shape = as.factor(h2), color=as.factor(Mc))) +
  facet_wrap(~setting.text) +
  scale_x_continuous(limits = c(0,1), expand=c(0,0.05)) +
  scale_y_continuous(limits = c(0,1), expand=c(0,0.05)) +
  theme_light() +
  theme(strip.background = element_rect(fill='gray30'),
        strip.text = element_text(color='white')) +
  labs(x = 'True Gini',
       y = 'Observed Gini post-GWAS/clumping',
       color = 'Polygenicity',
       shape = "Heritability")
ggsave('../generated_figures/gini.true_obs.png',
       width = 240, height=120, units='mm')
