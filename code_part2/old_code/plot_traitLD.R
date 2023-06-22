library(tidyverse)
library(data.table)

traits_table <- as_tibble(fread("../generated_data/traits_table.txt"))
pops <- c("AFR","AMR","CSA","EAS","EUR")

#
LD_long <- traits_table %>%
  select(prive_code, GWAS_trait_type, PGS_trait_type, all_of(c(paste0("traitLD_unadj_",pops)))) %>%
  pivot_longer(
    cols = paste0("traitLD_unadj_",pops),
    names_to = "pop",
    names_prefix = "traitLD_unadj_",
    values_to="traitLD_unadj")
LD_long %>%
  group_by(pop) %>%
  summarize(traitLD_unadj_mean = mean(traitLD_unadj),
            traitLD_unadj_sd = sd(traitLD_unadj)) %>%
  arrange(-traitLD_unadj_mean)

LD_long2 <- LD_long %>%
  left_join(
    traits_table %>% select(prive_code, short_label, n_sig_SNPs, portability_index,
                            traitLD_unadj_min,traitLD_unadj_max,traitLD_unadj_mean),
    by="prive_code"
  ) %>% filter(GWAS_trait_type=="quantitative", PGS_trait_type=="quantitative")
LD_long2 %>% group_by(pop) %>% summarize(traitLD_unadj = mean(traitLD_unadj)) %>% arrange(traitLD_unadj)

ggplot(LD_long2,aes(x=reorder(short_label, traitLD_unadj_mean), y = traitLD_unadj)) +
  geom_segment(aes(xend=reorder(short_label, traitLD_unadj_mean), y=traitLD_unadj_min, yend=traitLD_unadj_max)) +
  geom_point(aes(color=pop)) +
  geom_point(aes(y=traitLD_unadj_mean),color="black") +
  geom_line(aes(y=portability_index*-600), color="black") +
  #geom_text(aes(label=n_sig_SNPs, y=1.25), size=2.5) +
  #coord_flip() +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Trait (quantitative only)") +
  ylab("Trait LD-Score (not MAF-adjusted)") +
  labs()


# gini plot
gini_long <- traits_table %>%
  select(prive_code, GWAS_trait_type, PGS_trait_type, n_sig_SNPs,
         all_of(paste0("gini_",c(pops)))) %>%
  pivot_longer(
    cols = paste0("gini_",c(pops)),
    names_to = "pop",
    names_prefix = "gini_",
    values_to="gini") %>%
  filter(GWAS_trait_type == "quantitative") %>%
  left_join(traits_table %>% select(prive_code, gini_panUKB), by="prive_code") %>%
  group_by(prive_code) %>%
  mutate(gini_min = min(gini),
         gini_max = max(gini))

ggplot(gini_long,aes(x=reorder(prive_code, gini_panUKB), y = gini)) +
  geom_segment(aes(xend=reorder(prive_code, gini_panUKB), y=gini_min, yend=gini_max)) +
  geom_point(aes(color=pop)) +
  geom_point(aes(y=gini_panUKB),color="black") +
  #geom_text(aes(label=n_sig_SNPs_allpops, y=0), size=2.5, ) +
  coord_flip() +
  theme_light() +
  xlab("Trait (quantitative only)") +
  ylab("Gini") +
  ylim(0,1)
