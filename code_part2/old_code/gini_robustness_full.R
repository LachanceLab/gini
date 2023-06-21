# gini_robustness_full.R
# Calculates gini using different parameters to check robustness
# more thorough comparison that is not used in final paper (not supplemental)

### Libraries and directories ####
library(tidyverse)
library(data.table)
source("../code_part1/helper_functions/helper_functions.R")
source("../code_part1/helper_functions/winners_curse_functions.R")

# sets working directory
setwd("./")

# directory where panUKB GWAS summary statistics were saved to
dir_sf <- "../generated_data/panUKB_sf/"
# sets location to traits_table.txt
loc_table <- "../generated_data/traits_table.txt"

# reads traits_table and extract quantitative trait codes
traits_table <- as_tibble(fread(loc_table))
qcodes <- (traits_table %>% filter(GWAS_trait_type=="quantitative",
                                   PGS_trait_type=="quantitative"))$prive_code

# parameter options
beta_sources <- c("meta/EUR","pop-specific")
#AF_sources <- c("AF_1kG_","af_")
AF_source <- "af_"
thresholds <- c(50, 100, 250, 500, 750, 1000, -1)
top_SNPs <- c("meta/EUR","pop-specific")
pval_cutoff <- 1E-5

# makes table for collecting data
robustness_tbl <- tibble(
  prive_code = as.character(),
  pop = as.character(),
  beta_source = as.character(),
  beta_pop = as.character(),
  AF_source = as.character(),
  threshold = as.numeric(),
  top_SNP = as.character(),
  top_SNP_pop = as.character(),
  pval_cutoff = as.numeric(),
  gini = as.numeric(),
  n_sig_SNPs = as.numeric(),
  n_sig_SNPs_pop = as.numeric(),
  sum_gvc_top = as.numeric(),
  sum_gvc_all = as.numeric()
)
ii <- 0
# mega-loop
for (code in qcodes) {
  slice <- traits_table %>% filter(prive_code==code)
  filename <- paste0(code,"_sf_indep.txt")
  loc_summary_file <- paste0(dir_sf,filename)
  sf <- as_tibble(fread(loc_summary_file))
  n_sig_SNPs <- nrow(sf)
  
  # gets list of populations GWASs were done on
  beta_pops <- substring(colnames(sf %>% select(starts_with("beta_"), -contains("meta2use"))),6)
  af_pops <- substring(colnames(sf %>% select(starts_with("af_cases_"), -contains("meta2use"))),10)
  trait_type <- "binary"
  if (length(af_pops) == 0) {
    trait_type <- "quantitative"
    af_pops <- substring(colnames(sf %>% select(starts_with("af_", ignore.case = FALSE), -contains("meta2use"))),4)
  }
  af_pops <- beta_pops
  # defines columns to use in gvc calculations later
  if ("meta_hq" %in% beta_pops) { meta2use <- "meta_hq"
  } else if ("meta" %in% beta_pops) { meta2use <- "meta"
  } else { meta2use <- "EUR" }
  # shuffles order of af_pops such that meta & meta_hq come last
  af_pops2 <- c(af_pops[!grepl("meta",af_pops)],af_pops[grepl("meta",af_pops)])
  # extracts sample size
  for (pop in af_pops2) {
    col_AF <- paste0("af_",pop)
    af_cases_pop <- paste0("af_cases_",pop)
    af_controls_pop <- paste0("af_controls_",pop)
    
    # for non-meta pops, fills in NA values with 1kG AFs
    if ((pop != "meta_hq") & (pop != "meta")) {
      col_1kG_AF <- paste0("AF_1kG_",pop)
      if (col_AF %in% colnames(sf)) {
        sf[is.na(sf[[col_AF]]),col_AF] <- sf[is.na(sf[[col_AF]]),col_1kG_AF]
      }
      if (af_controls_pop %in% colnames(sf)) {
        sf[is.na(sf[[af_controls_pop]]),af_controls_pop] <- sf[is.na(sf[[af_controls_pop]]),col_1kG_AF]
        sf[is.na(sf[[af_cases_pop]]),af_cases_pop] <- sf[is.na(sf[[af_cases_pop]]),col_1kG_AF]
      }
    }
    
    # uses meta cases/controls ratio for meta_hq
    n_pop <- switch(pop,
                    "meta" = "full_cohort_both_sexes",
                    "meta_hq" = "full_cohort_both_sexes", #"hq_cohort_both_sexes",
                    pop)
    
    # extracts number of cases
    n_cases <- slice[1,paste0("n_cases_",n_pop)][[1]]
    # extracts number of controls
    if (trait_type=="binary") {
      # for meta/meta_hq, sums all population controls
      if ((pop == "meta_hq") | (pop == "meta")) {
        n_controls <- rowSums(slice %>% select(starts_with("n_controls_")), na.rm=TRUE)
      } else { n_controls <- slice[1,paste0("n_controls_",n_pop)][[1]] }
    } else { n_controls <- 0 }
    # gets total number of cases/controls
    n_total <- n_cases + n_controls
    
    # adds pop-specific case+control AFs based on weighted average AFs
    if (!(col_AF %in% colnames(sf))) {
      sf[,col_AF] <- (n_cases * sf[,af_cases_pop] + n_controls * sf[,af_controls_pop]) / n_total
    }
    
  }
  # filters to SNPs with AF-data for every pop
  col_afs <- paste0("af_",af_pops)
  sf <- sf %>% filter(if_all(all_of(col_afs), ~!is.na(.)))
  n_sig_SNPs_allpops <- nrow(sf)
  
  # Winners Curse Correction for every beta
  for (pop in beta_pops) {
    col_AF <- paste0("af_",pop)
    col_beta <- paste0("beta_",pop)
    col_se <- paste0("se_",pop)
    if (pop=="meta_hq") {n_pop <- slice[1,"n_cases_hq_cohort_both_sexes"][[1]]
    } else if (pop=="meta") {n_pop <- slice[1,"n_cases_full_cohort_both_sexes"][[1]]
    } else {n_pop <- sum(slice[1,paste0("n_",c("cases","controls"),"_",pop)], na.rm=TRUE)
    }
    sf_WC <- sf %>% mutate(discovery.n = n_pop) %>%
      select(discovery.beta = !!enquo(col_beta),
             discovery.se = !!enquo(col_se),
             discovery.n, discovery.freq = !!enquo(col_AF))
    NA_mask <- is.na(sf_WC$discovery.beta)
    sf_WC <- as_tibble(correct_winners_curse(as.data.frame(sf_WC[!NA_mask,]), pval_cutoff))
    sf[!NA_mask,col_beta] <- sf_WC$debiased.beta.mle
  }
  
  # gets top SNPs according to meta-analysis
  col_beta <- paste0("beta_",meta2use)
  col_AF <- paste0("af_",meta2use)
  sf_meta <- sf
  sf_meta[,"gvc"] <- get_gvc(sf, col_beta, col_AF)
  sf_meta <- sf_meta %>% arrange(-gvc) %>% filter(!is.na(gvc))
  meta_SNPs <- sf_meta$SNP
  
  # calculates gvc, top SNPs, and gini
  for (pop in af_pops) {
    col_AF <- paste0("af_",pop)
    for (beta_source in beta_sources) {
      if (beta_source=="meta/EUR") {col_beta <- paste0("beta_",meta2use)
      } else {col_beta <- paste0("beta_",pop)}
      
      sf_pop <- sf
      sf_pop[,"gvc"] <- get_gvc(sf, col_beta, col_AF)
      n_sig_SNPs_pop <- nrow(sf_pop[!is.na(sf_pop$gvc),])
      sum_gvc_all <- sum(sf_pop$gvc, na.rm=TRUE)
      
      for (threshold in thresholds) {
        for (top_SNP in top_SNPs) {
          if (top_SNP=="meta/EUR") {
            top_SNP_pop <- meta2use
            if (threshold==-1) {threshold <- nrow(sf_meta)}
            # gets top n='threshold' SNPs from the meta
            sf_pop_top <- sf_pop %>% filter(SNP %in% meta_SNPs)
            sf_pop_top <- sf_pop_top[match(meta_SNPs, sf_pop_top$SNP),]
            sf_pop_top <- sf_pop_top %>% filter(!is.na(gvc)) %>%
              filter(row_number() <= threshold)
          } else {
            top_SNP_pop <- pop
            sf_pop_top <- sf_pop
            if (threshold==-1) { threshold <- nrow(sf_pop_top) }
          }
          # pads with zeros (including for NA values)
          gvc_list <- pad_zeros(sf_pop_top$gvc, threshold)
          sum_gvc_top <- sum(gvc_list)
          # calculates gini
          gini <- get_gini(gvc_list)
          
          if (ii %% 50 == 0) {
            print("ii code pop beta_source threshold top_SNPs gini")
          }
          ii <- ii + 1
          print(paste(ii,code,pop,beta_source,threshold,top_SNP,round(gini,4),
                      sep=" "))
          
          robustness_tbl <- robustness_tbl %>% add_row(
            prive_code = code,
            pop = pop,
            beta_source = beta_source,
            beta_pop = substr(col_beta,6, nchar(col_beta)),
            AF_source = AF_source,
            threshold = threshold,
            top_SNP = top_SNP,
            top_SNP_pop = top_SNP_pop,
            pval_cutoff = pval_cutoff,
            gini = gini,
            n_sig_SNPs = n_sig_SNPs,
            n_sig_SNPs_pop = n_sig_SNPs_pop,
            sum_gvc_top = sum_gvc_top,
            sum_gvc_all = sum_gvc_all
          )
        }
      }
    }
  }
}

loc_out <- "../generated_data/gini_robustness_data.txt"
fwrite(robustness_tbl, loc_out, sep="\t")
robustness_tbl <- as_tibble(fread(loc_out))

robustness_tbl <- robustness_tbl %>%
  mutate(prop_gvc_top = sum_gvc_top / sum_gvc_all) %>%
  drop_na()

# Number of traits per number of populations
robustness_tbl %>% group_by(prive_code) %>% summarize(n_pops=n()/20) %>%
  group_by(n_pops) %>% summarize(n_traits=n())

robustness_tbl %>%
  filter(beta_source=="meta/EUR",top_SNP=="meta/EUR", pop==beta_pop, threshold %in% thresholds) %>%
  group_by(threshold) %>%
  mutate(gini_rank = order(order(gini, decreasing=FALSE))) %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
ggplot(aes(x=as.factor(threshold),y=gini_rank)) +
  geom_line(aes(group=prive_code, color=group_consolidated))

# pop + top_SNP
robustness_tbl %>%
  filter(beta_source=="meta/EUR",threshold==100, (pop!="meta" & pop!="meta_hq")) %>%
  group_by(top_SNP,pop) %>%
  mutate(gini_rank = order(order(gini, decreasing=FALSE))) %>%
ggplot(aes(x=as.factor(top_SNP),y=gini)) + # y=gini_rank
  geom_line(aes(group=prive_code, color=pop)) +
  facet_wrap(~pop)

# pop + beta_source
robustness_tbl %>%
  filter(top_SNP=="meta/EUR",threshold==100, (pop!="meta" & pop!="meta_hq")) %>%
  group_by(beta_source,pop) %>%
  mutate(gini_rank = order(order(gini, decreasing=FALSE))) %>%
ggplot(aes(x=as.factor(beta_source),y=gini)) +
  geom_line(aes(group=prive_code, color=pop)) +
  facet_wrap(~pop)

###
library(GGally)
pops <- c("AFR","AMR","CSA","EAS","EUR")
pops2use <- c("EUR","AMR","CSA","AFR","EAS","meta")

# pop-comparison, meta beta, meta top SNPs
robustness_tbl %>%
  filter(top_SNP=="meta/EUR",threshold==100, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
ggpairs(mapping=aes(color=group_consolidated),
        columns = paste0("gini_",pops2use))

# pop-comparison, meta beta, pop-specific top SNPs, threshold=100
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==100, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use))

# pop-comparison, pop-specific beta, pop-specific top SNPs
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==100, pop %in% pops2use,
         beta_source=="pop-specific") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use))
# pop-comparison, pop-specific beta, meta top SNPs
robustness_tbl %>%
  filter(top_SNP=="meta/EUR",threshold==100, pop %in% pops2use,
         beta_source=="pop-specific") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use))


##
# pop-comparison, meta beta, pop-specific top SNPs, threshold=50
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==50, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = 50")
# pop-comparison, meta beta, pop-specific top SNPs, threshold=100
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==100, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = 100")
# pop-comparison, meta beta, pop-specific top SNPs, threshold=250
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==250, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = 250")
# pop-comparison, meta beta, pop-specific top SNPs, threshold=500
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==500, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = 500")
# pop-comparison, meta beta, pop-specific top SNPs, threshold=750
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==750, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = 750")
# pop-comparison, meta beta, pop-specific top SNPs, threshold=1000
robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold==1000, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = 1000")
# pop-comparison, meta beta, pop-specific top SNPs, threshold=ALL
robustness_tbl %>%
  filter(top_SNP=="pop-specific",!threshold %in% thresholds, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  select(prive_code, pop, gini) %>%
  pivot_wider(names_from=pop, values_from=gini,names_prefix="gini_") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
  ggpairs(mapping=aes(color=group_consolidated),
          columns = paste0("gini_",pops2use)) +
  labs(title="threshold = ALL")

robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold %in% thresholds, pop %in% pops2use,
         beta_source=="meta/EUR") %>%
  group_by(threshold) %>%
  summarize(gini_mean = mean(gini),gini_sd = sd(gini),
            gini_min = min(gini), gini_max = max(gini),
            prop_gvc_top_mean = mean(prop_gvc_top))

robustness_tbl %>%
  filter(top_SNP=="pop-specific",threshold == 500, pop %in% pops2use,
         beta_source=="meta/EUR") %>% #filter(pop == "meta") %>%
  left_join(traits_table[,c("prive_code","group_consolidated")],by="prive_code") %>%
ggplot(aes(x=prop_gvc_top, y=gini)) +
  geom_point(aes(color=group_consolidated), alpha=0.5) +
  geom_smooth(method='lm') +
  facet_wrap(~pop) +
  xlim(0,1) + ylim(0,1) +
  labs(title="gini vs prop_gvc_top for different populations")
