# Setup
source("R/utils/000_setup.R")
load(
  "tests/data/test_phyloseq.rda"
)

# Test data set
## Making it small to it can run fast

test_small_phyloseq <- test_phyloseq %>%
  prune_taxa(taxa_sums(.) > 5000 & taxa_sums(.) < 50000, .) %>%
  prune_samples(sample_sums(.) >= 1, .)

# test_large_phyloseq <- prune_samples(sample_sums(test_phyloseq) >= 20000, test_phyloseq)

# save(test_phyloseq, file = "tests/data/test_small_phyloseq.rda")


# test multi_rarefy.R function 
otu_table_rare <-
    multi_rarefy(physeq = test_small_phyloseq,
                 depth_level = 5000,
                 num_iter = 99)

# Recreate the phyloseq object with the rarefied otu_table
test_phyloseq_rare <-
    phyloseq(
        otu_table(
            otu_table_rare %>%
                column_to_rownames("SampleID") %>%
                t() %>%
                as.matrix() %>%
                as.data.frame(),
            taxa_are_rows = TRUE
        ),
        test_small_phyloseq@sam_data,
        test_small_phyloseq@tax_table
    ) %>%
    prune_taxa(taxa_sums(x = .) > 0, x = .) %>%
    prune_samples(sample_sums(x = .) > 0, x = .)


# Extract the 'spatial' core microbiome across all sites. The 'Var' in the extract_core is 'site'.

spatial_core <- extract_core(
    test_small_phyloseq,
    Var = "site",
    method = "increase",
    increase_value = 2
) # Minimum seq depth was ~10,000 reads.


# Plot Bray-Curtis Dissimilarity Curve:
max <- 100 # Number of ranked-OTUs to plot
BC_ranked_max <- BC_ranked[1:max, ]


BC_ranked_max %>%
  ggplot(aes(x = rank[1:max], y = proportionBC)) +
  geom_point(
    pch = 21,
    col = "black",
    alpha = 0.85,
    size = 3.5
  ) +
  geom_vline(
    xintercept = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
    lty = 4,
    col = "darkred",
    cex = 0.5
  ) +
  scale_x_discrete(
    limits = BC_ranked$rank[1:max], # making the x axis more readable
    breaks = BC_ranked$rank[1:max][seq(1, length(BC_ranked$rank[1:max]),
      by =
        10
    )]
  ) +
  xlab("Ranked OTUs") +
  ylab("Contribution to Bray-Curtis Similarity") +
  annotate(
    geom = "text",
    x = last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1.02)])),
    y = 0.18,
    label = paste("Last 2% increase\n(", last(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >=
      1.02)])), " OTUs)", sep = ""),
    color = "darkred",
    size = 4,
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(.5, 1, .5, .5), "cm")
  )

# Plot Abundance Occupancy curve
occ_abun <- spatial_core[[4]]

occ_abun %>%
  ggplot(aes(
    x = log10(otu_rel),
    y = otu_occ,
    fill = fill
  )) +
  scale_fill_npg(
    name = "Core Membership",
    labels = c("Core Taxa", "Non Core Taxa")
  ) +
  geom_point(pch = 21, alpha = 1, size = 2.5) +
  labs(x = "log10(Mean Relative Abundance)", y = "Occupancy") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(.5, 1, .5, .5), "cm")
  )

# Fit Abundance-Occupancy Distribution to a Neutral Model
# Fit neutral model
taxon <- spatial_core[[7]]
spp <- t(spatial_core[[5]])
names(occ_abun)[names(occ_abun) == "otu"] <- "OTU_ID"
# source community pool
meta <- spatial_core[[6]]

# fitting model
obs.np <- sncm.fit(spp, taxon, stats = FALSE, pool = NULL)
sta.np <- sncm.fit(spp, taxon, stats = TRUE, pool = NULL)

obs.np$fit_class <- "As predicted"
obs.np[which(obs.np$freq < obs.np$pred.lwr), "fit_class"] <-
  "Below prediction"
obs.np[which(obs.np$freq > obs.np$pred.upr), "fit_class"] <-
  "Above prediction"
obs.np[which(is.na(obs.np$freq)), "fit_class"] <- "NA"

obs.np <- tibble::rownames_to_column(obs.np, "OTU_ID")
as.data.frame(left_join(occ_abun, obs.np, by = "OTU_ID")) -> fit_table
#
sta.np$above.pred <-
  sum(obs.np$freq > (obs.np$pred.upr), na.rm = TRUE) / sta.np$Richness
sta.np$below.pred <-
  sum(obs.np$freq < (obs.np$pred.lwr), na.rm = TRUE) / sta.np$Richness
fit_res <- as.data.frame(sta.np)
rownames(fit_res) <- "Value"
fit_res
obs2 <- as.data.frame(list_tab[[1]])

obs1 <- obs1 %>%
  mutate(fill_fit_class = paste0(fill, ":", fit_class))

obs1 <- obs1 %>%
  mutate(across(
    "fill_fit_class",
    str_replace,
    "no:Below prediction",
    "Non Core Taxa"
  )) %>%
  mutate(across(
    "fill_fit_class",
    str_replace,
    "no:Above prediction",
    "Non Core Taxa"
  )) %>%
  mutate(across(
    "fill_fit_class",
    str_replace,
    "no:As predicted",
    "Non Core Taxa"
  ))


obs1 %>%

  ggplot(aes(x = log10(otu_rel), y = otu_occ)) +
  scale_fill_npg(
    name = "Core Membership: Model Predictions",
    labels = c(
      "Core: Above Prediction",
      "Core: As Predicted",
      "Core: Below Prediction",
      "Non-Core Taxa"
    )
  ) +
  geom_point(
    aes(fill = fill_fit_class),
    pch = 21,
    alpha = 0.75,
    size = 2.2
  ) +
  geom_line(
    color = "red",
    data = obs1,
    size = 0.8,
    aes(y = obs1$freq.pred, x = log10(obs1$p)),
    alpha = 0.55
  ) +
  geom_line(
    color = "black",
    lty = "twodash",
    size = 0.9,
    data = obs1,
    aes(y = obs1$pred.upr, x = log10(obs1$p)),
    alpha = 0.55
  ) +
  geom_line(
    color = "black",
    lty = "twodash",
    size = 0.9,
    data = obs1,
    aes(y = obs1$pred.lwr, x = log10(obs1$p)),
    alpha = 0.55
  ) +
  labs(x = "Log10(mean abundance)", y = "Occupancy") +
  annotate(
    "text",
    -Inf,
    Inf,
    label = paste("italic(R)^2 ==", round(obs2$Rsqr, 3)),
    parse = TRUE,
    size = 4.8,
    hjust = -0.2,
    vjust = 1.2
  ) +
  annotate(
    "text",
    -Inf,
    Inf,
    label = paste("italic(m) ==", round(obs2$m, 3)),
    parse = TRUE,
    size = 4.8,
    hjust = -0.2,
    vjust = 3.2
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14),
    title = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.margin = unit(c(.5, 1, .5, .5), "cm")
  )

    ggplot(aes(x = log10(otu_rel), y = otu_occ)) +
    scale_fill_npg(
        name = "Core Membership: Model Predictions",
        labels = c(
            "Core: Above Prediction",
            "Core: As Predicted",
            "Core: Below Prediction",
            "Non-Core Taxa" )
    ) +
    geom_point(
        aes(fill = fill_fit_class),
        pch = 21,
        alpha = 0.75,
        size = 2.2
    ) +
    geom_line(
        color = "red",
        data = obs1,
        size = 0.8,
        aes(y = obs1$freq.pred, x = log10(obs1$p)),
        alpha = 0.55
    ) +
    geom_line(
        color = "black",
        lty = "twodash",
        size = 0.9,
        data = obs1,
        aes(y = obs1$pred.upr, x = log10(obs1$p)),
        alpha = 0.55
    ) +
    geom_line(
        color = "black",
        lty = "twodash",
        size = 0.9,
        data = obs1,
        aes(y = obs1$pred.lwr, x = log10(obs1$p)),
        alpha = 0.55
    ) +
    labs(x = "Log10(mean abundance)", y = "Occupancy") +
    annotate(
        "text",
        -Inf,
        Inf,
        label = paste("italic(R)^2 ==", round(obs2$Rsqr, 3)),
        parse = TRUE,
        size = 4.8,
        hjust = -0.2,
        vjust = 1.2
    ) +
    annotate(
        "text",
        -Inf,
        Inf,
        label = paste("italic(m) ==", round(obs2$m, 3)),
        parse = TRUE,
        size = 4.8,
        hjust = -0.2,
        vjust = 3.2
    ) +
    theme_bw() +
    theme(
        axis.title.x = element_text(size = 14),
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.margin = unit(c(.5, 1, .5, .5), "cm")
    )

