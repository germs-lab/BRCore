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


BC_ranked <- spatial_core[[2]]
BC_ranked$core <- factor(ifelse(BC_ranked$IncreaseBC > 1 + 0.01*(increase_value), "Core", "Non Core Taxa"))

# Plot IncreaseBC with proportionBC
BC_ranked %>%
  ggplot(aes(x=proportionBC, y = (IncreaseBC -1)*100, fill = core)) +
  geom_point(pch = 21, alpha = 1, size = 2.5) +
  geom_hline(
    yintercept = 2,
    lty = 4,
    col = "darkred",
    cex = 0.5
  ) +
  annotate(
    geom = "text",
    x = max(BC_ranked$proportionBC - 0.1),
    y = max((BC_ranked$IncreaseBC -1)*100),
    label = paste("Last 2% increase\n(", length(as.numeric(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1 + 0.01*(increase_value))])), " OTUs)", sep = ""),
    color = "darkred",
    size = 4,
  ) +
  xlab("Contribution to Bray-Curtis Dissimilarity") +
  ylab("% Increase in Bray-Curtis Dissimilarity") +
  scale_fill_npg(
    name = "Core Membership",
    labels = c("Core Taxa", "Non Core Taxa")
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
model_fit <- sncm.fit(spp, taxon, pool = NULL)
model_statistics <- model_fit[[1]]
model_predictions <- model_fit[[2]]



model_predictions$fit_class <- "As predicted"
model_predictions[which(model_predictions$freq < model_predictions$pred.lwr), "fit_class"] <-
  "Below prediction"
model_predictions[which(model_predictions$freq > model_predictions$pred.upr), "fit_class"] <-
  "Above prediction"
model_predictions[which(is.na(model_predictions$freq)), "fit_class"] <- "NA"

model_predictions <- tibble::rownames_to_column(model_predictions, "OTU_ID")

occ_abun_model_pred <- as.data.frame(left_join(occ_abun, model_predictions, by = "OTU_ID"))

occ_abun_model_pred <- occ_abun_model_pred %>%
  mutate(fill_fit_class = paste0(fill, ":", fit_class))

occ_abun_model_pred <- occ_abun_model_pred %>%
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


occ_abun_model_pred %>%
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
    size = 0.8,
    aes(y = freq.pred, x = log10(p)),
    alpha = 0.55
  ) +
  geom_line(
    color = "black",
    lty = "twodash",
    size = 0.9,
    aes(y = pred.upr, x = log10(p)),
    alpha = 0.55
  ) +
  geom_line(
    color = "black",
    lty = "twodash",
    size = 0.9,
    aes(y = pred.lwr, x = log10(p)),
    alpha = 0.55
  ) +
  labs(x = "Log10(mean abundance)", y = "Occupancy") +
  annotate(
    "text",
    -Inf,
    Inf,
    label = paste("italic(R)^2 ==", round(model_statistics$Rsqr, 3)),
    parse = TRUE,
    size = 4.8,
    hjust = -0.2,
    vjust = 1.2
  ) +
  annotate(
    "text",
    -Inf,
    Inf,
    label = paste("italic(m) ==", round(model_statistics$m, 3)),
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

