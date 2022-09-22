# Study Comparison ADNOB ---------------------------------------------------------

phago_data <-
  read.csv(file = "samplecheck2nepal.csv",
           fileEncoding = "UTF-8-BOM",
           stringsAsFactors = T)

phago_data <- phago_data[, -1]

phago_data <-
  phago_data %>% mutate(Visit = fct_relevel(Visit, "Pre-Vac", "Post-Vac"))

my_comparisons <- list(c("Pre-Vac", "Post-Vac"))


p <- phago_data %>% filter(Rep == "1",
                           Sample_ID != "PP",
                           Sample_ID != "NS",
                           Sample_ID != "mean_OB") %>%
  ggplot(aes(Visit, mean_OB)) +
  geom_boxplot(fill = c("#4045D5", "#40a9d5", "#4045D5", "#40a9d5")) +
  geom_point() +
  geom_line(aes(group = Sample_ID)) + facet_wrap ("Study") +
  labs(x = "Visit", y = "\nOxidative Burst Score\n", title = "ADNOB\n") +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "wilcox.test",
    paired = TRUE
  )  + theme_pubr(border = TRUE) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"),
    axis.title = element_text(face = "bold")
  )

tiff(
  file = "ADNOB_study_comp.tiff",
  res = 720,
  width = 4900,
  height = 5153,
  compression = "lzw"
)
p
dev.off()

