
# Set up -----------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(tidylog)
library(ggpubr)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(forcats)
library(ggtext)
library(ggprism)
library(rstatix)

setwd("~/ADNP/TyVAC_data")
load("tyvac.RData")

my_comparisons <- list(c("Pre-Vac", "Post-Vac"))

# Fig 1 ------------------------------------------------------------------------
phago_data %>%
  filter(Sample_ID != "PP", Sample_ID != "NS", Sample_ID != "NB") %>%
  ggplot(aes(Visit, mean_phago)) + xlab("Visit") + ylab ("Phagocytic Score\n") +
  geom_boxplot(fill = c("#6A4C7F", "#CE7DA5"),
               outlier.shape = NA)  +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.signif") +
  geom_jitter(size = 1.2,
              width = 0.2,
              height = 0.1) + theme_pubr(border = TRUE) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"),
    axis.title = element_text(face = "bold")
    )

# Fig 2 ------------------------------------------------------------------------
load(file = "Fig2_Vast_TyVAC.RData")

new_fig2 %>%
  ggplot(aes(Visit, mean_phago)) +
  geom_boxplot(fill = c("#6A4C7F", "#CE7DA5", "#6A4C7F", "#CE7DA5")) +
  geom_point() +
  my_theme +
  geom_line(aes(group = Sample_ID)) + facet_wrap ("Study") +
  labs(x = "Visit", y = "\n Phagocytic Score", title = "ADNP") + theme(panel.background = element_rect(color = "#6a6b65")) +
  ylim(0, 10500) + stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    test = "wilcoxon.test",
    paired = TRUE
  ) + theme_pubr(border = TRUE) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"),
    axis.title = element_text(face = "bold")
  )


# Fig 3 ------------------------------------------------------------------------
#load(file = "tyvac_paper.Rdata")

my_comparisons <- list( c("0-4", "5-9"), c("0-4", "10-15"), c("5-9", "10-15"))
#PreVAC
pre <-
  phago_data %>%
  filter(Visit == "Pre-Vac", !is.na(Age_Group)) %>%
  ggplot(aes(x = Age_Group, y = mean_phago)) +
  geom_boxplot(fill= c("#00AFBB", "#E7B800", "#FC4E07"), outlier.shape = NA) +
  geom_jitter(size = 1.2,
              width = 0.2,
              height = 0.1) +
  xlab("Age (years)") +
  ylab ("Phagocytic Score\n")  +
  ggtitle("Pre-vaccination") +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     vjust = 0.5, size = 5) +
  theme_pubr(border = TRUE) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"),
    axis.title = element_text(face = "bold"), legend.position = "none"
  )

# Post Vac
post <-
  phago_data %>%
  filter(Visit == "Post-Vac", !is.na(Age_Group)) %>%
  ggplot(aes(x = Age_Group, y = mean_phago)) +
  geom_boxplot(fill= c("#00AFBB", "#E7B800", "#FC4E07"), outlier.shape = NA) +
  geom_jitter(size = 1.2,
              width = 0.2,
              height = 0.1) +
  xlab("Age (years)") +
  ylab ("Phagocytic Score\n")  +
  ggtitle("Post-vaccination") +
  theme_pubr(border = TRUE) +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"),
    axis.title = element_text(face = "bold"), legend.position = "none"
  )

# Fold Change 

# Just plot sig comparison brackets
my_comparisons <- list( c("0-4", "5-9"), c("0-4", "10-15"))
Fc <- 
  phago_data %>%
  filter(Visit == "Post-Vac", !is.na(Age_Group)) %>%
  ggplot(aes(x = Age_Group, y = Foldchange)) +
  geom_boxplot(fill= c("#00AFBB", "#E7B800", "#FC4E07"), outlier.shape = NA) +
  geom_jitter(size = 1.2,
              width = 0.2,
              height = 0.1) +
  xlab("Age (years)") +
  ylab ("Fold change \n")  +
  ggtitle("ADNP fold-change") +
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     vjust = 0.6, size = 4) +
  theme_pubr(border = TRUE) + 
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"),
    axis.title = element_text(face = "bold"), legend.position = "none"
  )

#regroup plots
all_plots <-
  (pre | post / Fc) + plot_annotation(tag_levels = "a")

all_plots



ggsave(filename = "age_phago_box.png")


