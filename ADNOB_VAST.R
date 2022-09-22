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

#Load data with Rstudio import from onedrive/documents/ADNOB_VAST.xlsx // ADNOB_ragon_NS.csv
#setwd("C:/Users/mjohnson/Documents/useful graphs/ADNP_OB")

setwd("~/ADNP/VAST_data")

AD <- read.csv(file = "ADNOB_Ragon.csv", stringsAsFactors = T)

AD <- rename(AD, NS = ADNOB)
levels(AD$Vaccine) <- c("Vi-PS", "Vi-TCV")

AD <- AD %>% 
  mutate(Time = fct_relevel(Time, "PV","28"))
                            

# Calculate Fold change --------------------------------------------------------

AD <- AD %>%
  group_by(Lab_ID) %>%
  mutate(FC =
           (NS[Time == "28"] - NS[Time == "PV"])
         / (NS[Time == "PV"]))

options(scipen = 999)

AD28 <- AD %>%
  filter(Time == "28")

AD28$FC <- 
  round(AD28$FC, digits = 3)

# Plot fold change and ADNOB ---------------------------------------------------

AD %>% group_by(Time) %>% ggplot(aes(x = Diagnosis, y = NS)) + geom_boxplot() + facet_grid(Vaccine ~ Time)

# By diagnosis
AD %>% filter(Time == "PV") %>% ggplot(aes(x = Diagnosis, y = FC)) + geom_boxplot(outlier.shape = NA) + ylim(-50, 50) + facet_wrap(~ Vaccine)

# Grouped
AD %>% filter(Time == "PV") %>% ggplot(aes(x = Diagnosis, y = FC, fill = Diagnosis)) + geom_boxplot(outlier.shape = NA) + ylim(-50, 50) +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) + labs(y = "ADNOB Fold Change\n") + theme_bw()

# Vaccine split 
AD %>% filter(Time == "PV") %>% ggplot(aes(x = Diagnosis, y = FC, fill = Diagnosis)) + geom_boxplot(outlier.shape = NA) + ylim(-50, 50) +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) + labs(y = "ADNOB Fold Change\n") + facet_wrap(~ Vaccine) + theme_bw()


# Day 0 v D28

AD %>%
  ggplot(aes(x = Time, y = NS, fill = Time)) +
  geom_jitter() +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) +
  labs(y = "ADNOB MFI") +
  facet_grid( ~ Vaccine) + theme_bw()

# Wilcoxon tests ---------------------------------------------------------------

my_comparisons = list(c("nTD", "TD"))

w_test <-
  compare_means(data = AD28,
                method = "wilcox.test",
                formula = FC ~ Diagnosis)
w_test # p = 0.0141 *

# Test Vi-PS and Vi-TCV stat differences separately 
# Make time point_vax factor

AD <-
  AD %>% mutate(Vax_Time = as.factor(str_c(Vaccine, Time, sep = "_")))

AD_rm <- na.omit(AD) # Removes samples with missing obvs for the w.test

AD_PS <- filter(AD_rm, Vaccine == "Vi-PS")
AD_TCV <- filter(AD_rm, Vaccine == "Vi-TCV")

PS.test <-
  wilcox.test(NS ~ Time,
              data = AD_PS,
              exact = FALSE,
              paired = TRUE)
PS.test # p-value = 0.8065

TCV.test <-
  wilcox.test(NS ~ Time,
              data = AD_TCV,
              exact = FALSE,
              paired = TRUE)

PS.test # p = 0.6 # p-value = 0.8065
TCV.test # p-value = 0.00147


# ADNOB 28 by vaccine

vax.test <- 
  wilcox.test(NS ~ Vaccine,
              data = AD28,
              exact = FALSE,
              paired = FALSE)

# FC = NS
vax.test = p-value = 0.04263

# Plot ADNOB ----------------------------------------------------------

# Tidy levels for plots
levels(AD_rm$Time)
levels(AD_rm$Time) <- c("Pre-Vac", "Post-Vac")
levels(AD$Time) <- c("Pre-Vac", "Post-Vac")
levels(AD28$Time) <- c("Pre-Vac", "Post-Vac")


# Outliers removed from plot with limits, one point both in PS and TCV
# Lets u see the distribution  better, but see what supervisors say

 ggboxplot(
  AD,
  x = "Time",
  y = "NS",
  fill = "Time",
  palette = c("#4045D5", "#40a9d5"),
  ylab = "ADNOB (MFI)\n",
  xlab = "Visit",
  outlier.shape = NA,
  add = "jitter",
  facet.by = "Vaccine") +
  grids(linetype = "dashed", color = "lightgrey") +
  theme(legend.position = "FALSE") +
  theme(axis.title = element_text(face = "bold"))

# I like this plot :)
# error.plot = "pointrange", double check what this iswidth = 0.25

# ggplot version with stat compare means
 
 # Visit ~ ADNOB ~ Vaccine -----------------
 

 my_comparisons <- list(c("Pre-Vac", "Post-Vac"))
 
 #4195 x 5953
 tiff(file = "ADNOB_vax.tiff",
      res = 720,
      width = 4300, height = 4000,
      compression = "lzw")
 AD %>% ggplot(aes(Time, NS)) +
   geom_boxplot(fill = c("#4045D5", "#40a9d5", "#4045D5", "#40a9d5"), outlier.shape = NA) +
   geom_jitter(width = 0.25, alpha = 0.9) +
   facet_wrap ("Vaccine") +
   labs(x = "Visit", y = "\n Oxidative Burst Score\n") +
   stat_compare_means(
     comparisons = my_comparisons,
     label = "p.signif",
     method = "wilcox.test",
     paired = TRUE
   )  + theme_pubr(border = TRUE) +
   theme(panel.grid.major = element_line(linetype = "dashed", color = "lightgrey"), axis.title = element_text(face ="bold"))
 dev.off()
 


 # geom_line(aes(group = Lab_ID))
 # Line plot some samples go down post-vac, technical - could be background in NS control (was a corner well rip)
 # So they probably don't do anything, or could be a biological thing
 

 
 
 
 
 b <- AD28 %>%
  ggboxplot(
    x = "Diagnosis",
    y = "FC",
    fill = "Diagnosis",
    palette = c("seagreen3", "sandybrown"),
    outlier.shape = NA,
    ylab = "\nADNOB (Fold Change)\n"
  ) +
  ylim(-20, 30) +
  grids(linetype = "dashed", color = "lightgrey") +
  ggtitle("") +
  theme(
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA),
    plot.title = element_textbox(
      hjust = 0.5,
      vjust = 1,
      margin = margin(t = 5, b = 5),
      size = 12
    )
  )  +
  geom_rect(
    xmin = 0,
    xmax = 100,
    ymin = 30,
    ymax = 35,
    fill = "grey92"
  ) +
  geom_hline(yintercept = 30, colour = "black") +
  theme(axis.title = element_text(face = "bold"))

(a | b) + plot_annotation(tag_levels = "a")



# Manual P-val adding -----------------------------------------------------


# Issue: when plotting sig bar is at 1000 fold because of outliers
# If I just remove these outliers I think the sig might change
# potential 
# Doesn't plot even manually if ylim is used
# Ignore outliers in box plot first


ps.test <-
  wilcox.test(NS ~ Time,
              data = AD_PS,
              exact = FALSE,
              paired = TRUE)$p.value

tcv.test <-
  wilcox.test(NS ~ Time,
              data = AD_TCV,
              exact = FALSE,
              paired = TRUE)$p.value



# (B) FC Diagnosis - combined

fc.test <-
  wilcox.test(FC ~ Diagnosis,
              data = AD28,
              alternative = "two.sided")$p.value

# (c) Fc diagnosis, split by vaccine

ps.test.fc <-
  wilcox.test(FC ~ Diagnosis,
              data = AD_PS,
              exact = FALSE,
              paired = FALSE)$p.value

tcv.test.fc <-
  wilcox.test(FC ~ Diagnosis,
              data = AD_TCV,
              exact = FALSE,
              paired = FALSE)$p.value


# Round all test values
names(.GlobalEnv)
Pattern1 <- grep("\\.test", names(.GlobalEnv), value = TRUE)
Pattern1_list <- do.call("list", mget(Pattern1))

# Write over and modify all p vals in list
for (i in seq_along(Pattern1_list)) {
  Pattern1_list[[i]] <- signif(Pattern1_list[[i]], digits = 3)
}

p_vals <- as.data.frame(Pattern1_list)
# To do: add if else to add asterix

df_p_val <- data.frame(
  group1 = "1", #name of x var
  group2 = "2", #name of y var
  label = result,# test name
  y.position = 6 # set as nul
)


AD28 %>%
  ggboxplot(
    x = "Diagnosis",
    y = "FC",
    fill = "Diagnosis",
    palette = c("sandybrown", "seagreen3"),
    outlier.shape = NA,
    ylab = "\nADNOB (Fold Change)\n",
    group = "Diagnosis",
    error.plot = "linerange"
  ) +
  ylim(-30, 50) +
  add_pvalue(
    p_val,
    xmin = "group1",
    xmax = "group2",
    label = "label",
    y.position = 45,
    tip.length = 2,
  ) 


