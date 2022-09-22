library(dplyr)
library(ggplot2)
library(data.table)
library(tidylog)
library(ggpubr)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(patchwork)

#Load data with Rstudio import from onedrive/documents/ADNOB_VAST.xlsx // ADNOB_ragon_NS.csv


setwd("C:/Users/mjohnson/Documents/useful graphs/ADNP_OB")

AD <- read.csv(file = "ADNOB_Ragon.csv")

#AD$Vaccine <- as.factor(AD$Vaccine)
#AD$Diagnosis <- as.factor(AD$Diagnosis)
AD$Time <- as.factor(AD$Time)
#NS = ADNOB - NS control


AD <- AD %>% group_by(Lab_ID) %>%
  mutate(FC =
           (NS[Time == "28"] - NS[Time == "PV"])
         / (NS[Time == "PV"]))

#Plot fold change ####

AD %>% group_by(Time) %>% ggplot(aes(x = Diagnosis, y = ADNOB)) + geom_boxplot() + facet_grid(Vaccine ~ Time)

AD %>% filter(Time == "PV") %>% ggplot(aes(x = Diagnosis, y = FC)) + geom_boxplot(outlier.shape = NA) + ylim(-50, 50) + facet_wrap( ~
                                                                                                                                      Vaccine)


AD %>% filter(Time == "PV") %>% ggplot(aes(x = Diagnosis, y = FC, fill = Diagnosis)) + geom_boxplot(outlier.shape = NA) + ylim(-50, 50) +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) + labs(y = "ADNOB Fold Change\n")

#Wilcoxon test ####

my_comparisons = list(c("nTD", "TD"))


#+ stat_compare_means( comparisons = my_comparisons,method = "wilcox.test", label = "p.signif")


#regular stat_compare_means wouldn't adjust bracket position

AD28 <- AD %>% filter(Time == "28")
options(scipen = 999)

AD28$FC <- round(AD28$FC, digits = 3)

w_test <-
  compare_means(data = AD28,
                method = "wilcox.test",
                formula = FC ~ Diagnosis)

w_test <- w_test %>%
  mutate(y.position = c(45), x = "Diagnosis")

AD28 %>% ggplot(aes(x = Diagnosis, y = FC, fill = Diagnosis)) + geom_boxplot(outlier.shape = NA) + ylim(-50, 50) +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) + labs(y = "ADNOB Fold Change\n") + stat_pvalue_manual(w_test)

#Error in FUN(X[[i]], ...) : object 'Diagnosis' not found #Works if I remove fill for some reason
# Make plots without, Add in via powerpoint for now

levels(AD$Vaccine) <- c("Vi-PS", "Vi-TCV")
levels(AD28$Vaccine) <- c("Vi-PS", "Vi-TCV")
levels(AD_rm$Vaccine) <- c("Vi-PS", "Vi-TCV")


b <- AD28 %>% ggplot(aes(x = Diagnosis, y = FC, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) + ylim(-50, 50) +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) +
  labs(y = "\nADNOB Fold Change") + theme_bw() #+ theme(legend.position = "none") ,title = "ADNOB Fold Change
b

c <- AD %>% filter(Time == "28") %>%
  ggplot(aes(x = Diagnosis, y = FC, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) + ylim(-50, 50) +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) +
  labs(y = "ADNOB Fold Change") +
  facet_wrap( ~ Vaccine) + theme_bw()

c

w_test_vax <-
  compare_means(
    data = AD28,
    method = "wilcox.test",
    formula = FC ~ Diagnosis,
    group.by = group
  )
group <- c("Sample.ID", "Vaccine")

#Test diff in PV and Day 28 ####

#make time point vax factor
AD$Vaccine <- as.character(AD$Vaccine)
AD$Time <- as.character(AD$Time)
AD <-
  AD %>% mutate(Vax_Time = as.factor(str_c(Vaccine, Time, sep = "_")))

group2 <- c("Sample.ID", "Time")
AD$Time <- as.factor(AD$Time)
w_test_Time <- AD %>%
  compare_means(method = "wilcox.test",
                formula = NS ~ Time,
                group.by = group2)


#Attempt 2 ####


AD_rm <- na.omit(AD)


AD_rm %>%
  ggplot(aes(x = Time, y = NS, fill = Time)) +
  geom_jitter() +
  scale_fill_manual(values = c("seagreen3", "sandybrown")) +
  labs(y = "ADNOB MFI") +
  facet_grid( ~ Vaccine) + theme_bw()

AD_PS <- filter(AD_rm, Vaccine == "Vi-PS")
AD_TCV <- filter(AD_rm, Vaccine == "Vi-TT")
PS.test <-
  wilcox.test(NS ~ Time,
              data = AD_PS,
              exact = FALSE,
              paired = TRUE)
TCV.test <-
  wilcox.test(NS ~ Time,
              data = AD_TCV,
              exact = FALSE,
              paired = TRUE)

PS.test # p = 0.6
TCV.test # p = 0.000008


levels(AD_rm$Time) <- c("Day 0", "Day -28")
AD_rm$Time <-
  factor(AD_rm$Time,
         order = T,
         levels = c("Day -28", "Day 0"))
levels(AD_rm$Time) <- c("Pre-Vac", "Post-Vac")


a <- ggboxplot(
  AD_rm,
  x = "Time",
  y = "NS",
  color = "Time",
  palette = c("chocolate2", "#E7B800"),
  ylab = "ADNOB (MFI)\n",
  add = "jitter",
  add.params = list(size = 1.4, jitter = 0.2),
  facet.by = "Vaccine"
) + grids(linetype = "dashed", color = "lightgrey")  + ylim(-100, 8500) + theme(legend.position = "FALSE") + scale_fill_manual("chocolate2", "#E7B800")  +
  scale_alpha_manual(values = c(1, 0.7)) + labs(title = "ADNOB 28-days following Typhoid Vaccination")

a

pdf(file = "ADNOB_VAST.pdf")
a / (b | c) + plot_annotation(tag_levels = "a")
dev.off()

# Outliers removed from plot with limits, one point both in PS and TCV
# Lets u see the distribution  better, but see what supervisors say

a2 <- ggboxplot(
  AD_rm,
  x = "Time",
  y = "NS",
  fill = "Time",
  palette = c("#4045D5", "#4090D5"),
  ylab = "ADNOB (MFI)\n",
  outlier.shape = NA,
  facet.by = "Vaccine") +
  grids(linetype = "dashed", color = "lightgrey") +
  ylim(-1000, 3500) +
  theme(legend.position = "FALSE") +
  theme(axis.title = element_text(face = "bold"))

#+ labs(title = "ADNOB 28-days following Typhoid Vaccination")


a2

b <- AD28 %>%
  ggboxplot(
    x = "Diagnosis",
    y = "FC",
    fill = "Diagnosis",
    palette = c("seagreen3", "sandybrown"),
    outlier.shape = NA,
    ylab = "\nADNOB (Fold Change)\n"
  ) +
  ylim(-50, 50) +
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
    ymin = 50,
    ymax = 55,
    fill = "ivory2"
  ) +
  geom_hline(yintercept = 50, colour = "black") +
  theme(axis.title = element_text(face = "bold"))

(a2 | b) + plot_annotation(tag_levels = "a")



# Statistical significance - do nTD TD first with stars -----------

# TD v nTD ANDOB ####

AD28 %>%
  ggboxplot(
    x = "Diagnosis",
    y = "FC",
    fill = "Diagnosis",
    palette = c("seagreen3", "sandybrown"),
    outlier.shape = NA,
    ylab = "\nADNOB (Fold Change)\n",
    error.plot = "errorbar", group = "Diagnosis")  + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif", label.y.npc = "right", label.x.npc = "centre", na.rm = TRUE) +
  ylim(-50,100)


AD28 %>%
  ggplot(aes(x = Diagnosis, y = FC, group = Diagnosis)) + geom_boxplot(outlier.shape = NA)  + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif", label.y = 50, na.rm = TRUE) 
         
    #fill = "Diagnosis",
    palette = c("seagreen3", "sandybrown"),
    outlier.shape = NA,
    ylab = "\nADNOB (Fold Change)\n",
    error.plot = "errorbar", group = "Diagnosis")  + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif", label.y.npc = "right", label.x.npc = "centre", na.rm = TRUE) +
  ylim(-50,100)

#+
  ylim(-50,100) 

# I think it could be a standard form issue
# Convert FC out of SF options(scipen = 999)
# label.y = 5000, na.rm = TRUE) + ylim(-50,5000)


my_comparisons = list(c("nTD", "TD"))
# Issue: when plotting sig bar is at 1000 fold because of outliers
# If I just remove these outliers I think the sig might change
# potential 
# Doesn't plot even manually if ylim is used
# Ignore outliers in box plot first
?ggboxplot
?stat_compare_means

inherit.aes = TRUE


AD28 %>%
  ggboxplot(
    x = "Diagnosis",
    y = "FC",
    fill = "Diagnosis",
    palette = c("seagreen3", "sandybrown"),
    outlier.shape = NA,
    ylab = "\nADNOB (Fold Change)\n",
    error.plot = "errorbar", group = "Diagnosis")  + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test", label = "p.signif", label.y.npc = "right", label.x.npc = "centre", na.rm = TRUE, inherit.aes = FALSE) +
  ylim(-50,100)






#+
  geom_hline(yintercept = 50, colour = "black") +
  theme(axis.title = element_text(face = "bold")) +
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
    ymin = 50,
    ymax = 55,
    fill = "ivory2"
  ) +
  geom_hline(yintercept = 50, colour = "black") +
  theme(axis.title = element_text(face = "bold")) 


save.image(file = "ADNOB_sept22.RData")













# Test

b + 
  ggtitle("setosa")  + 
  coord_cartesian(clip = 'off', expand = 0.05) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)
        
        # ? element_textbox
        # ?  ? panel.border
        # ? geom_rect
        # ? theme

b

library(ggtext)

pdf(file = "ADNOB_VAST2.pdf")
a2/(b|c) + plot_annotation(tag_levels = "a")
dev.off()

a2/(b|c) + plot_annotation(tag_levels = "a")
ggsave(filename = "ADNOB_VAST.tiff")

(a2|b) + plot_annotation(tag_levels = "a")
ggsave(filename = "ADNOB_VAST2.tiff")
