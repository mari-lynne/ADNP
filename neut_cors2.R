# Packages and themes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(dplyr)
library(ggplot2)
library(data.table)
library(tidylog)
library(ggpubr)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(janitor)


# Theme 
white_cor <- theme(axis.text=element_text(size=12),
                   axis.title = element_text(size = 12.5),  #face="bold"
                   panel.background = element_rect(fill = "white", colour = "grey",
                                                   size = 0.5, linetype = "solid"),
                   panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey88"),
                   panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey88"))


# Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("~/ADNP/VAST_data")

neut <- read.csv("neut_titres.csv")
funs <- read.csv("Ragon_cors.csv")
funs <- funs[,-c(24:32)]
vac <- funs %>% dplyr::filter(Time.Point == "28")
join <- left_join(neut, vac, by = "Lab.ID")
join <- clean_names(join)


# Cor sections #
# ADNOB with ragon data - fc functions
# Neut counts with ADNOB
# NM ratio with ADNOB
# Neuts with time/bacteria quant

# ADNP cors (loop) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y <- join$adnp

# Get string of just names with fc functions 

for(i in fcnames){
  plt <- ggplot(join, aes_string(x=i, y = Y, color = "vaccine")) +
    geom_point(size=3.6, alpha=0.8) +
    stat_cor(method = "spearman")+
    scale_color_manual(values =  c("black"))+
    geom_hline(yintercept=0, size=0.06, color="black") + 
    geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
    yscale("log10", .format = TRUE)+
    xscale("log10", .format = TRUE)+
    ylab ("ADNP")+
    xlab(gsub("(\\_|\\-)", " ", i))+
    ggtitle("28 days post-vaccination")+
    white_cor
  print(plt)
  Sys.sleep(0.5)
}


# ADNP 
a <- ggplot(join, aes(x=fcg_r3a_vi, y=adnp)) +
  geom_point(color="#d56c40", size=2.5, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Phagocytic Score\n",
       x=expression(Fc*gamma~R3*alpha~binding~(MFI)))+
  xscale("log10", .format = TRUE)+
  white_cor

a


b <- ggplot(join, aes(x=fcg_r3b_vi, y=adnp)) +
  geom_point(color="#d58640", size=2.5, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Phagocytic Score Score\n",
       x=expression(Fc*gamma~R3*beta~binding~(MFI)))+
  xscale("log10", .format = TRUE)+
  white_cor


b

# ADNOB - FcR correlations (loop) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

colnames <- names(select_if(join, is.numeric))
join$vaccine <- as.factor(join$vaccine_x)


Y <- join$adnob
for(i in colnames){
  plt <- ggplot(join, aes_string(x=i, y = Y, color = "vaccine")) +
    geom_point(size=3.6, alpha=0.8) +
    stat_cor(method = "spearman")+
    scale_color_manual(values =  c("#00AFBB", "#E7B800"))+
    geom_hline(yintercept=0, size=0.06, color="black") + 
    geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
    yscale("log10", .format = TRUE)+
    xscale("log10", .format = TRUE)+
    ylab ("ADNOB")+
    xlab(gsub("(\\_|\\-)", " ", i))+
    ggtitle("28 days post-vaccination")+
    white_cor
  print(plt)
  Sys.sleep(0.5)
}

c <- ggplot(join, aes(x=fcg_r3a_vi, y=adnob)) +
  geom_point(color="#4045d5", size=2.5, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Oxidative Burst Score\n",
       x=expression(Fc*gamma~R3*alpha~binding~(MFI)))+
  xscale("log10", .format = TRUE)+
  white_cor

c

d <- ggplot(join, aes(x=fcg_r3b_vi, y=adnob)) +
  geom_point(color="#408fd5", size=2.5, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Oxidative Burst Score\n",
       x=expression(Fc*gamma~R3*beta~binding~(MFI)))+
  xscale("log10", .format = TRUE)+
  white_cor

d

(a|b)/(c|d) + plot_annotation(tag_levels = "a", title = "28 days post-vaccination")
ggsave(filename = "ADNP_OB_fcy_cors.tiff", device = "tiff", dpi = 720)


#ADNOB cor w Time to positive BC (both)
#ADNOB cor w FcGR3b (both)
#ADNOB cor w FcgR3a (Vi-PS, TCV ns)
#ADNOB neg cor with eosinophils -0.8 0.0018, Vi-PS
#ADNOB cor w n_m ration, both, but only sig in Vi-TCV
#ADNOB cor w Vi Neut count only trend seen in Vi-TCV

#Just assocs in TCV supports idea that neutrophil mediated vaccine immunity
#induced in Vi-TCV vaccination, but less so by Vi-PS
#Diffs in vaccine compisition/presentation etc. 

# FcyR binding cors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fc3a <- join[((!is.na(join$fcg_r3a_vi))&(!is.na(join$adnob))),]
fc3b <- join[((!is.na(join$fcg_r3b_vi))&(!is.na(join$adnob))),]

a <- ggplot(join, aes(x=fcg_r3a_vi, y=adnob)) +
  geom_point(color="black", size=3.2, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Oxidative Burst Score\n",
       x=expression(Fc*gamma~R3*alpha~binding~(MFI)))+
  xscale("log10", .format = TRUE)+
  white_cor

a

b <- ggplot(join, aes(x=fcg_r3b_vi, y=adnob)) +
  geom_point(color="black", size=3.2, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Oxidative Burst Score\n",
       x=expression(Fc*gamma~R3*beta~binding~(MFI)))+
  xscale("log10", .format = TRUE)+
  white_cor

b


(a|b) + plot_annotation(tag_levels = "a", title = "28 days post-vaccination")
ggsave(filename = "ADNOB_fcy_cors.tiff", device = "tiff", dpi = 720)


cor %>%
  ggplot(aes(x=adnob, y = time_to_positive_bc_hours)) +
  geom_point(size=3.2, alpha=0.8, color="black") +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  xscale("log10", .format = TRUE)+
  xlab ("Oxidative Burst Score")+
  ylab("Time for +blood culture (Hhrs)\n")+
  ggtitle("28 days post-vaccination")+
  white_cor


cor <- join[((!is.na(join$time_to_positive_bc_hours))&(!is.na(join$adnob))),
            c("lab_id","adnob","time_to_positive_bc_hours")]
cor.test(
  cor$adnob, cor$time_to_positive_bc_hours,
  method = "spearman",
)
cor.test(
  fc3a$adnob, fc3a$fcg_r3a_vi,
  method = "spearman", exact = FALSE,
)  # Still very significant p-value = 0.0001826 # 0.4479726 

cor.test(
  fc3a$adnob, fc3a$fcg_r3b_vi,
  method = "spearman", exact = FALSE,
) # 4.308e-05, 0.4845717 



a2 <- ggplot(join, aes(x=fcg_r3a_vi, y=adnob)) +
  geom_point(color="black", size=3.2, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Oxidative Burst Score\n",
       x=expression(Fc*gamma~R3*alpha~binding~(MFI)))+
  white_cor



b2 <- ggplot(fc3b, aes(x=fcg_r3b_vi, y=adnob)) +
  geom_point(color="black", size=3.2, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Oxidative Burst Score\n",
       x=expression(Fc*gamma~R3*beta~binding~(MFI)))+
  white_cor


(a2|b2) + plot_annotation(tag_levels = "a", title = "28 days post-vaccination")
ggsave(filename = "ADNOB_fcy_cor_nolog.tiff", device = "tiff", dpi = 720)


# ADNOB ~ CFU and time to diagnosis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ADNOB r- 0.68, p=0.05 (both, small n)

# Filter na obvs
cor2 <- join[((!is.na(join$bacterial_quantification_cfu_ml))&(!is.na(join$adnob))),]

a <- cor2 %>% 
  ggplot(aes(x=adnob, y = time_to_positive_bc_hours))+
  geom_point(color="#4045d5", size=3, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.2, color="black", fill="black")+
  xlab ("Oxidative Burst Score")+
  ylab("Time for +blood culture (Hhrs)\n")+
  white_cor


b <- cor2 %>% 
  ggplot(aes(x=adnob, y = bacterial_quantification_cfu_ml))+
  geom_point(color="#4045d5", size=3, alpha=0.7)+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.2, color="black", fill="black")+
  xlab ("Oxidative Burst Score")+
  ylab("Bacterial CFU/ml\n")+
  white_cor

# Furthermore, in those who were diagnosed with typhoid fever, higher oxidative burst scores following vaccination were associated with longer time to a positive blood culture, and lower bacterial load 

(a|b) + plot_annotation(tag_levels = "a")
ggsave(filename = "ADNOB_symptoms.tiff", device = "tiff", dpi = 720)



#Neutrophil counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot(
  join,aes(y=time_to_positive_bc_hours, x = neutrophils_x10_9_l))+
  geom_point(size=3.6, alpha=0.8)+
  scale_color_manual(values =  c("#00AFBB", "#E7B800"))+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black")+ 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(x=expression(Neutrophils~x10^9~"/"~L),
       y="Time for +blood culture (Hhrs)\n")+
  white_cor
#only correlated in the TCV group 0.65, p0.02

ggplot(join, aes(y=bacterial_quantification_cfu_ml, x = neutrophils_x10_9_l,
                 color = vaccine))+
  geom_point(size=3.6, alpha=0.8)+
  scale_color_manual(values =  c("#00AFBB", "#E7B800"))+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black")+ 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(x=expression(Neutrophils~x10^9~"/"~L),
       y="bacterial CFU/ml")+
  white_cor
#Only correlated in the TCV group -0.59, p=0.044


#N/M ratio 
ggplot(
  join,aes(y=time_to_positive_bc_hours, x = n_m_ratio))+
  geom_point(size=3.6, alpha=0.8)+
  scale_color_manual(values =  c("#00AFBB", "#E7B800"))+
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black")+ 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  labs(y="Time for +blood culture (Hhrs)\n",
       x="N/M ratio")+
  white_cor



#bacterial_quantification_cfu_ml
#nothing interesting/sig with ratios


# Cytokines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cyt <- read.csv(file = "38plex.csv")
cyt <- clean_names(cyt)
cyt <- filter(cyt, vaccine != "Control")

#Pre-challenge cytokine cors with ADNOB
#day 0 cors with ADNOB
#Post-challenge cytokine cors with protection

pre <- left_join(cyt, join, by = "lab_id")
pre <- pre[,!str_ends(colnames(pre), ".y")] #clean cols
colnames(pre) <- str_replace_all(colnames(pre), ".x", "")
pre <- pre[,-44] #remove dup vaccine column

#CD40L FC ####
CD40 <- select(pre, lab_id, time_point, vaccine, diagnosis, s_cd40l, adnob)

#Fold change #
CD40 <- 
  CD40 %>%
  group_by(lab_id) %>%
  mutate(
    FC_27 = ((s_cd40l[time_point == "-27"] - s_cd40l[time_point == "PV"])/
               s_cd40l[time_point =="PV"]),
    FC_25 = ((s_cd40l[time_point == "-25"] - s_cd40l[time_point == "PV"])/
               s_cd40l[time_point =="PV"]))

# Correlation graphs
CD40 %>%
  filter(time_point == "-25") %>%
  ggplot(aes(x=FC_27, y = adnob)) +
  geom_point(color="#9b42f5", size=3.6, alpha=0.7) +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  ylab("ADNOB D28\n")+
  xlab("sCD40L FC D3")+
  ggtitle("Post-vaccination")+
  white_cor



# Day 0 cytokines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chal <- read.csv(file = "post-challenge.csv")
chal <- clean_names(chal)
#chal <- filter(chal, vaccine != "Control")

#Join chal with ADNOB join table first then merge
chal2 <- inner_join(chal, join, by = "lab_id")
  #drop_na(adnob) %>% drop_na(s_cd40l) 
#Clean cols
chal2 <- chal2[,!str_ends(colnames(chal2), ".y")] #clean cols
colnames(chal2) <- str_replace_all(colnames(chal2), ".x", "")
chal2 <- chal2[,-14]

table(chal2$time_point)


#Fold change #
chal2 <- 
  chal2 %>%
  group_by(lab_id) %>%
  mutate(
    FC_D0 = ((s_cd40l[time_point == "D0"] - s_cd40l[time_point == "D-28"])/
               s_cd40l[time_point =="D-28"]),
    FC_D1 = ((s_cd40l[time_point == "D1"] - s_cd40l[time_point == "D-28"])/
               s_cd40l[time_point =="D-28"]))

#CD40L D0 cors ##

chal2 %>%
  filter(time_point == "D0") %>%
  ggplot(aes(x=FC_D1, y = adnob)) +
  geom_point(color="#9b42f5", size=3.6, alpha=0.7) +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  ylab("ADNOB D28\n")+
  xlab("sCD40L FC D28")+
  ggtitle("Post-vaccination")+
  white_cor




#loop cor at day 0 (28 days post-vac)
chal2$il1b <- as.numeric(chal2$il1b)
colnames <- colnames(chal2[,5:13]) #Just select cytokine
Y <- chal2 %>% filter(time_point == "D0") %>% select(adnob)
Y = Y$adnob

for(i in colnames){
  plt <- chal2 %>% filter(time_point == "D0") %>%
    ggplot(aes_string(x=i, y = Y)) +
    geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
    stat_cor(method = "spearman")+
    geom_hline(yintercept=0, size=0.06, color="black") + 
    geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
    ylab ("ADNOB")+
    xlab(gsub("(\\_|\\-)", " ", i))+
    ggtitle("28 days post-vaccination")+
    white_cor
  print(plt)
  Sys.sleep(2)
}

#GRO-a FC ####
#Fold change #
chal2 <- 
  chal2 %>%
  group_by(lab_id) %>%
  mutate(
    gFC_D28 = ((gr_oa[time_point == "D0"] - gr_oa[time_point == "D-28"])/
               gr_oa[time_point =="D-28"]),
    gFC_D1 = ((gr_oa[time_point == "D-27"] - gr_oa[time_point == "D-28"])/
               gr_oa[time_point =="D-28"]),
    gFC_D3 = ((gr_oa[time_point == "D-25"] - gr_oa[time_point == "D-28"])/
               gr_oa[time_point =="D-28"]))

colnames <- colnames(chal2[,str_detect(colnames(chal2), "gFC")]) #Just select cytokine
Y <- chal2 %>% filter(time_point == "D0") %>% select(adnob)
Y = Y$adnob

for(i in colnames){
  plt <- chal2 %>% filter(time_point == "D0") %>%
    ggplot(aes_string(x=i, y = Y)) +
    geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
    stat_cor(method = "spearman")+
    geom_hline(yintercept=0, size=0.06, color="black") + 
    geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
    ylab ("ADNOB")+
    xlab(gsub("(\\_|\\-)", " ", i))+
    ggtitle("28 days post-vaccination")+
    white_cor
  print(plt)
  Sys.sleep(2)
}


#Neut migration calcs ####

mig <-
  pre %>%
  group_by(lab_id) %>%
  mutate(
    gcsf_FC_D1 = ((g_csf[time_point == "-27"] - g_csf[time_point == "PV"])/
                g_csf[time_point =="PV"]),
    gcsf_FC_D3 = ((g_csf[time_point == "-25"] - g_csf[time_point == "PV"])/
                g_csf[time_point =="PV"]),
    IL17_FC_D1 = ((il_17a[time_point == "-27"] - il_17a[time_point == "PV"])/
                    il_17a[time_point =="PV"]),
    IL17_FC_D3 = ((il_17a[time_point == "-25"] - il_17a[time_point == "PV"])/
                   il_17a[time_point =="PV"]))
  
str(mig)
colnames <- colnames(mig[,str_detect(colnames(mig), "_FC")]) #Just select cytokine
Y <- mig %>% filter(time_point == "PV") %>% select(n_m_ratio)
Y = Y$n_m_ratio

for(i in colnames){
  plt <- mig %>% filter(time_point == "PV") %>%
    ggplot(aes_string(x=i, y = Y)) +
    geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
    stat_cor(method = "spearman")+
    geom_hline(yintercept=0, size=0.06, color="black") + 
    geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
    ylab ("Neutrophil Count")+
    xlab(gsub("(\\_|\\-)", " ", i))+
    ggtitle("28 days post-vaccination")+
    white_cor
  print(plt)
  Sys.sleep(2)
}


mig %>% filter(time_point == "PV") %>%
  ggplot(aes(x=il_17a, y = neutrophils10_9_l)) +
  geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  ylab ("Neutrophil Count")+
  xlab(gsub("(\\_|\\-)", " ", i))+
  ggtitle("28 days post-vaccination")+
  white_cor

mig %>% filter(time_point == "-25") %>%
  ggplot(aes(x=g_csf, y = neutrophils10_9_l)) +
  geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  ylab ("Neutrophil Count")+
  xlab(gsub("(\\_|\\-)", " ", i))+
  ggtitle("28 days post-vaccination")+
  white_cor


#Loop cor challenge cytokines w #time to diagnosis ####

#chal2

chal2 <- 
  chal2 %>%
  group_by(lab_id) %>%
  mutate(
    gFC_D1 = ((gr_oa[time_point == "D1"] - gr_oa[time_point == "D-28"])/
                gr_oa[time_point =="D-28"]),
    gFC_12h = ((gr_oa[time_point == "12h"] - gr_oa[time_point == "D-28"])/
                gr_oa[time_point =="D-28"]),
    gFC_D0 = ((gr_oa[time_point == "D0"] - gr_oa[time_point == "D-28"])/
                 gr_oa[time_point =="D-28"]))
    


  xvar = "cdFC_D1"
  yvar = "time_to_positive_bc_hours"
  name = "CXCL1"
  time = "12h"

chal2 %>% filter(time_point == "D0")%>%
  ggplot(aes_string(x=xvar, y =yvar))+
  geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  ylab ("time to positive bc (hrs)")+
    xlab(paste("Fold change", name, time))+
  ggtitle("Post-challenge")+
  xlim(-1,3)+
  white_cor



chal2 <- 
  chal2 %>%
  group_by(lab_id) %>%
  mutate(
    cdFC_D1 = ((s_cd40l[time_point == "D1"] - s_cd40l[time_point == "D0"])/
                s_cd40l[time_point =="D0"]),
    cdFC_12h = ((s_cd40l[time_point == "12h"] - s_cd40l[time_point == "D0"])/
                 s_cd40l[time_point =="D0"]),
    cdFC_D0 = ((s_cd40l[time_point == "D0"] - s_cd40l[time_point == "D-28"])/
                 s_cd40l[time_point =="D-28"]))


#Double check ADNOB cors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



xvar = "gFC_D0"
yvar = "adnob"
name = "CXCL1"
time = "D28"

chal2 %>% filter(time_point == "D0")%>%
  ggplot(aes_string(x=xvar, y =yvar))+
  geom_point(color="#9b42f5", size=3.6, alpha=0.8) +
  stat_cor(method = "spearman")+
  geom_hline(yintercept=0, size=0.06, color="black") + 
  geom_smooth(method=lm, alpha=0.25, color="black", fill="black")+
  ylab ("ADNOB D28")+
  xlab(paste("Fold change", name, time))+
  ggtitle("Post-Vaccination")+
  xlim(-1,1)+
  white_cor

#no correlations
