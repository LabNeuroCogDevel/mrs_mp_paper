#### Code for Paper ####
# 5/17/21 

library(tidyverse)
library(ggplot2)
library(readxl)
library(Hmisc)
library(lmerTest)
library(corrplot)
library(RColorBrewer)
library(data.table)
library(missMDA)
library(FactoMineR)
library(tidyr)

#### Get data and remove bad quality data ####
MRS_all <- read.csv("data/13MP20200207_LCMv2fixidx.csv")

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
lcm_qa <- read.table("data/lcm_bad_visual_qc.txt",header=FALSE) %>%
       separate(lcm, "V1", c("ld8", "junk","y","x"),extra="merge", sep = "[-.]") %>%
       select(lcm, -junk) %>%
       mutate(bad=TRUE)
MRS <- MRS %>% mutate(x=216+1-x,y=216+1-y) %>%
   merge(lcm_qa, by=c("ld8", "x", "y"), all=T)  %>%
   filter(MRS, is.na(bad)) %>%
   select(MRS, -bad)

#keep only visit 1 people
MRS <- MRS %>% filter(visitnum==1)
#keep people's correct coordinates
MRS <- MRS %>% filter(!is.na(roi))
#get rid of people who are actually visit 2 but for some reason aren't filtered out
MRS <- MRS %>% filter(ld8!="10195_20191205")

# get rid of junk data noticed recently 
MRS<- MRS %>% filter(Glu.Cr != 0)

# save out a dataframe to share data after this step

# Step 2 Outlier Detection - get rid of peole who have bad data for 3 major metabolite peaks - GPC+Cho, NAA+NAAG, Cr
MRS<- filter(MRS, GPC.Cho.SD <= 10 | is.na(GPC.Cho.SD))
MRS <- filter(MRS, NAA.NAAG.SD <= 10 | is.na(NAA.NAAG.SD))
MRS <- filter(MRS, Cr.SD <= 10 | is.na(Cr.SD))

# Step 3 Outlier Detection - get rid of people who have lots of macromolecule in their spectra, as that can create distortions
MRS <- filter(MRS, MM20.Cr <= 3 | is.na(MM20.Cr))

#make inverse age column
MRS$invage <- 1/MRS$age
#make age^2 column
MRS$age2 <- (MRS$age - mean(MRS$age))^2

z_thres = 2

#### Participant histogram ####
# doing this for the full sample prior to data quality exclusions
# pick one ROI as representative 
ROI1 <- MRS %>% filter(roi==1) # n = 144
agehist <- ggplot(ROI1, aes(x=age,fill=sex)) + geom_histogram(color="black", binwidth = 1) + theme_classic(base_size = 15) + xlab("Age (years)") + ylab("Count") + labs(fill = "Sex")

ggsave(agehist, filename = "0519_agehistogram.pdf", width = 7, height = 6, units = "in", dpi = 300)

#### Glutamate and Age ####

# Create dataframe with good quality Glutamate data 
MRS_glu <- MRS %>% filter(Glu.SD <=20)

# ROI 1 (R Anterior Insula) and 2 (L Anterior Insula)

ROI12_Glu <- MRS_glu %>% filter(roi == 1 | roi == 2)
ROI12_Glu <- ROI12_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))
ROI12_Glu <- ROI12_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

#ROI12_Glu %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI12_Glu_age <- lmer(data=ROI12_Glu, Glu.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_age)
ROI12_Glu_invage <- lmer(data=ROI12_Glu, Glu.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_invage) 
ROI12_Glu_quadage <- lmer(data=ROI12_Glu, Glu.Cr ~ age + age2 + label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_quadage) 

AIC(ROI12_Glu_age)
AIC(ROI12_Glu_invage) # best fit by a lot
AIC(ROI12_Glu_quadage)

# effect sizes
ROI12_zGlu <- lmer(data=ROI12_Glu, zscore_glu ~ zscore_invage + label + (1|ld8))
summary(ROI12_zGlu) # just glu ~ age effect

ROI12_zGlu_invage <- lmer(data=ROI12_Glu, zscore_glu ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI12_zGlu_invage) # all covariates for table
summary(ROI12_zGlu_invage)$coefficients

# test for interactions w/ age
ROI12_Glu_hemi_int <- lmer(data=ROI12_Glu, Glu.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI12_Glu_hemi_int) 
ROI12_Glu_sex_int <- lmer(data=ROI12_Glu, Glu.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI12_Glu_sex_int) 
ROI12_Glu_gmrat_int <- lmer(data=ROI12_Glu, Glu.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI12_Glu_gmrat_int) 

ggplot(ROI12_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
aiglusd_age <- lmer(data=ROI12_Glu, scale(Glu.SD) ~ scale(age) +(1|ld8))
summary(aiglusd_age)


# ROI 7 (ACC)
ROI7_Glu <- MRS_glu %>% filter(roi == 7)

ROI7_Glu <- ROI7_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI7_Glu <- ROI7_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

ROI7_Glu %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI7_Glu_age <- lm(data=ROI7_Glu, Glu.Cr ~ age + sex + GMrat)
summary(ROI7_Glu_age)
ROI7_Glu_invage <- lm(data=ROI7_Glu, Glu.Cr ~ invage + sex + GMrat)
summary(ROI7_Glu_invage)
ROI7_Glu_quadage <- lm(data=ROI7_Glu, Glu.Cr ~ age + age2 + sex + GMrat)
summary(ROI7_Glu_quadage)

AIC(ROI7_Glu_age) 
AIC(ROI7_Glu_invage)
AIC(ROI7_Glu_quadage)

#effect sizes
ROI7_zGlu <- lm(data=ROI7_Glu, zscore_glu ~ zscore_invage)
summary(ROI7_zGlu)
ROI7_zGlu_invage <- lm(data=ROI7_Glu, zscore_glu ~ zscore_invage + sex + zscore_gm)
summary(ROI7_zGlu_invage)

# test for interactions
ROI7_age_sex_int <- lm(data=ROI7_Glu, Glu.Cr ~ invage * sex + GMrat)
summary(ROI7_age_sex_int)
ROI7_gm_sex_int <- lm(data=ROI7_Glu, Glu.Cr ~ invage * GMrat + sex)
summary(ROI7_gm_sex_int)

ggplot(ROI7_Glu, aes(x=age, y=Glu.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ggplot(ROI7_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
accglusd_age <- lm(data=ROI7_Glu, scale(Glu.SD) ~ scale(age))
summary(accglusd_age)

# ROI 8 (MPFC)
ROI8_Glu <- MRS_glu %>% filter( roi == 8)

ROI8_Glu <- ROI8_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI8_Glu <- ROI8_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)

ROI8_Glu %>% select(ld8, label) %>% spread(label, -ld8)

ROI8_Glu_age <- lm(data=ROI8_Glu, Glu.Cr ~ age + sex + GMrat)
summary(ROI8_Glu_age)
ROI8_Glu_invage <- lm(data=ROI8_Glu, Glu.Cr ~ invage + sex + GMrat)
summary(ROI8_Glu_invage)
ROI8_Glu_quadage <- lm(data=ROI8_Glu, Glu.Cr ~ age + age2 + sex + GMrat)
summary(ROI8_Glu_quadage)

AIC(ROI8_Glu_age)
AIC(ROI8_Glu_invage)
AIC(ROI8_Glu_quadage)

ROI8_zGlu <- lm(data=ROI8_Glu, zscore_glu ~ zscore_invage)
summary(ROI8_zGlu)
ROI8_zGlu_invage <- lm(data=ROI8_Glu, zscore_glu ~ zscore_invage + sex + zscore_gm)
summary(ROI8_zGlu_invage)

# test for interactions
ROI8_age_sex_int <- lm(data=ROI8_Glu, Glu.Cr ~ invage * sex + GMrat)
summary(ROI8_age_sex_int)
ROI8_gm_sex_int <- lm(data=ROI8_Glu, Glu.Cr ~ invage * GMrat + sex)
summary(ROI8_gm_sex_int)

ggplot(ROI8_Glu, aes(x=age, y=Glu.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ggplot(ROI8_Glu, aes(x=age, y=Glu.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
mpfcglusd_age <- lm(data=ROI8_Glu, scale(Glu.SD) ~ scale(age))
summary(mpfcglusd_age)


# ROI 9 ( R DLPFC) and 10 (L DLPFC)

ROI910_Glu <- MRS_glu %>% filter(roi == 9 | roi == 10)
ROI910_Glu <- ROI910_Glu%>%
  mutate(zscore_glu=scale(Glu.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI910_Glu <- ROI910_Glu %>% 
  filter(abs(zscore_glu) <= z_thres)


ROI910_Glu %>% select(ld8,label) %>% gather(label)
ROI910_Glu %>% select(ld8, label) %>% spread(label, -ld8)

ROI910_Glu_age <- lmer(data=ROI910_Glu, Glu.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI910_Glu_age) 

ROI910_Glu_invage <- lmer(data=ROI910_Glu, Glu.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI910_Glu_invage) 

ROI910_Glu_quadage <- lmer(data=ROI910_Glu, Glu.Cr ~ age + age2+ label + sex + GMrat + (1|ld8))
summary(ROI910_Glu_quadage) 

AIC(ROI910_Glu_age)
AIC(ROI910_Glu_invage)
AIC(ROI910_Glu_quadage)

ROI910_Glu_invagez <- lmer(data=ROI910_Glu, zscore_glu ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI910_Glu_invagez) 
summary(ROI910_Glu_invagez)$coefficients

ROI910_zGlu <- lmer(data=ROI910_Glu, zscore_glu ~ zscore_invage + label + (1|ld8))
summary(ROI910_zGlu)

# test for interactions w/ age
ROI910_age_hemi_int <- lmer(data=ROI910_Glu, Glu.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI910_age_hemi_int) 
ROI910_age_sex_int <- lmer(data=ROI910_Glu, Glu.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI910_age_sex_int) 
ROI910_age_gmrat_int <- lmer(data=ROI910_Glu, Glu.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI910_age_gmrat_int) 


ggplot(ROI910_Glu, aes(x=age, y=Glu.SD, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()
glusd_age <- lm(data=ROI910_Glu, scale(Glu.SD) ~ scale(age))
summary(glusd_age)




#### GABA and Age ####

# Create dataframe with good quality GABA data 
MRS_GABA <- MRS %>% filter(GABA.SD <=20)
z_thres = 2

# ROI 1 (R Anterior Insula) and 2 (L Anterior Insula)

ROI12_GABA <- MRS_GABA %>% filter(roi == 1 | roi == 2)
ROI12_GABA <- ROI12_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI12_GABA <- ROI12_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

ROI12_GABA %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI12_GABA_age <- lmer(data=ROI12_GABA, GABA.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_age)
ROI12_GABA_invage <- lmer(data=ROI12_GABA, GABA.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_invage) 
ROI12_GABA_quadage <- lmer(data=ROI12_GABA, GABA.Cr ~ age + age2 + label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_quadage) 

AIC(ROI12_GABA_age)
AIC(ROI12_GABA_invage)
AIC(ROI12_GABA_quadage)

# effect sizes
ROI12_zGABA <- lmer(data=ROI12_GABA, zscore_GABA ~ zscore_invage + label + (1|ld8))
summary(ROI12_zGABA) # just GABA ~ age effect

ROI12_zGABA_invage <- lmer(data=ROI12_GABA, zscore_GABA ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI12_zGABA_invage) # all covariates for table
summary(ROI12_zGABA_invage)$coefficients

# test for interactions w/ age
ROI12_GABA_hemi_int <- lmer(data=ROI12_GABA, GABA.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI12_GABA_hemi_int) 
ROI12_GABA_hemi_int <- lmer(data=ROI12_GABA, zscore_GABA ~ zscore_invage * label + sex + zscore_gm + (1|ld8))
summary(ROI12_GABA_hemi_int) 
summary(ROI12_GABA_hemi_int)$coefficients
ROI12_GABA_sex_int <- lmer(data=ROI12_GABA, GABA.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI12_GABA_sex_int) 
ROI12_GABA_gmrat_int <- lmer(data=ROI12_GABA, GABA.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI12_GABA_gmrat_int) 


ggplot(ROI12_GABA, aes(x=age, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
aigabasd_age <- lmer(data=ROI12_GABA, scale(GABA.SD) ~ scale(age) +(1|ld8))
summary(aigabasd_age)

# ROI 7 (ACC)
ROI7_GABA <- MRS_GABA %>% filter(roi == 7)
ROI7_GABA <- ROI7_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI7_GABA <- ROI7_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

ROI7_GABA %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI7_GABA_age <- lm(data=ROI7_GABA, GABA.Cr ~ age + sex + GMrat)
summary(ROI7_GABA_age)
ROI7_GABA_invage <- lm(data=ROI7_GABA, GABA.Cr ~ invage + sex + GMrat)
summary(ROI7_GABA_invage)
ROI7_GABA_quadage <- lm(data=ROI7_GABA, GABA.Cr ~ age + age2 + sex + GMrat)
summary(ROI7_GABA_quadage)

AIC(ROI7_GABA_age) 
AIC(ROI7_GABA_invage) # best fit
AIC(ROI7_GABA_quadage)

#effect sizes
ROI7_zGABA <- lm(data=ROI7_GABA, zscore_GABA ~ zscore_invage)
summary(ROI7_zGABA)
ROI7_zGABA_invage <- lm(data=ROI7_GABA, zscore_GABA ~ zscore_invage + sex + zscore_gm)
summary(ROI7_zGABA_invage)

# test for interactions
ROI7_age_sex_int <- lm(data=ROI7_GABA, GABA.Cr ~ invage * sex + GMrat)
summary(ROI7_age_sex_int)
ROI7_gm_sex_int <- lm(data=ROI7_GABA, GABA.Cr ~ invage * GMrat + sex)
summary(ROI7_gm_sex_int)


ggplot(ROI7_GABA, aes(x=age, y=GABA.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()

ggplot(ROI7_GABA, aes(x=age, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
accgabasd_age <- lm(data=ROI7_GABA, scale(GABA.SD) ~ scale(age))
summary(accgabasd_age)


# ROI 8 (MPFC)
ROI8_GABA <- MRS_GABA %>% filter( roi == 8)
ROI8_GABA <- ROI8_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))

ROI8_GABA <- ROI8_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

ROI8_GABA %>% select(ld8, label) %>% spread(label, -ld8)

ROI8_GABA_age <- lm(data=ROI8_GABA, GABA.Cr ~ age + sex + GMrat)
summary(ROI8_GABA_age)
ROI8_GABA_invage <- lm(data=ROI8_GABA, GABA.Cr ~ invage + sex + GMrat)
summary(ROI8_GABA_invage)
ROI8_GABA_quadage <- lm(data=ROI8_GABA, GABA.Cr ~ age + age2 + sex + GMrat)
summary(ROI8_GABA_quadage)

AIC(ROI8_GABA_age)
AIC(ROI8_GABA_invage) # best fit
AIC(ROI8_GABA_quadage)

ROI8_zGABA <- lm(data=ROI8_GABA, zscore_GABA ~ zscore_invage)
summary(ROI8_zGABA)
ROI8_zGABA_invage <- lm(data=ROI8_GABA, zscore_GABA ~ zscore_invage + sex + zscore_gm)
summary(ROI8_zGABA_invage)

# test for interactions
ROI8_age_sex_int <- lm(data=ROI8_GABA, GABA.Cr ~ invage * sex + GMrat)
summary(ROI8_age_sex_int)
ROI8_gm_sex_int <- lm(data=ROI8_GABA, GABA.Cr ~ invage * GMrat + sex)
summary(ROI8_gm_sex_int)

ggplot(ROI8_GABA, aes(x=age, y=GABA.Cr)) + geom_point() + geom_smooth(method="lm") + theme_classic()
ggplot(ROI8_GABA, aes(x=age, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
mpfcgabasd_age <- lm(data=ROI8_GABA, scale(GABA.SD) ~ scale(age))
summary(mpfcgabasd_age)
ggplot(ROI8_GABA, aes(x=GABA.Cr, y=GABA.SD)) + geom_point() + geom_smooth(method="lm") + theme_classic()
mpfcgabacr_age <- lm(data=ROI8_GABA, scale(GABA.SD) ~ scale(GABA.Cr))
summary(mpfcgabacr_age)




# ROI 9 ( R DLPFC) and 10 (L DLPFC)
ROI910_GABA <- MRS_GABA %>% filter(roi == 9 | roi == 10)
ROI910_GABA <- ROI910_GABA%>%
  mutate(zscore_GABA=scale(GABA.Cr), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI910_GABA <- ROI910_GABA %>% 
  filter(abs(zscore_GABA) <= z_thres)

ROI910_GABA %>% select(ld8, label) %>% spread(label, -ld8)


ROI910_GABA_invage <- lmer(data=ROI910_GABA, GABA.Cr ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI910_GABA_invage) 
ROI910_GABA_age <- lmer(data=ROI910_GABA, GABA.Cr ~ age + label + sex + GMrat + (1|ld8))
summary(ROI910_GABA_age) 
ROI910_GABA_quadage <- lmer(data=ROI910_GABA, GABA.Cr ~ age + age2+ label + sex + GMrat + (1|ld8))
summary(ROI910_GABA_quadage) 


AIC(ROI910_GABA_age)
AIC(ROI910_GABA_invage)
AIC(ROI910_GABA_quadage)

ROI910_zGABA <- lmer(data=ROI910_GABA, zscore_GABA ~ zscore_invage + label + (1|ld8))
summary(ROI910_zGABA)

ROI910_GABA_invagez <- lmer(data=ROI910_GABA, zscore_GABA ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI910_GABA_invagez) 
summary(ROI910_GABA_invagez)$coefficients



# test for interactions w/ age
ROI910_age_hemi_int <- lmer(data=ROI910_GABA, GABA.Cr ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI910_age_hemi_int) 
ROI910_age_sex_int <- lmer(data=ROI910_GABA, GABA.Cr ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI910_age_sex_int) 
ROI910_age_gmrat_int <- lmer(data=ROI910_GABA, GABA.Cr ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI910_age_gmrat_int) 


ggplot(ROI910_GABA, aes(x=zscore_invage, y=GABA.Cr, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()
ggplot(ROI910_GABA, aes(x=age, y=GABA.SD, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()
gabasd_age <- lm(data=ROI910_GABA, scale(GABA.SD) ~ scale(age))
summary(gabasd_age)

#### Ratio and Age ####
# create dataframe keeping only people who have both good quality Glu and good GABA data to
# make a ratio out of 
MRS_Ratio <- MRS_glu %>% filter(GABA.SD <=20)

MRS_Ratio$Ratio <- MRS_Ratio$Glu.Cr/MRS_Ratio$GABA.Cr


# ROI 1 (R Anterior Insula) and 2 (L Anterior Insula)

ROI12_Ratio <- MRS_Ratio %>% filter(roi == 1 | roi == 2)
ROI12_Ratio <- ROI12_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI12_Ratio <- ROI12_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

ROI12_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI12_Ratio_age <- lmer(data=ROI12_Ratio, Ratio ~ age + label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_age)
ROI12_Ratio_invage <- lmer(data=ROI12_Ratio, Ratio ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_invage) 
ROI12_Ratio_quadage <- lmer(data=ROI12_Ratio, Ratio ~ age + age2 + label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_quadage) 

AIC(ROI12_Ratio_age)
AIC(ROI12_Ratio_invage)
AIC(ROI12_Ratio_quadage)

# effect sizes
ROI12_zRatio <- lmer(data=ROI12_Ratio, zscore_Ratio ~ zscore_invage + label + (1|ld8))
summary(ROI12_zRatio) # just Ratio ~ age effect

ROI12_zRatio_invage <- lmer(data=ROI12_Ratio, zscore_Ratio ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI12_zRatio_invage) # all covariates for table

# test for interactions w/ age
ROI12_Ratio_hemi_int <- lmer(data=ROI12_Ratio, Ratio ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI12_Ratio_hemi_int) 
ROI12_Ratio_sex_int <- lmer(data=ROI12_Ratio, Ratio ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI12_Ratio_sex_int) 
ROI12_Ratio_gmrat_int <- lmer(data=ROI12_Ratio, Ratio ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI12_Ratio_gmrat_int) 

ROI1 <- ROI12_Ratio %>% filter(roi==1)
ROI1_zRatio <- lm(data=ROI1, zscore_Ratio ~ zscore_invage + sex + GMrat)
summary(ROI1_zRatio) # just Ratio ~ age effect

ROI2 <- ROI12_Ratio %>% filter(roi==2)
ROI2_zRatio <- lm(data=ROI2, zscore_Ratio ~ zscore_invage)
summary(ROI2_zRatio) # just Ratio ~ age effect

# ROI 7 (ACC)
ROI7_Ratio <- MRS_Ratio %>% filter(roi == 7)
ROI7_Ratio <- ROI7_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI7_Ratio <- ROI7_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

ROI7_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

# pick best model
ROI7_Ratio_age <- lm(data=ROI7_Ratio, Ratio ~ age + sex + GMrat)
summary(ROI7_Ratio_age)
ROI7_Ratio_invage <- lm(data=ROI7_Ratio, Ratio ~ invage + sex + GMrat)
summary(ROI7_Ratio_invage)
ROI7_Ratio_quadage <- lm(data=ROI7_Ratio, Ratio ~ age + age2 + sex + GMrat)
summary(ROI7_Ratio_quadage)

AIC(ROI7_Ratio_age) 
AIC(ROI7_Ratio_invage) # best fit
AIC(ROI7_Ratio_quadage)

#effect sizes
ROI7_zRatio <- lm(data=ROI7_Ratio, zscore_Ratio ~ zscore_invage)
summary(ROI7_zRatio)
ROI7_zRatio_invage <- lm(data=ROI7_Ratio, zscore_Ratio ~ zscore_invage + sex + zscore_gm)
summary(ROI7_zRatio_invage)

# test for interactions
ROI7_age_sex_int <- lm(data=ROI7_Ratio, Ratio ~ invage * sex + GMrat)
summary(ROI7_age_sex_int)
ROI7_gm_sex_int <- lm(data=ROI7_Ratio, Ratio ~ invage * GMrat + sex)
summary(ROI7_gm_sex_int)


ggplot(ROI7_Ratio, aes(x=age, y=Ratio)) + geom_point() + geom_smooth(method="lm") + theme_classic()


# ROI 8 (MPFC)
ROI8_Ratio <- MRS_Ratio %>% filter( roi == 8)
ROI8_Ratio <- ROI8_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI8_Ratio <- ROI8_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

ROI8_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

ROI8_Ratio_age <- lm(data=ROI8_Ratio, Ratio ~ age + sex + GMrat)
summary(ROI8_Ratio_age)
ROI8_Ratio_invage <- lm(data=ROI8_Ratio, Ratio ~ invage + sex + GMrat)
summary(ROI8_Ratio_invage)
ROI8_Ratio_quadage <- lm(data=ROI8_Ratio, Ratio ~ age + age2 + sex + GMrat)
summary(ROI8_Ratio_quadage)

AIC(ROI8_Ratio_age)
AIC(ROI8_Ratio_invage) # best fit
AIC(ROI8_Ratio_quadage)

ROI8_zRatio <- lm(data=ROI8_Ratio, zscore_Ratio ~ zscore_invage)
summary(ROI8_zRatio)
ROI8_zRatio_invage <- lm(data=ROI8_Ratio, zscore_Ratio ~ zscore_invage + sex + zscore_gm)
summary(ROI8_zRatio_invage)

# test for interactions
ROI8_age_sex_int <- lm(data=ROI8_Ratio, Ratio ~ invage * sex + GMrat)
summary(ROI8_age_sex_int)
ROI8_gm_sex_int <- lm(data=ROI8_Ratio, Ratio ~ invage * GMrat + sex)
summary(ROI8_gm_sex_int)

ggplot(ROI8_Ratio, aes(x=age, y=Ratio)) + geom_point() + geom_smooth(method="lm") + theme_classic()


# ROI 9 ( R DLPFC) and 10 (L DLPFC)
ROI910_Ratio <- MRS_Ratio %>% filter(roi == 9 | roi == 10)

#ROI910_Ratio %>% select(ld8,label) %>% gather(label)
ROI910_Ratio <- ROI910_Ratio%>%
  mutate(zscore_Ratio=scale(Ratio), zscore_invage=scale(invage), zscore_gm=scale(GMrat))


ROI910_Ratio <- ROI910_Ratio %>% 
  filter(abs(zscore_Ratio) <= z_thres)

ROI910_Ratio %>% select(ld8, label) %>% spread(label, -ld8)

ROI910_Ratio_invage <- lmer(data=ROI910_Ratio, Ratio ~ invage + label + sex + GMrat + (1|ld8))
summary(ROI910_Ratio_invage) 
ROI910_Ratio_age <- lmer(data=ROI910_Ratio, Ratio ~ age + label + sex + GMrat + (1|ld8))
summary(ROI910_Ratio_age) 
ROI910_Ratio_quadage <- lmer(data=ROI910_Ratio, Ratio ~ age + age2+ label + sex + GMrat + (1|ld8))
summary(ROI910_Ratio_quadage) 
AIC(ROI910_Ratio_age)
AIC(ROI910_Ratio_invage)
AIC(ROI910_Ratio_quadage)


ROI910_zRatio <- lmer(data=ROI910_Ratio, zscore_Ratio ~ zscore_invage + label + (1|ld8))
summary(ROI910_zRatio)
ROI910_zRatio_invage <- lmer(data=ROI910_Ratio, zscore_Ratio ~ zscore_invage + label + sex + zscore_gm + (1|ld8))
summary(ROI910_zRatio_invage) 

# test for interactions w/ age
ROI910_age_hemi_int <- lmer(data=ROI910_Ratio, Ratio ~ invage * label + sex + GMrat + (1|ld8))
summary(ROI910_age_hemi_int) 
ROI910_age_sex_int <- lmer(data=ROI910_Ratio, Ratio ~ invage * sex + label + GMrat + (1|ld8))
summary(ROI910_age_sex_int) 
ROI910_age_gmrat_int <- lmer(data=ROI910_Ratio, Ratio ~ invage * GMrat + sex + label + (1|ld8))
summary(ROI910_age_gmrat_int) 


ggplot(ROI910_Ratio, aes(x=age, y=Ratio, color=label)) + geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + theme_classic()

#### Plot GABA and Glu on same graph ####
gaba12 <- ROI12_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat, label) %>% mutate(metabolite="GABA")
glu12 <- ROI12_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat, label) %>% mutate(metabolite="Glu")
gabaglu12 <- rbind(gaba12,glu12)
AI_gabaglu <- ggplot(gabaglu12) + aes(y=concentration, x=age, color=metabolite, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  ggtitle("Anterior Insula") +
  labs(color= "Metabolite", shape = "Hemisphere") + 
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(AI_gabaglu, filename = "0518_AI_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)


gaba910 <- ROI910_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat, label) %>% mutate(metabolite="GABA")
glu910 <- ROI910_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat, label) %>% mutate(metabolite="Glu")
gabaglu910 <- rbind(gaba910,glu910)
DLPFC_gabaglu <- ggplot(gabaglu910) + aes(y=concentration, x=age, color=metabolite, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  labs(color= "Metabolite", shape = "Hemisphere") + 
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67))) +
  ggtitle("Dorsolateral Prefrontal Cortex")
ggsave(DLPFC_gabaglu, filename = "0518_DLPFC_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)


gaba7 <- ROI7_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat) %>% mutate(metabolite="GABA")
glu7 <- ROI7_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat) %>% mutate(metabolite="Glu")
gabaglu7 <- rbind(gaba7,glu7)
ACC_gabaglu <- ggplot(gabaglu7) + aes(y=concentration, x=age, color=metabolite) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  labs(color= "Metabolite") + 
  ggtitle("ACC")+ 
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(ACC_gabaglu, filename = "0518_ACC_Glu_GABA.pdf", width = 6, height = 6, units = "in", dpi = 300)

gaba8 <- ROI8_GABA %>% select(concentration=GABA.Cr, age, sex, GMrat) %>% mutate(metabolite="GABA")
glu8 <- ROI8_Glu %>% select(concentration=Glu.Cr, age, sex, GMrat) %>% mutate(metabolite="Glu")
gabaglu8 <- rbind(gaba8,glu8)
MPFC_gabaglu <- ggplot(gabaglu8) + aes(y=concentration, x=age, color=metabolite) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Concentration (Metabolite/Cr)") + 
  labs(color= "Metabolite") + 
  ggtitle("MPFC") +
  scale_shape_manual(values=c(16, 1)) + 
  scale_color_manual(values=c("#FF6666", "#66B2FF"), guide = guide_legend(reverse = TRUE)) +
  guides(linetype=FALSE) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(MPFC_gabaglu, filename = "0518_MPFC_Glu_GABA.png", width = 5, height = 5, units = "in", dpi = 300)



#### Plot Ratio and Age ####
AI_ratio <- ggplot(ROI12_Ratio) + aes(y=Ratio, x=age, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Concentration") + 
  labs(shape = "Hemisphere") + 
  ggtitle("Anterior Insula") +
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(AI_ratio, filename = "0518_AI_ratio.pdf", width = 6, height = 6, units = "in", dpi = 300)

DLPFC_ratio <- ggplot(ROI910_Ratio) + aes(y=Ratio, x=age, linetype = label, shape=label) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Concentration") + 
  labs(shape = "Hemisphere") + 
  scale_shape_manual(values=c(16, 1)) + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(DLPFC_ratio, filename = "0518_DLPFC_Ratio.png", width = 5, height = 5, units = "in", dpi = 300)


ACC_ratio <- ggplot(ROI7_Ratio) + aes(y=Ratio, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Concentration") + 
  ggtitle("ACC") + 
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(ACC_ratio, filename = "0518_ACC_Ratio.png", width = 5, height = 5, units = "in", dpi = 300)

MPFC_ratio <- ggplot(ROI8_Ratio) + aes(y=Ratio, x=age) +
  geom_point() + geom_smooth(method="lm", formula = y~I(1/x), fullrange=T) + 
  theme_classic(base_size = 15) +xlab("Age (years)") + ylab("Glu/GABA Concentration") + 
  ggtitle("MPFC")+
  theme(legend.key = element_rect(fill = "white", colour = "black"), legend.title = element_text(size=rel(0.67)), legend.text = element_text(size=rel(0.67)))
ggsave(MPFC_ratio, filename = "0518_MPFC_Ratio.png", width = 5, height = 5, units = "in", dpi = 300)


#### Plot correlations as bar graphs ####
# make dataframe that has only good glutamate and gaba data for all ROIs
MRS_glu <- MRS_glu %>%
  group_by(roi) %>%
  mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 
MRS_corr <- MRS_glu %>%
  filter(GABA.SD <= 20) %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup %>%
  mutate(agegrp = cut(age,
                      breaks=c(0,16,22,Inf),
                      labels=c("10-16","17-22", "23-30")))


# comparing correlations using cocor package 
region_look <- list("1"="RAntInsula",
                    "2"="LAntInsula",
                    "7"='ACC',
                    "8"='MPFC',
                    "9"="RDLPFC",
                    "10"="LDLPFC")
keep_rois <- c(1,2,7,8,9,10)   # or keep_rois <- names(region_look)

region_avg <- MRS_corr %>%
  filter(roi %in% keep_rois) %>%          # remove any rois not in region_loopup
  select(ld8, agegrp,Glu.Cr,GABA.Cr, roi) %>%
  mutate(region=unlist(region_look[as.character(roi)])) %>%
  filter(!is.na(region)) %>%
  group_by(region,ld8, agegrp) %>%
  summarise_at(vars(Glu.Cr, GABA.Cr), mean)

cor_vals <- region_avg %>%
  group_by(region,agegrp) %>%
  summarise(GabaGlu_r=cor(Glu.Cr, GABA.Cr), 
            GabaGlu_p=cor.test(Glu.Cr, GABA.Cr, method=c("pearson"), use = "complete.obs")$p.value, 
            GabaGlu_lb=cor.test(Glu.Cr, GABA.Cr, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[1], 
            GabaGlu_ub=cor.test(Glu.Cr, GABA.Cr, method=c("pearson"), use = "complete.obs", conf.level=.95)$conf.int[2], 
            n=n())

corr_bar <- ggplot(cor_vals) +
  aes(x=agegrp, y=GabaGlu_r, fill=agegrp) +
  geom_bar(stat="identity") +
  facet_wrap(~region) +
  theme_classic(base_size = 15) + labs(x='Age',y='Glu GABA Corr (r)', fill = "Age Group") +
  scale_fill_brewer(palette = "viridis")

ggsave(corr_bar, filename = "0518_corr_bargraph.pdf", width = 8, height = 6, units = "in", dpi = 300)


#### Plot correlations as line graphs ####
MRS1 <- MRS_corr %>% filter(roi==1)
MRS2 <- MRS_corr %>% filter(roi==2)
MRS7 <- MRS_corr %>% filter(roi==7)
MRS8 <- MRS_corr %>% filter(roi==8)
MRS9 <- MRS_corr %>% filter(roi==9)
MRS10 <- MRS_corr %>% filter(roi==10)

MRS1 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS1 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS1 <- MRS1 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))
mrs1<- MRS1 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
ggsave(mrs1, filename = "0518_MRS1corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS2 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS2 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS2 <- MRS2 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))
mrs2<- MRS2 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
ggsave(mrs2, filename = "0518_MRS2corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS7 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS7 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS7<- MRS7 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))


mrs7<- MRS7 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
ggsave(mrs7, filename = "0518_MRS7corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS8 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS8 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS8<- MRS8 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))


ggplot(MRS8) + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +xlab("Glu.Cr") + ylab("GABA.Cr") + 
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")

mrs8<- MRS8 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
ggsave(mrs8, filename = "0518_MRS8corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


MRS9 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS9 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS9<- MRS9 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))

mrs9<- MRS9 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
ggsave(mrs9, filename = "0518_MRS9corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)

MRS10 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% split(.$agegrp) %>% sapply(function(x) cor(x$Glu.Cr, x$GABA.Cr))
MRS10 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age)) %>% group_by(agegrp) %>% tally
MRS10<- MRS10 %>% mutate(agegrp=cut(breaks=c(0,16,22,Inf), age))
mrs10<- MRS10 %>% 
  mutate(agegrp = fct_recode(agegrp, "10-16" = "(0,16]","17-22" = "(16,22]", "23-30" = "(22,Inf]")) %>% 
  ggplot() + aes(y=GABA.Cr, x=Glu.Cr) +
  geom_point() + geom_smooth(method="lm") + 
  theme_classic(base_size = 18) +
  theme(legend.key = element_rect(fill = "white", colour = "black")) + 
  facet_wrap(~ agegrp) +
  xlab("Glu/Cr") +
  ylab("GABA/Cr")
ggsave(mrs10, filename = "0518_MRS10corr_scatter.pdf", width = 9, height = 5, units = "in", dpi = 300)


#### Correlation matrices ####

# used cleaned data 
MRS_glu <- MRS_glu %>%
  group_by(roi) %>%
  mutate(zscore=scale(Glu.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 

MRS_GABA <- MRS_GABA %>%
  group_by(roi) %>%
  mutate(zscore=scale(GABA.Cr, center=T, scale=T)) %>%
  filter(abs(zscore) <= z_thres) %>% ungroup 

# How correlated is Glutamate between regions?

# Without age residualized out 
# Step 1 - take good Glu data in MRS_glu dataframe and change the format
MRS_glu_wide <- pivot_wider(MRS_glu, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide <- select(MRS_glu_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide <- MRS_glu_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

#ggsave(corrmat1, filename = "0519_GluCorrMat.pdf", width = 6, height = 6, units = "in", dpi = 300)
# issue with the ggsave; need to resave these out as high res figs; says it cant save a matrix

# With age residualized out 

# Step 1 - Residualize age out of the correlation
mk_age_resid_glu <- function(d) {d[,'Glu.Cr'] <- lm(Glu.Cr ~ invage, d)$residuals; return(d) }
resids_glu <- MRS_glu %>% split(MRS_glu$label) %>% lapply(mk_age_resid_glu)
resids_glu <- resids_glu %>% bind_rows
# Step 1 - take good Glu data in MRS_glu dataframe and change the format
MRS_glu_wide_resids <- pivot_wider(resids_glu, id_cols=ld8, names_from=label, values_from=Glu.Cr)
MRS_glu_wide_resids <- select(MRS_glu_wide_resids, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_glu_wide_resids <- MRS_glu_wide_resids %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
Glu_corr<- cor(MRS_glu_wide_resids %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_glu_wide_resids %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(Glu_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))


# How correlated is GABA between regions?

# without age residualized out 
MRS_GABA_wide <- pivot_wider(MRS_GABA, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide <- select(MRS_GABA_wide, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide <- MRS_GABA_wide %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))

# with age residualized out 
# Step 1 - take good GABA data in MRS_GABA dataframe and change the format
mk_age_resid_GABA <- function(d) {d[,'GABA.Cr'] <- lm(GABA.Cr ~ invage, d)$residuals; return(d) }
resids_GABA <- MRS_GABA %>% split(MRS_GABA$label) %>% lapply(mk_age_resid_GABA)
resids_GABA <- resids_GABA %>% bind_rows

MRS_GABA_wide_resids <- pivot_wider(resids_GABA, id_cols=ld8, names_from=label, values_from=GABA.Cr)
MRS_GABA_wide_resids <- select(MRS_GABA_wide_resids, -`R STS`, -`L STS`, -`R Thalamus`, -`L Caudate`, -`R Caudate`, -`L Posterior Insula`, -`R Posterior Insula`)
MRS_GABA_wide_resids <- MRS_GABA_wide_resids %>% relocate('ld8', 'ACC', 'MPFC', 'R Anterior Insula', 'L Anterior Insula', 'R DLPFC', 'L DLPFC')

#Step 2 - Make correlation matrix
GABA_corr<- cor(MRS_GABA_wide_resids %>% select_if(is.numeric), use="pairwise.complete.obs")

#Step 3 - Make significance matrix 
pvals <- cor.mtest(MRS_GABA_wide_resids %>% select_if(is.numeric), conf.level = .95)

#Step 4 - plot correlation matrix 
corrplot(GABA_corr, p.mat = pvals$p, sig.level = .05, method = "number", type="lower", tl.col = "black", tl.srt = 45, insig = "blank", col = rev(brewer.pal(n = 8, name = "Spectral")))




