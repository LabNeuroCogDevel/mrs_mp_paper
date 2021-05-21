library(glue)

# column and row are counted in opposite directions between MRRC and LNCD
# matrix is 216x216: 1 => 216, 216 => 1
reverse_index <- function(x, len=216) len+1-x

#### Get data and remove bad quality data ####
read_mrsi_roi <- function() {
   MRS_all <-
    read.csv("data/13MP20200207_LCMv2fixidx.csv") %>%
    mutate(
      x216=x, y216=y, # keep origian for comparison
      x=reverse_index(x), y=reverse_index(y) # but rewrite x and y to be what we expect
    )
}

# Step 1 Outlier Detection - visual inspection of LCModel fits/spectra
# create a list of who to remove and remove them
remove_bad_qc <- function(MRS_all){
   lcm_qc_fail <-
       read.table("data/lcm_bad_visual_qc.txt",header=FALSE) %>%
       separate(V1, c("ld8", "junk","y","x"),extra="merge", sep = "[-.]") %>%
       select(-junk) %>%
       mutate(bad=TRUE)
   # anti_join: any that are in lcm_qc should be removed 
   MRS_qc <- MRS_all %>% 
      merge(lcm_qc_fail, by=c("ld8", "x", "y"), all=T)  %>%
      filter(is.na(bad)) %>%
      select(-bad)
}

# want only visit1 and need to data (roi and Glu.Cr)
visit1_and_data <- function(MRS_qc) {
   MRS <- MRS_qc %>%
       filter(
        #keep only visit 1 people
        visitnum==1,
        # be very sure visit 2 isn't leakign in (id passes filter?)
        ! ld8 %in% c("10195_20191205"), # 20210520WF - visitnum for 10195_2019 is correct now, already excluded
        #keep people's correct coordinates
        !is.na(roi),
        # get rid of junk data noticed recently 
        Glu.Cr != 0)
}

add_inv_and_quad_age <- function(MRS) {
  #make inverse age column
  MRS$invage <- 1/MRS$age
  #make age^2 column
  MRS$age2 <- (MRS$age - mean(MRS$age))^2
  return(MRS)
}

# Step 2 Outlier Detection
# - get rid of peole who have bad data for 3 major metabolite peaks - GPC+Cho, NAA+NAAG, Cr
# - get rid of people who have lots of macromolecule in their spectra, as that can create distortions
# would keep NAs around for later (but there are none)
remove_major_met_outliers <- function(MRS) {
   MRS %>%
    filter(GPC.Cho.SD  <= 10 | is.na(GPC.Cho.SD),
           NAA.NAAG.SD <= 10 | is.na(NAA.NAAG.SD),
           Cr.SD       <= 10 | is.na(Cr.SD),
           MM20.Cr     <=  3 | is.na(MM20.Cr))
}


# change y column (e.g. "GABA.Cr") to be the residuals after linear model with invage (or a provided agecol)
mk_age_resid <- function(d, y, agecol='invage') {
  d[,y] <- lm(as.formula(glue("{y} ~ {agecol}")), d)$residuals
  return(d)
}
