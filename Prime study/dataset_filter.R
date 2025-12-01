library(dplyr)
library(rstanarm)
library(bayesplot)
library(ggplot2)
library(tidyverse)

setwd("~/Desktop/Thesis/code/prime/AllProvidedFiles_309/PDS_DSA_20050203")

#' ===============================================
#' Import and Select
#' ===============================================================

#' Event data frame
#' -----------------
kras<-haven::read_sas("biomark_pds2019.sas7bdat")
kras0<-kras |> select("SUBJID", "BMMTR1")

adsl<-haven::read_sas("adsl_pds2019.sas7bdat")
event.df0<-adsl |>
  left_join(kras0, by="SUBJID") |>
  filter(BMMTR1=="Wild-type") |>
  mutate(STUDY="1", EVENTYR=DTHDY/365.25, EVENTFL=DTH) |>
  select("SUBJID", "STUDY", "ATRT", "EVENTYR", "EVENTFL")

UID.event0<-unique(event.df0$SUBJID)
#' length(UID.event0)
#' 514


#' Biomarker (SLD) data frame
#' -----------------
adtr<-haven::read_sas("adls_pds2019.sas7bdat")
biom.df0<-adtr |>
  filter(LSCAT=="Target lesion", !is.na(LSSLD)) |>
  mutate(STUDY="1", BIOMYR=VISITDY/365.25, BIOMVAL=LSSLD) |>
  group_by(SUBJID, VISITDY) |> slice(1) |> ungroup() |>
  select("SUBJID", "BIOMYR", "BIOMVAL")

UID.biom0<-unique(biom.df0$SUBJID)
#' length(UID.biom0)
#' 488


#' Retain matching patients
#' -----------------
retainID<-intersect(UID.event0, UID.biom0)
#' 263

event.df<-event.df0 |> filter(SUBJID %in% retainID)
#' length(unique(event.df$SUBJID))
#' 263

desn0<-event.df0 |>
  filter(SUBJID %in% retainID) |>
  select(SUBJID, STUDY, ATRT)

biom.df<-biom.df0 |>
  filter(SUBJID %in% retainID) |>
  left_join(desn0, by="SUBJID")

print(paste("Number of subjects:", length(unique(biom.df$SUBJID))))
#------------------------------------------------------------------------------------

#' Display SLD spaghetti
ybreaks<-c(3, 30, 100, 300)
mycols<-c(rev(ghibli::ghibli_palettes$YesterdayMedium)[c(2,4)])
ggplot(data=biom.df, aes(x=BIOMYR, y=BIOMVAL))+
  geom_point(colour="grey33", alpha=0.3, size=0.9)+
  geom_line(aes(group=SUBJID, colour=as.factor(ATRT)), alpha=0.6)+
  facet_wrap(~ATRT)+
  scale_x_continuous("Year", breaks=0.5*(0:5))+
  scale_y_continuous("SLD (mm)", breaks=ybreaks)+
  scale_colour_manual(values=mycols, guide="none")+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())

# ------------------ filtering --------------------------------------------

# the treatment group to analyse
target_trt <- "Panitumumab + FOLFOX"

# Calculate the time of the first observation for each subject
first_obs_times <- biom.df %>%
  group_by(SUBJID) %>%
  summarise(first_obs_yr = min(BIOMYR))

# Filter to keep subjects whose first observation is not more than one week before time 0
a_week_in_years <- 7 / 365.25
subjects_within_week <- first_obs_times %>%
  filter(first_obs_yr >= -a_week_in_years & first_obs_yr <= 0)

# Now add the treatment filter BEFORE taking top 45
top_45_subjects <- biom.df %>%
  filter(
    ATRT == target_trt,
    SUBJID %in% subjects_within_week$SUBJID
  ) %>%
  count(SUBJID, sort = TRUE) %>%
  slice_head(n = 45)     # if fewer than 45 exist, you'll just get that smaller number

# Final filtered dataset: restrict to those subjects AND the treatment
biom.df.filtered <- biom.df %>%
  filter(
    SUBJID %in% top_45_subjects$SUBJID,
    ATRT == target_trt
  )

# Summary checks
cat("--- Data Filtering and Selection Summary ---\n")
cat("Treatment:", target_trt, "\n")
cat("Number of subjects after filtering:", length(unique(biom.df.filtered$SUBJID)), "\n")
cat("Total number of observations:", nrow(biom.df.filtered), "\n\n")

cat("Number of observations per treatment group (should be 1 level):\n")
print(table(biom.df.filtered$ATRT))
cat("------------------------------------------\n")

# Save the filtered dataset
setwd("~/Desktop/Thesis/code/prime")
saveRDS(biom.df.filtered, file = "biom.df.filtered.rds")
# Or CSV:
# write.csv(biom.df.filtered, "biom.df.filtered.panitumumab_top45.csv", row.names = FALSE)

# (Optional) Plot just this treatment group
ggplot(biom.df.filtered, aes(x = BIOMYR, y = BIOMVAL, group = SUBJID)) +
  geom_point(colour = "grey33", alpha = 0.35, size = 0.9) +
  geom_line(alpha = 0.6, colour = "#718a52") +
  scale_x_continuous("Year", breaks = 0.5 * (0:5)) +
  scale_y_continuous("SLD (mm)", breaks = c(3, 30, 100, 300)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("SLD Over Time: Panitumumab + FOLFOX (Top 45 by row count)")

