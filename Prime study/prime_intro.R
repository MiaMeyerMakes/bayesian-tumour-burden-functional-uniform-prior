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