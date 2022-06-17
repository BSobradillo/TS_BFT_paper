
##======================================================##
## SCRIPT : BFT TS PUBLICATION 2020
##
## Authors : Jon Uranga,  Guillermo Boyra & Bea Sobradillo
## Last update : 17-06-2022
## Description: 
##              Este script lee todas las exportaciones de TS de los lances 
##              registrados en la campaña de BFT Index y calcula el TS tipico del bluefin tuna 
##              
##======================================================##


## Load packages  -------------------------------------------------------------------------
library(data.table)
library(stringr)
library(BuoysAZTI) # Own package with necesary functions and tables to run the script
library(dplyr)
library(RPostgreSQL)
library(tidyr)
library(seewave)
library(ggplot2)
source("getmode.R")
library(ggridges)
library(geosphere)


# BIO: LOAD BIOLOGICAL SAMPLING----
muestreo<-fread("Datos/muestreo_2016_bft_index.csv")
muestreo <- muestreo[-6,]
muestreo_lon <- muestreo %>% 
  pivot_longer(8:30 ,values_to = "FL") %>% 
  group_by(calera) %>% 
  summarise(meanFL = mean(FL, na.rm = T),
            sdFL = sd(FL, na.rm = T),
            sbL = meanFL/5,
            RcTD = (21^2)/3.9,
            RcSB = ((sbL/2)^2)/3.9)


# AC: LOAD AND RBIND ALL ST FROM ALL TRACKS AND SETS----
idx.files <- list.files("Datos/exports/",pattern="(targets)");idx.files
idx.files<-as.data.frame(idx.files)
names(idx.files)<-"source"
idx.files$calera<-as.numeric(sub("\\D+","",idx.files$source))
idx.files <- idx.files %>% arrange(calera) 

all_ts<-c()
for(i in 1:dim(idx.files)[1]){
  # i<-1
  ft <- setDT(read.table(file=paste("Datos/exports/",idx.files$source[i],sep=""), sep=",",header = TRUE))
  
  ft <- ft %>% 
    mutate(Target_true_depth = Target_range + 3) %>% 
    group_by(Region_name) %>% 
    mutate(N_ST = n()) %>% 
    ungroup() %>% 
    filter(N_ST>2)#filter tracks from set1 with less than 3 targets
    
  # idx <- which(fish_tracks$id==fish_tracks_filter[i])
  
  id <- i
  Num_Fish_Tracks <- length(unique(ft$Region_name))
  FL <- muestreo$MEDIA[i]
  N_Sampled_Fish <- length(which(!is.na(muestreo[i,8:30])))
  ft<-cbind(ft, id, FL, N_Sampled_Fish, Num_Fish_Tracks)
  all_ts <- rbind(all_ts,ft)
  print(i)
  
}

setDT(all_ts)
table(all_ts$id)

# all_ts <- all_ts[-which(all_ts$id==6)]
# all_ts <- all_ts[-which(all_ts$id==9)]


all_ts <- all_ts %>% 
  mutate(track_id = paste(id,ï..Region_ID, sep="-"))
# number of tracks detected: 
length(unique(all_ts$track_id))

# Angle calculation (Toby) ----

# all_ts2<-as.data.frame(all_ts)
all_ts <- all_ts %>% 
  group_by(track_id) %>% 
  mutate(lat2 = lag(Target_latitude),
         lon2 = lag(Target_longitude),
         range2 = lag(Target_true_depth)) %>% 
  ungroup()

all_ts<- all_ts[which(!is.na(all_ts$lon2)),]

all_ts %>% 
  select(id, Target_longitude) %>% 
  group_by(id) %>% 
  summarise(min(Target_longitude), max(Target_longitude))

# analizamos sets
ggplot(all_ts, aes(x=Target_true_depth, y=TS_comp)) + 
  geom_point()+
  facet_wrap(.~id)


# number of total tracks
length(unique(all_ts$track_id))


# filtramos set sin GPS y nearfield
all_ts <- all_ts %>%
  filter(id!=7, Target_range>4.55) 
  # mutate(FL = case_when(id==6~58.8,TRUE~FL))
unique(all_ts$id)
length(unique(all_ts$track_id))

# analizamos sets
ggplot(all_ts, aes(x=Target_true_depth, y=TS_comp)) + 
  geom_point()+
  facet_wrap(.~id)



setDT(all_ts)
all_ts$x<-NA
all_ts$y<-NA
for(i in 1 : dim(all_ts)[1]){
  all_ts$x[i] <-as.numeric(distm(c(all_ts$Target_longitude[i], 
                                   all_ts$Target_latitude[i]),
                                 c(all_ts$lon2[i],
                                   all_ts$lat2[i]), fun = distHaversine))
  all_ts$y[i] <-as.numeric(all_ts$range2[i]-all_ts$Target_true_depth[i])
  print(i)
}

all_ts <- all_ts %>% 
  mutate(Tilt_angle = atan(y/x)*180/pi,Tilt_angle_abs = abs(Tilt_angle))

summary(all_ts)
all_ts<-all_ts[which(!is.na(all_ts$Tilt_angle)),]

# ggplot(all_ts,aes(x=Tilt_angle,y=TS_comp))+
  # geom_point()+
  # theme_light()  



####SUMMARY TABLES ----
names(all_ts)
unique(all_ts$id)
table(all_ts$track_id)
summary(all_ts$N_ST)


# SINGLE TARGETS/TRACK----
# 1. Average of single targets per track
bio_ft_summary <- all_ts %>% 
  mutate(FL = as.numeric(as.character(FL))) %>%
  group_by(track_id) %>%
  # group_by(id) %>% 
  summarise(medianTS = round(median(TS_comp),1),
            meanTS = round(meandB(TS_comp),1),
            N_ST = n(),
            FL = round(mean(FL),1),
            N_length = length(FL),
            sdev_TS = round(sddB(TS_comp ),1),
            Mean_depth = round(mean(Target_true_depth),0), 
            sdev_depth = round(sd(Target_true_depth),1),
            Num_Fish_Tracks = unique(Num_Fish_Tracks),
            N_Sampled_Fish = unique(N_Sampled_Fish),
            
            N_depth = length(Target_true_depth),
            tilt = mean(Tilt_angle),
            id = unique(id)) %>%
  ungroup() 




# DATA MINING: remove haul 7 because it does not have GPS data and will not be useful to tell the tilt relationship with TS. Ranges from TD below 12 are also removed (2xNF)
#el set 6 se pesca en el campo cercano por lo que lo descartamos
bio_ft_summary <-  bio_ft_summary %>%
# filter(Mean_depth > 9.5 ) %>%
mutate(id = case_when(id == 8 ~ 7,
id == 9 ~ 8,
id == 10 ~ 9, TRUE ~ as.numeric(id))) #para que no haya saltos en el numero de pesca

unique(bio_ft_summary$id)

 
 

bio_ft_summary_order <- bio_ft_summary %>% 
  arrange(FL) 

# 2. Model: average ST/track vs length
modelo <- lm(bio_ft_summary$meanTS ~ log10(bio_ft_summary$FL))
summary(modelo)
modelo$coefficients
TS.L.cor <- cor(bio_ft_summary$meanTS, log10(bio_ft_summary$FL)); TS.L.cor
# fwrite(bio_ft_summary,"Datos/Resumen_lances.csv")

# 3. Ponderamos los valores promedio de los tracks al numero de ST por track
# media = suma(sigma_ponderado)/suma(peso de ponderacion)
ponderate_data <- bio_ft_summary %>% 
  mutate(
    sigma = 10^(meanTS/10), 
    pon_weight = Num_Fish_Tracks/sum(Num_Fish_Tracks),
    sigma_pon = sigma*pon_weight,
    length_pon = FL*pon_weight) %>% 
  summarise(
    mean_sigma = sum(sigma_pon)/sum(pon_weight),
    TS_pon = 10*log10(mean_sigma),
    sdTS = sddB(TS_pon),
    length_pon = sum(length_pon)/sum(pon_weight),
    # sd_TS = sddB(meanTS),
    sd_length = sd(FL), 
    min_L = min(FL),
    max_L = max(FL),
    mean_z = mean(Mean_depth),
    sd_z = mean(sdev_depth))


# 4. FIGURE TS(ST)-L: free fitting

a<-setDT(summary(modelo)[4])[2,4]; a
my_text <- paste("Corr. coeff.:" ,
                 round(as.numeric(TS.L.cor),2), 
                 "\np-value: ",
                 round(as.numeric(a),7))
my_text2 <- paste("- Slope (a):" ,
                  round(as.numeric(setDT(summary(modelo)[4])[2,1]),2), 
                  "\nIntercept (b): ",
                  round(as.numeric(setDT(summary(modelo)[4])[1,1]),2))

ggplot(bio_ft_summary,aes(x=FL,y=meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=meanTS-(sdev_TS/2), ymax=meanTS+(sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  theme_minimal()+
  xlab("Log-Length (cm)")+
  ylab("TS (dB)")+
  annotate("text", x=60, y=-15, label=my_text, hjust=0) +
  annotate("text", x=60, y=-19, label=my_text2, hjust=0) +
  scale_x_log10() +
  theme_classic() 

# ggsave("Figuras/ts_log_length_ST.tiff", width=17, height=12, units="cm")


# FIGURE TS(ST)-L----

SE <- sd(bio_ft_summary$meanTS)/sqrt(length(bio_ft_summary$meanTS))
L_mean_slope  <-  20
flN  <-  transform(bio_ft_summary,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(meanTS) ~ 1 + offset(log10(flN$FL)*L_mean_slope), data = flN)
m1$coefficients
summary(m1)
b20cor <- cor(as.numeric(flN$meanTS),1 + offset(log10(flN$FL)*L_mean_slope)); b20cor
# slope chosen by lm(...)
m0 <- lm(as.numeric(meanTS) ~ log10(FL), data=flN)
m0$coefficients
summary(m0)


# b20 linear model
# my_text_b20 <- paste("Corr. coeff.:" ,
#    round(as.numeric(b20cor),2),
#    "\n- - Slope (a):" ,
# 20, 
# "\nIntercept (b): ",
# round(as.numeric(setDT(summary(m1)[4])[1,1]),2))
my_text_b20 <- paste("- - Slope (a):" ,
                     20, 
                     "\nIntercept (b): ",
                     round(as.numeric(setDT(summary(m1)[4])[1,1]),2))

ggplot(bio_ft_summary,aes(x=FL,y=meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  # geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=meanTS-(sdev_TS/2), ymax=meanTS+(sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  geom_abline(slope=20, intercept=coef(m1)[1], lty=2) +
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1]) +
  theme_minimal()+
  xlab("Length (cm)")+
  ylab("TS (dB)")+
  annotate("text", x=60, y=-10, label=my_text2, hjust=0) +
  annotate("text", x=60, y=-18, label=my_text_b20, hjust=0) +
  
  scale_x_log10() +
  theme_classic() 


ggsave("Figuras/TSmined_NFplus40mREALDEPTH-st-logL-b20.tiff", width=17, height=12, units="cm")

## TRACKS/SET----
# 1. Average of tracks per set 
tt_summary <- bio_ft_summary %>% 
  mutate(length = as.numeric(as.character(FL))) %>%
  group_by(id) %>% 
  summarise(tt_meanTS = round(meandB(meanTS),1),
            tt_N_TS = n(),
            tt_length = round(mean(FL),1),
            tt_N_length = length(tt_length),
            tt_sdev_TS = round(sddB(meanTS),1),
            tt_Mean_depth = round(mean(Mean_depth),0),
            tt_sdev_depth = round(sd(Mean_depth),1),
            tt_N_depth = length(tt_Mean_depth), 
            tt_Num_Fish_Tracks = unique(Num_Fish_Tracks),
            tt_N_Sampled_Fish = unique(N_Sampled_Fish)) %>%
  ungroup()

# 2. Model: average tracks/set vs length---
tt_modelo <- lm(tt_summary$tt_meanTS ~ log10(tt_summary$tt_length))
summary(tt_modelo)
tt_modelo$coefficients
tt_TS.L.cor <- cor(tt_summary$tt_meanTS, log10(tt_summary$tt_length)); tt_TS.L.cor
# fwrite(tt_summary,"Datos/Resumen_lances.csv")

tt_summary_order <- tt_summary %>% arrange(tt_length) ; tt_summary_order



# 3. Ponderamos los valores promedio de cada set al numero de tracks por set
# media = suma(sigma_ponderado)/suma(peso de ponderacion)
tt_ponderate_data <- tt_summary %>% 
  mutate(
    sigma = 10^(tt_meanTS/10), 
    pon_weight = tt_Num_Fish_Tracks/sum(tt_Num_Fish_Tracks),
    sigma_pon = sigma*pon_weight,
    length_pon = FL*pon_weight) %>% 
  summarise(
    mean_sigma = sum(sigma_pon)/sum(pon_weight),
    tt_TS_pon = 10*log10(mean_sigma),
    tt_sdTS = mean(tt_sdev_TS),
    length_pon = sum(length_pon)/sum(pon_weight),
    # sd_TS = sddB(meanTS),
    sd_length = sd(FL), 
    min_L = min(FL),
    max_L = max(FL),
    mean_z = mean(tt_Mean_depth),
    sd_z = mean(tt_sdev_depth))



# 4. FIGURE TS(tracks)-L: free fitting
a<-setDT(summary(tt_modelo)[4])[2,4]; a
tt_my_text <- paste("Corr. coeff.:" ,
                 round(as.numeric(tt_TS.L.cor),2), 
                 "\np-value: ",
                 round(as.numeric(a),7))
tt_my_text2 <- paste("- Slope (a):" ,
                  round(as.numeric(setDT(summary(tt_modelo)[4])[2,1]),2), 
                  "\nIntercept (b): ",
                  round(as.numeric(setDT(summary(tt_modelo)[4])[1,1]),2))
ggplot(tt_summary,aes(x=tt_length,y=tt_meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=tt_meanTS-(tt_sdev_TS/2), ymax=tt_meanTS+(tt_sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  theme_minimal()+
  xlab("Log-Length (cm)")+
  ylab("TS (dB)")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=60, y=-15, label=tt_my_text, hjust=0) +
  annotate("text", x=60, y=-19, label=tt_my_text2, hjust=0) +
  scale_x_log10() +
  theme_classic()


# 5. FIGURE TS(Tracks)-L: fixed slope

SE <- sd(tt_summary$tt_meanTS)/sqrt(length(tt_summary$tt_meanTS))
L_mean_slope  <-  20
flN  <-  transform(tt_summary,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(tt_meanTS) ~ 1 + offset(log10(flN$tt_length)*L_mean_slope), data = flN)
m1$coefficients
summary(m1)
tt_b20cor <- cor(as.numeric(flN$tt_meanTS),1 + offset(log10(flN$tt_length)*L_mean_slope)); tt_b20cor
# slope chosen by lm(...)
m0 <- lm(as.numeric(tt_meanTS) ~ log10(tt_length), data=flN)
m0$coefficients
summary(m0)


# b20 linear model
# my_text_b20 <- paste("Corr. coeff.:" ,
#    round(as.numeric(b20cor),2),
#    "\n- - Slope (a):" ,
# 20, 
# "\nIntercept (b): ",
# round(as.numeric(setDT(summary(m1)[4])[1,1]),2))
tt_my_text_b20 <- paste("- - Slope (a):" ,
                     20, 
                     "\nIntercept (b): ",
                     round(as.numeric(setDT(summary(m1)[4])[1,1]),2))

ggplot(tt_summary,aes(x=tt_length,y=tt_meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  # geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=tt_meanTS-(tt_sdev_TS/2), ymax=tt_meanTS+(tt_sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  geom_abline(slope=20, intercept=coef(m1)[1], lty=2) +
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1]) +
  theme_minimal()+
  xlab("Length (cm)")+
  ylab("TS (dB)")+
  annotate("text", x=60, y=-15, label=tt_my_text2, hjust=0) +
  annotate("text", x=60, y=-20, label=tt_my_text_b20, hjust=0) +
  
  scale_x_log10() +
  theme_classic() 


ggsave("Figuras/TSmined-tracks-logL-b20.tiff", width=17, height=12, units="cm")






# TS histogram facets----

ggplot(bio_ft_summary, aes(x=meanTS)) +
  geom_histogram(aes(y = stat(density)), fill="white",binwidth = 2,color="black") +
  geom_density(alpha = 0.2, fill = "white",color="black", size=0.75) +
  theme(legend.position = "none") +
  # geom_vline(aes(xintercept = -45),color = "#FC4E07",linetype = "dashed", size = 0.1)+
  # ggtitle ("TS freq. distr. for all sets") + 
  xlab ("TS (dB)") + 
  ylab ("Frequency") +
  theme(text = element_text(size = 15)) +
  facet_wrap(~id,nrow=3,scales = "free_y") +
  theme_classic()

ggsave("Figuras/ggridges_tracks_mined_plus40mREALDEPTH_distrib_ALL_facet.tiff", width = 27, height=12,units="cm", dpi = 300)
       # width=27, height=17, units="cm")





# FIGURA RIDGES DEL TS MEDIO DE LOS FT EN MUESTREOS BIOLOGICOS----
# CON ETIQUETAS DE NUM DE FT

theme_set(theme_minimal())


ggplot(all_ts, aes(x = TS_comp, y = as.factor(FL), fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.975)
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  )


# Los ridges de los sets con la misma talla los junto para evitar la mala visibilidad del set con 1 solo track

ggplot(data=bio_ft_summary, aes(x = meanTS, y = as.factor(round(FL,2)))) +
  geom_density_ridges_gradient(panel_scaling = T) +
  # scale_fill_viridis_c(name="N fish tracks") +
  theme_ridges(center_axis_labels = T) + 
  theme(legend.position = "none") +
  # geom_vline(aes(xintercept = -45),color = "#FC4E07",linetype = "dashed", size = 0.1)+ 
  # ggtitle ("TS freq. distr. and fish lengths for all sets ") + 
  xlab ("Target Strength (dB)") +
  ylab ("Mean fish length (cm)")  +
  # theme(axis.title = element_text(size = 12)) +
  geom_text(aes(x = -5, y=0.64,label="FT(N)"),color="grey30",)+
  # geom_text(aes(x = -1, y=0.64,label="Set"),color="blue")+
  geom_text(aes(x = -5, y=1.3,label=tt_summary_order$tt_Num_Fish_Tracks[1]),color="grey30")+
  geom_text(aes(x = -5, y=2.3,label=tt_summary_order$tt_Num_Fish_Tracks[2]+1),color="grey30")+
  geom_text(aes(x = -5, y=3.3,label=tt_summary_order$tt_Num_Fish_Tracks[4]),color="grey30")+
  geom_text(aes(x = -5, y=4.3,label=tt_summary_order$tt_Num_Fish_Tracks[5]),color="grey30")+
  geom_text(aes(x = -5, y=5.3,label=tt_summary_order$tt_Num_Fish_Tracks[6]),color="grey30")+
  geom_text(aes(x = -5, y=6.3,label=tt_summary_order$tt_Num_Fish_Tracks[7]),color="grey30")+
  geom_text(aes(x = -5, y=7.3,label=tt_summary_order$tt_Num_Fish_Tracks[8]),color="grey30")+
  geom_text(aes(x = -5, y=8.3,label=tt_summary_order$tt_Num_Fish_Tracks[9]),color="grey30")
  # geom_text(aes(x = -5, y=9.3,label=tt_summary_order$tt_Num_Fish_Tracks[]),color="grey30")
  # geom_text(aes(x = -5, y=10.3,label=bio_ft_summary_order$Num_Fish_Tracks[10]),color="grey30")+
  # geom_text(aes(x = -5, y=11.3,label=bio_ft_summary_order$Num_Fish_Tracks[11]),color="grey30")
  # geom_text(aes(x = -1, y=1.3,label=bio_ft_summary_order$id[1]),color="blue")+
  # geom_text(aes(x = -1, y=2.3,label=bio_ft_summary_order$id[2]),color="blue")+
  # geom_text(aes(x = -1, y=3.3,label=bio_ft_summary_order$id[3]),color="blue")+
  # geom_text(aes(x = -1, y=4.3,label=bio_ft_summary_order$id[4]),color="blue")+
  # geom_text(aes(x = -1, y=5.3,label=bio_ft_summary_order$id[5]),color="blue")+
  # geom_text(aes(x = -1, y=6.3,label=bio_ft_summary_order$id[6]),color="blue")+
  # geom_text(aes(x = -1, y=7.3,label=bio_ft_summary_order$id[7]),color="blue")+
  # geom_text(aes(x = -1, y=8.3,label=bio_ft_summary_order$id[8]),color="blue")+
  # geom_text(aes(x = -1, y=9.3,label=bio_ft_summary_order$id[9]),color="blue")+
  # geom_text(aes(x = -1, y=10.3,label=bio_ft_summary_order$id[10]),color="blue")

ggsave("figuras/ggridges_mined_plus40.tiff", width=15, height=15, units="cm")
# ggsave("Figuras/ggridges_por_BIO_ken6.tiff", width=17, height=17, units="cm")



# FIGURA  TS ~ TILT ----
library(magick)
# img1 <- image_read("Figuras/fish_tilt_img/tuna_edit1.PNG")
img2 <- image_read("Figuras/fish_tilt_img/tuna_edit2.PNG")
img3 <- image_read("Figuras/fish_tilt_img/tuna_edit3.PNG")
# img4 <- image_read("Figuras/fish_tilt_img/tuna_edit4.PNG")
img5 <- image_read("Figuras/fish_tilt_img/tuna_edit5.PNG")
sb_tilt <- image_read("Figuras/fish_tilt_img/sb_tilt.PNG")

# img_rot1 <- image_rotate(img,70)
# img_rot2 <- image_rotate(img,15)
# img_rot3 <- image_rotate(img,0)
# img_rot4 <- image_rotate(img,-70)
# img_rot5 <- image_rotate(img,-30)
# print(img_rot)

model_TS_alpha <- loess(TS_comp~Tilt_angle, all_ts)
summary(model_TS_alpha)

  ggplot(all_ts, aes(x=(Tilt_angle),y=TS_comp)) +
  geom_point(size=1.3, alpha=0.2) +
  geom_smooth(method = "loess",span = 0.8) +
  # geom_smooth(method = "gam") +
  theme_classic() +
  xlab("Tilt angle (degrees)") +
  ylab("TS (dB)") +
    # theme(axis.line = element_line(arrow = arrow(type='closed', length = unit(10,'pt')))) +
  # annotation_raster(img1, ymin = -35,ymax= -25,xmin = -90,xmax = -55) + 
    annotation_raster(img2, ymin = -30,ymax= -23,xmin = -50,xmax = -15) + 
    annotation_raster(img3, ymin = -35,ymax= -30,xmin = -15,xmax = 15) + 
    # annotation_raster(img4, ymin = -40,ymax= -30,xmin = 20,xmax = 55) +
    annotation_raster(img5, ymin = -40,ymax= -31,xmin = 20,xmax = 53) +
    annotation_raster(sb_tilt, ymin = -17,ymax= -10,xmin = -95,xmax = -45)
    
  ggsave("Figuras/TS-mined-tilt-tuna.tiff", width = 17, height=12,units="cm", dpi = 300)

  
  #  FIGURA TILT HISTOGRAM----
  
  ggplot(all_ts, aes(x=Tilt_angle)) +
    geom_histogram(aes(y = stat(density)), fill="white",binwidth = 2,color="black") +
    geom_density(alpha = 0.2, fill = "white",color="black", size=0.75) +
    theme(legend.position = "none") +
    xlab ("Tilt angle (deg)") + 
    ylab ("Frequency") +
    theme(text = element_text(size = 12)) +
    theme_classic()
  
  ggsave("Figuras/tilt_histogram.tiff", width = 17, height=12,units="cm", dpi = 300)
  

  # FIGURA TS ~ L + TILT----

  modelo_TS_L_tilt <- lm(bio_ft_summary$meanTS ~ log10(bio_ft_summary$FL) + bio_ft_summary$tilt))
  summary(modelo_TS_L_tilt)

# ggsave("ts_tilt_FT_loess.tiff", width=18, height=18, units="cm")
# ggsave("figuras/ts_tilt_FT_gam.tiff", width=21, height=17, units="cm")
# ggsave("ts_tilt_FT_gam_sin_085.tiff", width=18, height=18, units="cm")





# FIGURA  LENGTH ~ DEPTH----

# Modelo optimo:
modelo_L_Z <- lm(bio_ft_summary$FL ~ bio_ft_summary$Mean_depth)
  summary(modelo_L_Z)
a<-setDT(summary(modelo_L_Z)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo_L_Z)[8]),3),
                 "\np-value: ",
                 round(as.numeric(a),3))



ggplot(bio_ft_summary,aes(x=Mean_depth,y=as.numeric(FL)))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  # geom_errorbar(aes(ymin=meanTS-sdev, ymax=meanTS+sdev), width=.0051,col="grey",alpha=0.49) + 
  # scale_y_continuous(breaks=c(0,-20,-40,-60), labels=c("0","20","40","60")) +
  theme_minimal()+
  ylab("Length (cm)")+
  xlab("Mean depth (m)")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=40, y=120, label=my_text, hjust=0) +
  theme_classic()
  # theme(axis.line = element_line(arrow = arrow(type='closed', length = unit(10,'pt')))) 

# ggsave("depth_length.tiff", width=18, height=18, units="cm")
ggsave("Figuras/depth_length_ken6.tiff", width=17, height=12, units="cm")




# FIGURA  TS ~ DEPTH----

aa <- bio_ft_summary$meanTS
# aa <- bio_ft_summary$meanTS
z <- log10(1+bio_ft_summary$Mean_depth/10)
# z <- bio_ft_summary$Mean_depth

# Modelo optimo:
modelo <- lm(aa~z)
summary(modelo)
a<-setDT(summary(modelo)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo)[8]),3), 
                 "p-value: ",
                 round(as.numeric(a),3))


ggplot(bio_ft_summary,aes(x=Mean_depth,y=meanTS))+
 
# ggplot(aes(x=Target_range,y=TS_comp, color=as.factor(id)))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm", color="black") +
  theme_minimal()+
  ylab("TS (dB)")+
  xlab("Depth (m)")+
  # ggtitle("RELACION del TS ~ Z (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  # annotate("text", x=30, y=0, label=my_text) +
  theme_classic()

ggsave("Figuras/ts_st_depth.tiff", width=18, height=18, units="cm")
# ggsave("ts_tilt_all.tiff", width=18, height=18, units="cm")




#  FIGURA TS ~ DEPTH + LENGTH----
modelo <- lm (data = bio_ft_summary, meanTS ~ log10(FL) + log10(1+Mean_depth/10))
summary(modelo)

#  FIGURA TS ~ DEPTH + LENGTH + TILT----
modelo <- lm (data = bio_ft_summary, meanTS ~ log10(FL) + log10(1+Mean_depth/10) + tilt)
summary(modelo)

# FIGURA DEPTH ~ TILT----
all_ts2 %>% 
  filter(Target_range<40) %>% 
ggplot(aes(x=Tilt_angle,y=-range2, color=id )) +
  geom_point(size=1.3, alpha=0.2) +
  geom_smooth(method = "loess",span = 0.8) +
  # geom_smooth(method = "gam") +
  theme_classic() +
  xlab("Tilt angle (degrees)") +
  ylab("Depth (m)") 

# TS-LENGTH + TILT----
summary(all_ts$id)

ggplot(all_ts, aes(x=as.factor(id), y=tilt)) +
  geom_boxplot()

# HISTOGRAMA + DENSITY PLOT + NUM FT (GRIS)----

ggplot(all_ts,aes( x = TS_comp)) +
  geom_histogram(aes(y = ..density..),binwidth = 2,
                 alpha = 0.8, show.legend = FALSE) +
  geom_density() +
  facet_wrap(~id,ncol=5) +
  theme_classic() +
  ggtitle ("TS freq. distr. for all the fish tracks") + 
  xlab ("Target Strength (dB)") + 
  ylab ("Count") +
  geom_text(aes(-10, .31, label=as.character(Num_Fish_Tracks)))

# ggsave("figuras/ts_distrib_fishtracks.tiff", width=17, height=17, units="cm")
# ggsave("Figuras/ts_distrib_fishtracks_ken6.tiff", width=17, height=17, units="cm")




# ridges DENSITY PLOT COLORES----


ggplot(all_ts, aes(x = TS_comp, y = as.factor(length))) +
  geom_density_ridges(scale=1) +
  theme_ridges() + 
  # theme(legend.position = "none") +
  # geom_vline(aes(xintercept = -45),
  # color = "#FC4E07",
  # linetype = "dashed", size = 0.1)+ 
  # ggtitle ("TS freq. distr. for all the fish tracks") + 
  xlab ("Target Strength (dB)") + 
  ylab ("Mean length (cm)") +
  # facet_wrap(~tipo)+ 
  theme(legend.position=c(.75,.25))+  
  # annotate("segment", x=-50, xend=-33, y=3, yend=3)+ 
  # annotate("segment", x=-31, xend=-5, y=3, yend=3)+
  annotate("rect", xmin=-33, xmax=-31, ymin=3, ymax=3.9, alpha=.5,)+ 
  annotate("segment", x=-33, xend=-33, y=3, yend=3.9)+ 
  annotate("segment", x=-31, xend=-31, y=3, yend=3.9)+ 
  annotate("segment", x=-33, xend=-31, y=3.9, yend=3.9)


# ggsave("Figuras/ggridges_fish_track_distrib_ALL_facet.tiff", width=17, height=17, units="cm")






























# ///////////////////
# OLD o PRUEBA///////----



# MODELOS LINEALES DE REGRESION-----


model_data<-all_ts %>% 
  group_by(id) %>% 
  filter(TS_mean>=-45&tipo=="BIO"&N>1)%>% 
  summarise(ts_media=meandB(TS_mean),ts_mediana=median(TS_mean),ts_moda=getmode(TS_mean),st_dev=sd(TS_mean),
            length=as.numeric(as.character(unique(length))),
            tipo=unique(tipo),
            N=unique(N),
            order=unique(order)) %>% 
  ungroup() %>% 
  mutate(b20_media=ts_media -(20*log10(as.numeric(length))),
         b20_mediana=ts_mediana -(20*log10(as.numeric(length))),
         b20_moda=ts_moda -(20*log10(as.numeric(length))))%>% 
  mutate(fuente=1)
str(model_data) 




# Modelo  TS (MEDIANA)------
aa <- as.numeric(as.character(model_data$ts_mediana))
bb <- as.numeric(as.character(model_data$length))
table<-as.data.frame(cbind(aa,bb))
table$id<-as.numeric(as.character(model_data$id)) 
table$sd<-as.numeric(as.character(model_data$st_dev)) 
# Modelo optimo:
modelo_TS <- lm(aa~log10(bb))
# model1<-lm(y ~ x) 
# model2 <- lm(y ~ 1 + offset(1.5*x))
summary(modelo_TS)
# Modelo ajustado a la pendiente 20
?offset
SE <- sd(table$aa)/sqrt(length(table$aa))
L_mean_slope  <-  20
flN  <-  transform(table,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(aa) ~ 1 + offset(log10(bb)*L_mean_slope), data = flN)
summary(m1)
# slope chosen by lm(...)
m0 <- lm(as.numeric(aa) ~ log10(bb), data=flN)
summary(m0)
a<-summary(m0)[4]

mediana_bio<-as.data.frame(
  c(coef(m1)[1],coef(m0)[1],coef(m0)[2],summary(m0)[8],a$coefficients[2,4],summary(m0)[6]))

names(mediana_bio)<-c("Forced b20 Intercept","Optimun Intercept","Optimun Slope","R Squared","P value","Res Stdrd Dev")


ggplot(flN,aes(x=log10(bb),y=aa)) +
  geom_errorbar(aes(ymin=aa-sd, ymax=aa+sd), width=.0051,col="grey",alpha=0.99) + 
  geom_point()+
  ggtitle("Target Strength (median) VS log Length ") + 
  xlab("Length (m) at log10 scale") + 
  ylab("Target Strength (dB)") + 
  theme_minimal()+ 
  geom_abline(slope=20, intercept=coef(m1)[1], color="tomato",size=1,lty="dashed")+ 
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1],color="darkblue",size=1)+
  geom_label(label=paste("Y = ",round(coef(m0)[1],2)," + ",round(coef(m0)[2],2),"X",sep=""),x=1.85,y=-20,label.padding = unit(0.7125, "lines"),label.size = 1.,color = "darkblue",fill="grey96")+
  geom_label(label=paste("Y = ",round(coef(m1)[1],2)," + 20X",sep=""),x=1.84,y=-22,label.padding = unit(0.7125, "lines"),label.size = 1,color = "tomato",fill="grey96")
ggsave("FIGURAS/LM_TS_MEDIANA_FT.tiff", width=17, height=17, units="cm")



aa <- as.numeric(as.character(model_data$b20_mediana))
bb <- as.numeric(as.character(model_data$length))
table<-as.data.frame(cbind(aa,bb))
table$id<-as.numeric(as.character(model_data$id)) 
table$sd<-as.numeric(as.character(model_data$st_dev)) 
# Modelo optimo:
modelo_TS <- lm(aa~log10(bb))
# model1<-lm(y ~ x) 
# model2 <- lm(y ~ 1 + offset(1.5*x))
summary(modelo_TS)
# Modelo ajustado a la pendiente 20
?offset
SE <- sd(table$aa)/sqrt(length(table$aa))
L_mean_slope  <-  20
flN  <-  transform(table,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(aa) ~ 1 + offset(log10(bb)*L_mean_slope), data = flN)
# slope chosen by lm(...)
m0 <- lm(as.numeric(aa) ~ log10(bb), data=flN)
a<-summary(m0)[4]
mediana_b20_bio<-as.data.frame(c(coef(m1)[1],coef(m0)[1],coef(m0)[2],
                                 summary(m0)[8],
                                 a$coefficients[2,4],
                                 summary(m0)[6]))
names(mediana_b20_bio)<-c("Forced b20 Intercept","Optimun Intercept","Optimun Slope","R Squared","P value","Res Stdrd Dev")
mediana_bio_all <-cbind(mediana_bio,mediana_b20_bio)

ggsave("FIGURAS/LM_TS_B20_MEDIANA_FT.tiff", width=17, height=17, units="cm")
ggplot(flN,aes(x=log10(bb),y=aa)) +
  geom_point()+
  geom_errorbar(aes(ymin=aa-sd, ymax=aa+sd), col="grey",alpha=0.99) + 
  ggtitle("Target Strength (median b20) VS log Length ") + 
  xlab("Length (m) at log10 scale") + 
  ylab("Target Strength (dB)") + 
  theme_minimal()+ 
  geom_abline(slope=20, intercept=coef(m1)[1], color="tomato",size=1,lty="dashed")+
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1],color="darkblue",size=1)+
  geom_label(label=paste("Y = ",round(coef(m0)[1],2)," + ",round(coef(m0)[2],2),"X",sep=""),x=1.85,y=-63,label.padding = unit(0.7125, "lines"),label.size = 1.,color = "darkblue",fill="grey96")+
  geom_label(label=paste("Y = ",round(coef(m1)[1],2)," + 20X",sep=""),x=1.84,y=-65,label.padding = unit(0.7125, "lines"),label.size = 1,color = "tomato",fill="grey96")
dev.off()
# -----

# Modelo  TS (MEDIA)------
aa <- as.numeric(as.character(model_data$ts_media))
bb <- as.numeric(as.character(model_data$length))
table<-as.data.frame(cbind(aa,bb))
table$id<-as.numeric(as.character(model_data$id)) 
table$sd<-as.numeric(as.character(model_data$st_dev)) 
# Modelo optimo:
modelo_TS <- lm(aa~log10(bb))
# model1<-lm(y ~ x) 
# model2 <- lm(y ~ 1 + offset(1.5*x))
summary(modelo_TS)
# Modelo ajustado a la pendiente 20
?offset
SE <- sd(table$aa)/sqrt(length(table$aa))
L_mean_slope  <-  20
flN  <-  transform(table,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(aa) ~ 1 + offset(log10(bb)*L_mean_slope), data = flN)
# slope chosen by lm(...)
m0 <- lm(as.numeric(aa) ~ log10(bb), data=flN)
a<-summary(m0)[4]
media_bio<-as.data.frame(c(coef(m1)[1],coef(m0)[1],coef(m0)[2],
                           summary(m0)[8],
                           a$coefficients[2,4],
                           summary(m0)[6]))
names(media_bio)<-c("Forced b20 Intercept","Optimun Intercept","Optimun Slope","R Squared","P value","Res Stdrd Dev")

ggsave("FIGURAS/LM_TS_MEDIA_FT.tiff", width=17, height=17, units="cm")
ggplot(flN,aes(x=log10(bb),y=aa)) +
  geom_errorbar(aes(ymin=aa-sd, ymax=aa+sd), width=.0051,col="grey",alpha=0.99) + 
  geom_point()+
  ggtitle("Target Strength (mean) VS log Length ") + 
  xlab("Length (m) at log10 scale") + 
  ylab("Target Strength (dB)") + ylim(-40, -15)+
  theme_minimal()+ 
  geom_abline(slope=20, intercept=coef(m1)[1], color="tomato",size=1,lty="dashed")+ 
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1],color="darkblue",size=1)+
  geom_label(label=paste("Y = ",round(coef(m0)[1],2)," + ",round(coef(m0)[2],2),"X",sep=""),x=1.84,y=-20,label.padding = unit(0.7125, "lines"),label.size = 1.,color = "darkblue",fill="grey96")+
  geom_label(label=paste("Y = ",round(coef(m1)[1],2)," + 20X",sep=""),x=1.83,y=-22,label.padding = unit(0.7125, "lines"),label.size = 1,color = "tomato",fill="grey96")
dev.off()



aa <- as.numeric(as.character(model_data$b20_media))
bb <- as.numeric(as.character(model_data$length))
table<-as.data.frame(cbind(aa,bb))
table$id<-as.numeric(as.character(model_data$id)) 
table$sd<-as.numeric(as.character(model_data$st_dev)) 
# Modelo optimo:
modelo_TS <- lm(aa~log10(bb))
# model1<-lm(y ~ x) 
# model2 <- lm(y ~ 1 + offset(1.5*x))
summary(modelo_TS)
# Modelo ajustado a la pendiente 20
?offset
SE <- sd(table$aa)/sqrt(length(table$aa))
L_mean_slope  <-  20
flN  <-  transform(table,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(aa) ~ 1 + offset(log10(bb)*L_mean_slope), data = flN)
# slope chosen by lm(...)
m0 <- lm(as.numeric(aa) ~ log10(bb), data=flN)
a<-summary(m0)[4]
media_b20_bio<-as.data.frame(c(coef(m1)[1],coef(m0)[1],coef(m0)[2],
                               summary(m0)[8],
                               a$coefficients[2,4],
                               summary(m0)[6]))
names(media_b20_bio)<-c("Forced b20 Intercept","Optimun Intercept","Optimun Slope","R Squared","P value","Res Stdrd Dev")
media_bio_all <-cbind(media_bio,media_b20_bio)



ggsave("FIGURAS/LM_TS_B20_MEDIA_FT.tiff", width=17, height=17, units="cm")
ggplot(flN,aes(x=log10(bb),y=aa)) +
  geom_errorbar(aes(ymin=aa-sd, ymax=aa+sd),col="grey",alpha=0.99) + 
  geom_point()+
  ggtitle("Target Strength (mean b20) VS log Length ") + 
  xlab("Length (m) at log10 scale") + 
  ylab("Target Strength (dB)") + 
  theme_minimal()+ 
  geom_abline(slope=20, intercept=coef(m1)[1], color="tomato",size=1,lty="dashed")+ 
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1],color="darkblue",size=1)+
  geom_label(label=paste("Y = ",round(coef(m0)[1],2)," + ",round(coef(m0)[2],2),"X",sep=""),x=1.85,y=-62,label.padding = unit(0.7125, "lines"),label.size = 1.,color = "darkblue",fill="grey96")+
  geom_label(label=paste("Y = ",round(coef(m1)[1],2)," + 20X",sep=""),x=1.84,y=-65,label.padding = unit(0.7125, "lines"),label.size = 1,color = "tomato",fill="grey96")
dev.off()

# Modelo  TS (MODA)------
aa <- as.numeric(as.character(model_data$ts_moda))
bb <- as.numeric(as.character(model_data$length))
table<-as.data.frame(cbind(aa,bb))
table$id<-as.numeric(as.character(model_data$id)) 
table$sd<-as.numeric(as.character(model_data$st_dev)) 
# Modelo optimo:
modelo_TS <- lm(aa~log10(bb))
# model1<-lm(y ~ x) 
# model2 <- lm(y ~ 1 + offset(1.5*x))
summary(modelo_TS)
# Modelo ajustado a la pendiente 20
?offset
SE <- sd(table$aa)/sqrt(length(table$aa))
L_mean_slope  <-  20
flN  <-  transform(table,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(aa) ~ 1 + offset(log10(bb)*L_mean_slope), data = flN)
# slope chosen by lm(...)
m0 <- lm(as.numeric(aa) ~ log10(bb), data=flN)
a<-summary(m0)[4]
moda_bio<-as.data.frame(c(coef(m1)[1],coef(m0)[1],coef(m0)[2],
                          summary(m0)[8],
                          a$coefficients[2,4],
                          summary(m0)[6]))
names(moda_bio)<-c("Forced b20 Intercept","Optimun Intercept","Optimun Slope","R Squared","P value","Res Stdrd Dev")


ggsave("FIGURAS/LM_TS_MODA_FT.tiff", width=17, height=17, units="cm")
ggplot(flN,aes(x=log10(bb),y=aa)) +
  geom_errorbar(aes(ymin=aa-sd, ymax=aa+sd), width=.0051,col="grey",alpha=0.99) + 
  geom_point()+
  ggtitle("Target Strength (mode) VS log Length ") + 
  xlab("Length (m) at log10 scale") + 
  ylab("Target Strength (dB)") + 
  theme_minimal()+ 
  geom_abline(slope=20, intercept=coef(m1)[1], color="tomato",size=1,lty="dashed")+ 
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1],color="darkblue",size=1)+
  geom_label(label=paste("Y = ",round(coef(m0)[1],2)," + ",round(coef(m0)[2],2),"X",sep=""),x=1.85,y=-20,label.padding = unit(0.7125, "lines"),label.size = 1.,color = "darkblue",fill="grey96")+
  geom_label(label=paste("Y = ",round(coef(m1)[1],2)," + 20X",sep=""),x=1.84,y=-24,label.padding = unit(0.7125, "lines"),label.size = 1,color = "tomato",fill="grey96")
dev.off()



aa <- as.numeric(as.character(model_data$b20_moda))
bb <- as.numeric(as.character(model_data$length))
table<-as.data.frame(cbind(aa,bb))
table$id<-as.numeric(as.character(model_data$id)) 
table$sd<-as.numeric(as.character(model_data$st_dev)) 
# Modelo optimo:
modelo_TS <- lm(aa~log10(bb))
# model1<-lm(y ~ x) 
# model2 <- lm(y ~ 1 + offset(1.5*x))
summary(modelo_TS)
# Modelo ajustado a la pendiente 20
?offset
SE <- sd(table$aa)/sqrt(length(table$aa))
L_mean_slope  <-  20
flN  <-  transform(table,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(aa) ~ 1 + offset(log10(bb)*L_mean_slope), data = flN)
# slope chosen by lm(...)
m0 <- lm(as.numeric(aa) ~ log10(bb), data=flN)

a<-summary(m0)[4]
moda_b20_bio<-as.data.frame(c(coef(m1)[1],coef(m0)[1],coef(m0)[2],
                              summary(m0)[8],
                              a$coefficients[2,4],
                              summary(m0)[6]))
names(moda_b20_bio)<-c("Forced b20 Intercept","Optimun Intercept","Optimun Slope","R Squared","P value","Res Stdrd Dev")
moda_bio_all <-cbind(moda_bio,moda_b20_bio)

ggsave("FIGURAS/LM_TS_B20_MODA_FT.tiff", width=17, height=17, units="cm")
ggplot(flN,aes(x=log10(bb),y=aa)) +
  geom_errorbar(aes(ymin=aa-sd, ymax=aa+sd), width=.0051,col="grey",alpha=0.99) + 
  geom_point()+
  ggtitle("Target Strength (mode b20) VS log Length ") + 
  xlab("Length (m) at log10 scale") + 
  ylab("Target Strength (dB)") + 
  theme_minimal()+ 
  geom_abline(slope=20, intercept=coef(m1)[1], color="tomato",size=1,lty="dashed")+ 
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1],color="darkblue",size=1)+
  geom_label(label=paste("Y = ",round(coef(m0)[1],2)," + ",round(coef(m0)[2],2),"X",sep=""),x=1.85,y=-62,label.padding = unit(0.7125, "lines"),label.size = 1.,color = "darkblue",fill="grey96")+
  geom_label(label=paste("Y = ",round(coef(m1)[1],2)," + 20X",sep=""),x=1.84,y=-65,label.padding = unit(0.7125, "lines"),label.size = 1,color = "tomato",fill="grey96")
dev.off()

# Tabla resumen----

resumen <- rbind(moda_bio_all,mediana_bio_all,media_bio_all) ;rownames(resumen)<-c("moda","mediana","media")
fwrite(resumen,file="Resultados/resumen_bio.csv",row.names = T)




