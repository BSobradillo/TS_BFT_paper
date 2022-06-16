
##======================================================##
## SCRIPT : BFT TS PUBLICATION 2020
##
## Authors : Jon Uranga,  Guillermo Boyra & Bea Sobradillo
## Last update : 02-03-2022
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


# Biological sampling----
muestreo<-fread("Datos/muestreo_2016_bft_index.csv")
muestreo<-muestreo[-6,]


# Fish tracks----
# Cargamos los SINGLE TARGET DE LOS FISH TRACKS NEW TOBY

idx.files <- list.files("Datos/exports/",pattern="(targets)");idx.files
idx.files<-as.data.frame(idx.files)
names(idx.files)<-"source"
idx.files$calera<-as.numeric(sub("\\D+","",idx.files$source))
idx.files <- idx.files %>% arrange(calera) 

all_ts<-c()
for(i in 1:dim(idx.files)[1]){
  
  ft <- setDT(read.table(file=paste("Datos/exports/",idx.files$source[i],sep=""), sep=",",header = TRUE))
  
  # idx <- which(fish_tracks$id==fish_tracks_filter[i])
  
  id <- i
  Num_Fish_Tracks <- length(unique(ft$Region_name))
  length <- muestreo$MEDIA[i]
  N_Sampled_Fish <- length(which(!is.na(muestreo[i,8:30])))
  ft<-cbind(ft, id, length, N_Sampled_Fish, Num_Fish_Tracks)
  all_ts <- rbind(all_ts,ft)
  print(i)
  
}
setDT(all_ts)
table(all_ts$id)

# all_ts <- all_ts[-which(all_ts$id==6)]
# all_ts <- all_ts[-which(all_ts$id==9)]



# ////   RESUMEN

library(seewave)

bio_ft_summary <- all_ts %>% 
  mutate(length = as.numeric(as.character(length))) %>%
  group_by(id) %>% 
  summarise(meanTS = round(meandB(TS_comp),1),
            N_TS = length(TS_comp),
            length = round(mean(length),1),
            N_length = length(length),
            sdev_TS = round(sd(TS_comp ),1),
            sdev_depth = round(sd(Target_range   ),1),
            Num_Fish_Tracks = unique(Num_Fish_Tracks),
            N_Sampled_Fish = unique(N_Sampled_Fish),
            Mean_depth = round(mean(Target_range),0), 
            N_depth = length(Target_range)) %>%
  ungroup() 

# FIGURA TS-FL: (BS: añado log10 de L)----
modelo <- lm(bio_ft_summary$meanTS ~ log10(bio_ft_summary$length))
summary(modelo)
modelo$coefficients
# fwrite(bio_ft_summary,"Datos/Resumen_lances.csv")

bio_ft_summary_order <- bio_ft_summary %>% arrange(length) ; bio_ft_summary_order

# Bea
  
  
meanTS <- meandB(bio_ft_summary$meanTS); meanTS
medianTS <- median(bio_ft_summary$meanTS); medianTS
mean_sigma_w <- mean(10^(bio_ft_summary$meanTS/10)*bio_ft_summary$Num_Fish_Tracks/sum(bio_ft_summary$Num_Fish_Tracks)); mean_sigma_w
mean_TS_w <- 10*log10(mean_sigma_w); mean_TS_w

# Despues de hablar con Guille vemos que es un TS demasiado bajo. Hay que hacer la ponderacion un poco diferente de manera que no disminuya tanto el valor: 
# media = suma(sigma_ponderado)/suma(peso de ponderacion)
# Lo mismo hay que aplicar a la talla

ponderate_data <- bio_ft_summary %>% 
  mutate(
    sigma = 10^(meanTS/10), 
    pon_weight = Num_Fish_Tracks/sum(Num_Fish_Tracks),
    sigma_pon = sigma*pon_weight,
    length_pon = length*pon_weight) %>% 
  summarise(
    mean_sigma = sum(sigma_pon)/sum(pon_weight),
    TS_pon = 10*log10(mean_sigma),
    sdTS = mean(sdev_TS),
    length_pon = sum(length_pon)/sum(pon_weight),
    # sd_TS = sddB(meanTS),
    
    sd_length = sd(length), 
    min_L = min(length),
    max_L = max(length),
    mean_z = mean(Mean_depth),
    sd_z = mean(sdev_depth))


a<-setDT(summary(modelo)[4])[2,4]; a
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo)[8]),2), 
                 "\np-value: ",
                 round(as.numeric(a),7))
my_text2 <- paste("- Slope (a):" ,
                  round(as.numeric(setDT(summary(modelo)[4])[2,1]),2), 
                  "\nIntercept (b): ",
                  round(as.numeric(setDT(summary(modelo)[4])[1,1]),2))

ggplot(bio_ft_summary,aes(x=length,y=meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=meanTS-(sdev_TS/2), ymax=meanTS+(sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  theme_minimal()+
  xlab("Log-Length (cm)")+
  ylab("TS (dB)")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=60, y=-15, label=my_text, hjust=0) +
  annotate("text", x=60, y=-19, label=my_text2, hjust=0) +
  scale_x_log10() +
  theme_classic() +
# theme(axis.line = element_line(arrow = arrow(type='closed', length = unit(10,'pt')))) +
  geom_text(label=bio_ft_summary$id,mapping = aes(x=length+0.75,y=meanTS+0.75))

# ggsave("Figuras/ts_log_length_ken6.tiff", width=17, height=12, units="cm")


# slope fixed to 20
?offset

SE <- sd(bio_ft_summary$meanTS)/sqrt(length(bio_ft_summary$meanTS))
L_mean_slope  <-  20
flN  <-  transform(bio_ft_summary,w=1/as.numeric(SE)^2)

# slope fixed
m1 <- lm(as.numeric(meanTS) ~ 1 + offset(log10(flN$length)*L_mean_slope), data = flN)
m1$coefficients
summary(m1)
# slope chosen by lm(...)
m0 <- lm(as.numeric(meanTS) ~ log10(length), data=flN)
m0$coefficients
summary(m0)


# b20 linear model
my_text_b20 <- paste("- - Slope (a):" ,
                  20, 
                  "\nIntercept (b): ",
                  round(as.numeric(setDT(summary(m1)[4])[1,1]),2))

ggplot(bio_ft_summary,aes(x=length,y=meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  # geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=meanTS-(sdev_TS/2), ymax=meanTS+(sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  geom_abline(slope=20, intercept=coef(m1)[1], lty=2) +
  geom_abline(slope=coef(m0)[2], intercept=coef(m0)[1]) +
  theme_minimal()+
  xlab("Length (cm)")+
  ylab("TS (dB)")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=60, y=-15, label=my_text, hjust=0) +
  annotate("text", x=60, y=-19, label=my_text2, hjust=0) +
  annotate("text", x=120, y=-35, label=my_text_b20, hjust=0) +

  scale_x_log10() +
  theme_classic() +
  # theme(axis.line = element_line(arrow = arrow(type='closed', length = unit(10,'pt')))) +
  geom_text(label=bio_ft_summary$id,mapping = aes(x=length+0.75,y=meanTS+0.75))

ggsave("Figuras/TS-logL-b20.tiff", width=17, height=12, units="cm")




# TS histogram facets----

ggplot(all_ts, aes(x=TS_comp)) +
  geom_histogram(aes(y = stat(density)), fill="white",binwidth = 2,color="black") +
  geom_density(alpha = 0.2, fill = "white",color="black", size=0.75) +
  theme(legend.position = "none") +
  # geom_vline(aes(xintercept = -45),color = "#FC4E07",linetype = "dashed", size = 0.1)+
  # ggtitle ("TS freq. distr. for all sets") + 
  xlab ("TS (dB)") + 
  ylab ("Frequency") +
  theme(text = element_text(size = 12)) +
  facet_wrap(~id,nrow=2,scales = "free_y") 
# +
  # theme_classic()

ggsave("Figuras/ggridges_fish_track_distrib_ALL_facet.tiff", width = 27, height=12,units="cm", dpi = 300)
       # width=27, height=17, units="cm")





# FIGURA RIDGES DEL TS MEDIO DE LOS FT EN MUESTREOS BIOLOGICOS----
# CON ETIQUETAS DE N Y NUM DE FT

library(ggplot2)
library(ggridges)
theme_set(theme_minimal())


ggplot(all_ts, aes(x = TS_comp, y = as.factor(length), fill = factor(stat(quantile)))) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    quantiles = c(0.025, 0.975)
  ) +
  scale_fill_manual(
    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
  )


ggplot(data=all_ts, aes(x = TS_comp, y = as.factor(round(length,2)))) +
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
  geom_text(aes(x = -5, y=1.3,label=bio_ft_summary_order$Num_Fish_Tracks[1]),color="grey30")+
  geom_text(aes(x = -5, y=2.3,label=bio_ft_summary_order$Num_Fish_Tracks[2]),color="grey30")+
  geom_text(aes(x = -5, y=3.3,label=bio_ft_summary_order$Num_Fish_Tracks[3]),color="grey30")+
  geom_text(aes(x = -5, y=4.3,label=bio_ft_summary_order$Num_Fish_Tracks[4]),color="grey30")+
  geom_text(aes(x = -5, y=5.3,label=bio_ft_summary_order$Num_Fish_Tracks[5]),color="grey30")+
  geom_text(aes(x = -5, y=6.3,label=bio_ft_summary_order$Num_Fish_Tracks[6]),color="grey30")+
  geom_text(aes(x = -5, y=7.3,label=bio_ft_summary_order$Num_Fish_Tracks[7]),color="grey30")+
  geom_text(aes(x = -5, y=8.3,label=bio_ft_summary_order$Num_Fish_Tracks[8]),color="grey30")+
  geom_text(aes(x = -5, y=9.3,label=bio_ft_summary_order$Num_Fish_Tracks[9]),color="grey30")+
  geom_text(aes(x = -5, y=10.3,label=bio_ft_summary_order$Num_Fish_Tracks[10]),color="grey30")+
  geom_text(aes(x = -5, y=11.3,label=bio_ft_summary_order$Num_Fish_Tracks[11]),color="grey30")
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

ggsave("figuras/ggridges_por_BIO.tiff", width=15, height=15, units="cm")
# ggsave("Figuras/ggridges_por_BIO_ken6.tiff", width=17, height=17, units="cm")




# Angle calculation (Toby) ----

all_ts2<-as.data.frame(all_ts)
all_ts2 <- all_ts2 %>% 
  group_by(Region_name,id) %>% 
  mutate(lat2 = lag(Target_latitude),
         lon2 = lag(Target_longitude),
         range2 = lag(Target_range)) %>% 
  ungroup()

all_ts2<- all_ts2[which(!is.na(all_ts2$lon2)),]



setDT(all_ts2)
all_ts2$x<-NA
all_ts2$y<-NA
for(i in 1 : dim(all_ts2)[1]){
  all_ts2$x[i] <-as.numeric(distm(c(all_ts2$Target_longitude[i], 
                                    all_ts2$Target_latitude[i]),
                                  c(all_ts2$lon2[i],
                                    all_ts2$lat2[i]), fun = distHaversine))
  all_ts2$y[i] <-as.numeric(all_ts2$range2[i]-all_ts2$Target_range[i])
  print(i)
}

all_ts2 <- all_ts2 %>% 
  mutate(Tilt_angle = atan(y/x)*180/pi,Tilt_angle_abs = abs(Tilt_angle))

proba <- all_ts2[which(!is.na(all_ts2$Tilt_angle)),]
all_ts3<-all_ts2[append(which(all_ts2$Tilt_angle > 0.85 ), which(all_ts2$Tilt_angle < (-0.85))),]

# kopia<-all_ts2
# all_ts2<-kopia
# all_ts2<- all_ts2[which(all_ts2$Tilt_angle_abs>2),]
summary(all_ts2)
all_ts2<-all_ts2[which(!is.na(all_ts2$Tilt_angle)),]

ggplot(all_ts2,aes(x=as.factor(length),y=TS_comp))+
  geom_point()+
  theme_light()  



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

model_TS_alpha <- loess(TS_comp~Tilt_angle, all_ts2)
summary(model_TS_alpha)

  ggplot(all_ts2, aes(x=Tilt_angle,y=TS_comp )) +
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
    
  ggsave("Figuras/TS-tilt-tuna.tiff", width = 17, height=12,units="cm", dpi = 300)


  

  # ggtitle("TS~Tilt for all sets")



# ggsave("ts_tilt_FT_loess.tiff", width=18, height=18, units="cm")
# ggsave("figuras/ts_tilt_FT_gam.tiff", width=21, height=17, units="cm")
# ggsave("ts_tilt_FT_gam_sin_085.tiff", width=18, height=18, units="cm")





# FIGURA  LENGTH ~ DEPTH----

# Modelo optimo:
modelo_L_Z <- lm(all_ts$length ~ (all_ts$Target_range))
  summary(modelo_L_Z)
a<-setDT(summary(modelo_L_Z)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo_L_Z)[8]),3),
                 "\np-value: ",
                 round(as.numeric(a),3))



ggplot(all_ts,aes(x=as.numeric(all_ts$length),y=-Target_range))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  # geom_errorbar(aes(ymin=meanTS-sdev, ymax=meanTS+sdev), width=.0051,col="grey",alpha=0.49) + 
  scale_y_continuous(breaks=c(0,-20,-40,-60), labels=c("0","20","40","60")) +
  theme_minimal()+
  xlab("Length (cm)")+
  ylab("Mean depth (m)")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=120, y=-55, label=my_text, hjust=0) +
  theme_classic()
  # theme(axis.line = element_line(arrow = arrow(type='closed', length = unit(10,'pt')))) 

# ggsave("depth_length.tiff", width=18, height=18, units="cm")
ggsave("Figuras/depth_length_ken6.tiff", width=17, height=12, units="cm")




# FIGURA  TS ~ DEPTH----

# aa <- all_ts$TS_comp
aa <- bio_ft_summary$meanTS
# z <- all_ts$Target_range
z <- bio_ft_summary$Mean_depth

# Modelo optimo:
modelo <- lm(aa~z)
summary(modelo)
a<-setDT(summary(modelo)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo)[8]),3), 
                 "p-value: ",
                 round(as.numeric(a),3))


ggplot(bio_ft_summary,aes(x=meanTS,y=-Mean_depth))+
# ggplot(all_ts,aes(x=Target_range,y=TS_comp))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm", color="black") +
  theme_minimal()+
  xlab("TS (dB)")+
  ylab("Depth (m)")+
  # ggtitle("RELACION del TS ~ Z (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  # annotate("text", x=30, y=0, label=my_text) +
  theme_classic()

ggsave("Figuras/ts_depth.tiff", width=18, height=18, units="cm")
# ggsave("ts_tilt_all.tiff", width=18, height=18, units="cm")

model_depth <- lm(bio_ft_summary$Mean_depth ~ bio_ft_summary$meanTS)
summary(model_depth)


#  FIGURA TS ~ DEPTH + LENGTH----
modelo <- lm (data = bio_ft_summary, meanTS ~ log10(length) + log10(Mean_depth))
summary(modelo)

# FIGURA DEPTH ~ TILT
all_ts2 %>% 
  filter(Target_range<40) %>% 
ggplot(aes(x=Tilt_angle,y=-range2 )) +
  geom_point(size=1.3, alpha=0.2) +
  geom_smooth(method = "loess",span = 0.8) +
  # geom_smooth(method = "gam") +
  theme_classic() +
  xlab("Tilt angle (degrees)") +
  ylab("Depth (m)") 

# HISTOGRAMA + DENSITY PLOT + NUM FT (GRIS)

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




# ridges DENSITY PLOT COLORES


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




