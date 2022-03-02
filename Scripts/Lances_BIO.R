
##======================================================##
## SCRIPT : BFT TS PUBLICATION 2020
##
## Authors : Jon Uranga &  Guillermo Boyra
## Last update : 20-10-2020
## Description: 
##              Este script lee todas las exportaciones de TS d elos lances 
##              registrados en la campaña de BFT Index y calcula el TS tipico del bluefin tuna 
##              
##======================================================##


## Load packages  -------------------------------------------------------------------------
library(data.table)
library(stringr)
library(BuoysAZTI) # Own package with neccesary functions and tables to run the script
library(dplyr)
library(RPostgreSQL)
library(tidyr)
library(seewave)
library(ggplot2)
source("getmode.R")
library(ggridges)
library(geosphere)

# PREPARACION DE datos

# //////////////    BIO    ////////////----


# Muestreo biologico
muestreo<-fread("muestreo_2016_bft_index.csv")
muestreo<-muestreo[-6,]


# Cargamos los SINGLE TARGET DE LOS FISH TRACKS NEW TOBY-----

idx.files <- list.files("exports/",pattern="(targets)");idx.files
idx.files<-as.data.frame(idx.files)
names(idx.files)<-"source"
idx.files$calera<-as.numeric(sub("\\D+","",idx.files$source))
idx.files <- idx.files %>% arrange(calera) 

all_ts<-c()
for(i in 1:dim(idx.files)[1]){
  
  ft <- setDT(read.table(file=paste("exports/",idx.files$source[i],sep=""), sep=",",header = TRUE))
  
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
            length = round(mean(length),1), 
            sdev_TS = round(sd(TS_comp ),1),
            sdev_depth = round(sd(Target_range   ),1),
            Num_Fish_Tracks = unique(Num_Fish_Tracks),
            N_Sampled_Fish = unique(N_Sampled_Fish),
            Mean_depth = round(mean(Target_range),0)) %>% 
  ungroup() 

# Modelo optimo:
modelo <- lm(bio_ft_summary$meanTS ~ (bio_ft_summary$length))
summary(modelo)
fwrite(bio_ft_summary,"Resumen_lances.csv")

bio_ft_summary_order <- bio_ft_summary %>% arrange(length) ; bio_ft_summary_order
a<-setDT(summary(modelo)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo)[8]),2), 
                 "p-value: ",
                 round(as.numeric(a),7))
my_text2 <- paste("Slope (a):" ,
                  round(as.numeric(setDT(summary(modelo)[4])[2,1]),3), 
                  "Intercept (b): ",
                  round(as.numeric(setDT(summary(modelo)[4])[1,1]),3))

ggplot(bio_ft_summary,aes(x=length,y=meanTS))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin=meanTS-(sdev_TS/2), ymax=meanTS+(sdev_TS/2)), width=.0051,col="grey30",alpha=0.49) + 
  theme_minimal()+
  xlab("Log (Length)")+
  ylab("Target Strength")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=80.5, y=-13, label=my_text,col="maroon") +
  annotate("text", x=77, y=-15, label=my_text2,col="maroon") +
  scale_x_log10() +
  theme_classic()+geom_text(label=bio_ft_summary$id,mapping = aes(x=length+0.75,y=meanTS+0.75))

ggsave("figuras/ts_log_length_ken6.tiff", width=23, height=18, units="cm")



# HISTOGRAMAS FACET

ggplot(all_ts, aes(x=TS_comp)) +
  geom_histogram(aes(y = stat(density)), fill="darkslategray4",binwidth = 2) +
  geom_density(alpha = 0.2, fill = "grey50",color="grey30") +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = -45),color = "#FC4E07",linetype = "dashed", size = 0.1)+
  ggtitle ("TS freq. distr. for all sets") + 
  xlab ("Target Strength (dB)") + 
  ylab ("Frequency") +
  facet_wrap(~id,nrow=2,scales = "free")+theme_classic()

ggsave("FIGURAS/ggridges_fish_track_distrib_ALL_facet.tiff", width=27, height=17, units="cm")





# FIGURA RIDGES DEL TS MEDIO DE LOS FT EN MUESTREOS BIOLOGICOS
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
  geom_density_ridges(panel_scaling = T) +
  theme_ridges() + 
  theme(legend.position = "none") +
  # geom_vline(aes(xintercept = -45),color = "#FC4E07",linetype = "dashed", size = 0.1)+ 
  ggtitle ("TS freq. distr. and fish lengths for all sets ") + 
  xlab ("Target Strength (dB)") + 
  ylab ("Mean fish length (cm)")  +
  geom_text(aes(x = -5, y=0.64,label="FT(N)"),color="maroon",)+
  geom_text(aes(x = -1, y=0.64,label="Set"),color="blue")+
  geom_text(aes(x = -5, y=1.3,label=bio_ft_summary_order$Num_Fish_Tracks[1]),color="maroon")+
  geom_text(aes(x = -5, y=2.3,label=bio_ft_summary_order$Num_Fish_Tracks[2]),color="maroon")+
  geom_text(aes(x = -5, y=3.3,label=bio_ft_summary_order$Num_Fish_Tracks[3]),color="maroon")+
  geom_text(aes(x = -5, y=4.3,label=bio_ft_summary_order$Num_Fish_Tracks[4]),color="maroon")+
  geom_text(aes(x = -5, y=5.3,label=bio_ft_summary_order$Num_Fish_Tracks[5]),color="maroon")+
  geom_text(aes(x = -5, y=6.3,label=bio_ft_summary_order$Num_Fish_Tracks[6]),color="maroon")+
  geom_text(aes(x = -5, y=7.3,label=bio_ft_summary_order$Num_Fish_Tracks[7]),color="maroon")+
  geom_text(aes(x = -5, y=8.3,label=bio_ft_summary_order$Num_Fish_Tracks[8]),color="maroon")+
  geom_text(aes(x = -5, y=9.3,label=bio_ft_summary_order$Num_Fish_Tracks[9]),color="maroon")+
  geom_text(aes(x = -5, y=10.3,label=bio_ft_summary_order$Num_Fish_Tracks[10]),color="maroon")+
  # geom_text(aes(x = -5, y=11.3,label=bio_ft_summary_order$Num_Fish_Tracks[11]),color="maroon")+
  geom_text(aes(x = -1, y=1.3,label=bio_ft_summary_order$id[1]),color="blue")+
  geom_text(aes(x = -1, y=2.3,label=bio_ft_summary_order$id[2]),color="blue")+
  geom_text(aes(x = -1, y=3.3,label=bio_ft_summary_order$id[3]),color="blue")+
  geom_text(aes(x = -1, y=4.3,label=bio_ft_summary_order$id[4]),color="blue")+
  geom_text(aes(x = -1, y=5.3,label=bio_ft_summary_order$id[5]),color="blue")+
  geom_text(aes(x = -1, y=6.3,label=bio_ft_summary_order$id[6]),color="blue")+
  geom_text(aes(x = -1, y=7.3,label=bio_ft_summary_order$id[7]),color="blue")+
  geom_text(aes(x = -1, y=8.3,label=bio_ft_summary_order$id[8]),color="blue")+
  geom_text(aes(x = -1, y=9.3,label=bio_ft_summary_order$id[9]),color="blue")+
  geom_text(aes(x = -1, y=10.3,label=bio_ft_summary_order$id[10]),color="blue")

# ggsave("figuras/ggridges_por_BIO.tiff", width=17, height=17, units="cm")
ggsave("figuras/ggridges_por_BIO_ken6.tiff", width=17, height=17, units="cm")




# CALCULAMOS EL ANGULO ----

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



# FIGURA  TS ~ TILT    TOBY

ggplot(all_ts2,aes(x=Tilt_angle,y=TS_comp ))+
  geom_point(size=1.3, alpha=0.2)+
  geom_smooth(method = "loess",span = 0.8) +
  # geom_smooth(method = "gam") +
  theme_minimal()+
  xlab("Tilt")+
  ylab("Target Strength")+
  ggtitle("TS~Tilt for all sets")



# ggsave("ts_tilt_FT_loess.tiff", width=18, height=18, units="cm")
ggsave("figuras/ts_tilt_FT_gam.tiff", width=21, height=17, units="cm")
# ggsave("ts_tilt_FT_gam_sin_085.tiff", width=18, height=18, units="cm")

# ----



# FIGURA  LENGTH ~ DETPH

# Modelo optimo:
modelo <- lm(all_ts$length ~ (all_ts$Target_range))
a<-setDT(summary(modelo)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo)[8]),3),
                 "p-value: ",
                 round(as.numeric(a),3))



ggplot(all_ts,aes(x=as.numeric(all_ts$length),y=Target_range))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  # geom_errorbar(aes(ymin=meanTS-sdev, ymax=meanTS+sdev), width=.0051,col="grey",alpha=0.49) + 
  theme_minimal()+
  xlab("Length (cm)")+
  ylab("Mean depth (m)")+
  # ggtitle("RELACION del TS ~ LENGTH (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  annotate("text", x=94, y=8, label=my_text) +
  theme_classic()

# ggsave("depth_length.tiff", width=18, height=18, units="cm")
ggsave("depth_length_ken6.tiff", width=18, height=18, units="cm")




# FIGURA  TS ~ DEPTH

# aa <- all_ts$TS_comp
aa <- bio_ft_summary$meanTS
# z <- all_ts$Target_range
z <- bio_ft_summary$Mean_depth

# Modelo optimo:
modelo <- lm(aa~z)
a<-setDT(summary(modelo)[4])[2,4]
my_text <- paste("Multiple R-squared:" ,
                 round(as.numeric(summary(modelo)[8]),3), 
                 "p-value: ",
                 round(as.numeric(a),3))


ggplot(bio_ft_summary,aes(x=meanTS,y=-Mean_depth))+
# ggplot(all_ts,aes(x=Target_range,y=TS_comp))+
  geom_point(size=1.3, alpha=0.5)+
  geom_smooth(method = "lm") +
  theme_minimal()+
  xlab("TS")+
  ylab("Z")+
  ggtitle("RELACION del TS ~ Z (FISH TRACKS)")+
  # theme(plot.title = element_text(hjust = 0.5))+
  # annotation_custom(my_grob)+
  # annotate("text", x=30, y=0, label=my_text) +
  theme_classic()

ggsave("ts_tilt_mean.tiff", width=18, height=18, units="cm")
# ggsave("ts_tilt_all.tiff", width=18, height=18, units="cm")









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
ggsave("figuras/ts_distrib_fishtracks_ken6.tiff", width=17, height=17, units="cm")




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


ggsave("FIGURAS/ggridges_fish_track_distrib_ALL_facet.tiff", width=17, height=17, units="cm")






























# ///////////////////
# OLD o PRUEBA///////



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
# -----

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
# -----

# Tabla resumen

resumen <- rbind(moda_bio_all,mediana_bio_all,media_bio_all) ;rownames(resumen)<-c("moda","mediana","media")
fwrite(resumen,file="Resultados/resumen_bio.csv",row.names = T)




