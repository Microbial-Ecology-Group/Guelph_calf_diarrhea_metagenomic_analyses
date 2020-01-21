### final figures
library(tidyverse)
library(RColorBrewer)
library(ggthemes)
#source('scripts/NCBA2_anosim.r')

### Start of code for figures
setkey(amr_melted_raw_analytic,ID) 
setkey(amr_melted_analytic,ID) 

setkey(microbiome_melted_analytic,ID)

# Set keys for both metadata files
setkey(metadata,ID)
setkey(microbiome_metadata,ID)

microbiome_melted_analytic <- microbiome_melted_analytic[microbiome_metadata]
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
amr_melted_analytic <- amr_melted_analytic[metadata]

#### AMR figures
###
##
# Heatmap at group level - by sample
AMR_group_samples <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, ID)]#[order(-sum_class )]

group_annotations <- annotations[,-1]
group_annotations
group_annotations <- unique(group_annotations,by = c("mechanism"))

setkey(group_annotations, mechanism)
setkey(AMR_group_samples , Name)

AMR_group_samples  <- group_annotations[AMR_group_samples ]
setkey(AMR_group_samples , mechanism)
AMR_group_samples  <- metadata[AMR_group_samples]

ggplot(data = AMR_group_samples , aes(x = ID, y = class)) +
  geom_tile(aes(fill = log_Normalized_Count)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=2.5) 



#write.csv(AMR_group_samples, "16S_norm_group_melted_sample_counts_.csv")


## Heatmap at Mechanism level - by sample
AMR_mech_samples <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, ID)]#[order(-sum_class )]

mech_annotations <- annotations[,-1]
mech_annotations <- unique(mech_annotations,by = c("mechanism"))

setkey(mech_annotations, mechanism)
setkey(AMR_mech_samples , Name)

AMR_mech_samples  <- mech_annotations[AMR_mech_samples ]
setkey(AMR_mech_samples , ID)
AMR_mech_samples  <- metadata[AMR_mech_samples]
#write.csv(AMR_mech_samples, "16S_norm_melted_sample_counts.csv")


#
##
### Heatmap at Mechanism level
##
#
AMR_mech_sample_ID <- amr_melted_raw_analytic[Level_ID=="Mechanism"][Normalized_Count > 0, .(num_samples= .N, Normalized_Count= sum(Normalized_Count),log_Normalized_Count= log(sum(Normalized_Count))),by=.(Name, ID)]#[order(-sum_class )]
mech_annotations <- annotations[,-c(1,4)]
mech_annotations <- unique(mech_annotations,by = c("mechanism"))
setkey(mech_annotations, mechanism)
setkey(AMR_mech_sample_ID, Name)

AMR_mech_sample_ID <- mech_annotations[AMR_mech_sample_ID]
#write.csv(AMR_mech_sample_groups, "16S_norm_melted_Group_counts.csv")
ggplot(data = AMR_mech_sample_ID , aes(x = ID, y = mechanism)) +
  geom_tile(aes(fill = log_Normalized_Count)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(label = num_samples), size=2.5) 


## Stacked bars - average sum of AMR class
AMR_class_mean <- amr_melted_analytic[Level_ID=="Class", .(sum_class= mean(Normalized_Count)),by=.(Group, Name)][order(-sum_class )]
ggplot(AMR_class_mean, aes(x = Group, y = sum_class, fill = Name)) + 
  geom_bar(stat = "identity")

#
##
###
#### Stacked 100% bars
###
##
#

## By ID, faceted by packaging
AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID, Name, Farm, Time)]#[order(-Packaging )]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(ID, Name) ]
AMR_class_sum$Class <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions

ggplot(AMR_class_sum, aes(x = ID, y = percentage, fill = Class)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Farm + Time, scales='free',ncol = 2) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=10),
    legend.title=element_text(size=14),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sampling Farm and Time') +
  ylab('Relative abundance\n') 
ggsave("~/Dropbox/Projects/MISC_projects/Weese/AMR_relabundance_Class_byFarmTime.jpeg", width = 30, height = 20, units = "cm")



## By Treatment, faceted by packaging
AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(Treatment, Name, Packaging)][order(-Packaging )]
AMR_class_sum$Name = droplevels(AMR_class_sum$Name)
AMR_class_sum[,total:= sum(sum_class), by=.(Treatment)]
AMR_class_sum[,percentage:= sum_class/total ,by=.(Treatment, Name) ]
AMR_class_sum$Class <- AMR_class_sum$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(Treatment, Name) ] removes some with low proportions
ggplot(AMR_class_sum, aes(x = Treatment, y = percentage, fill = Class, order = -as.numeric(Class))) + 
  geom_bar(stat = "identity")+
  #facet_grid( ~ Packaging, scales='free') +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample ID') +
  ylab('Relative abundance\n') 
ggsave("/home/enriquedoster/Dropbox/Projects/FC_retail_meat/final_figures/Figure2_AMR_barplot_byTreatment.jpeg", width = 30, height = 20, units = "cm")






### Figure 1 - AMR Class stacked graph ####

AMR_class_sum <- amr_melted_raw_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(Group, Name)][order(-sum_class )]
AMR_class_sum$Name = droplevels(AMR_class_sum$Name)
AMR_class_sum <- AMR_class_sum[with(AMR_class_sum, order(-sum_class)), ]
#AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=sort(levels(AMR_class_sum$Name ), FALSE))
#AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=c("Fluoroquinolones","Phenicol","Bacitracin","Cationic antimicrobial peptides", "Multi-drug resistance" ,"Aminoglycosides",   "betalactams" ,"MLS" ,"Tetracyclines"))
AMR_class_sum[,total:= sum(sum_class), by=.(Group)]
AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(Group, Name) ]
AMR_class_sum[percentage < .025,percentage:=as.character('')]
#AMR_class_sum[percentage > 2.5,percentage:= paste0(percentage,"%")]

#myColors <- rev(brewer.pal(14,"Set1"))
#names(myColors) <- levels(AMR_class_sum$Name)
#colScale <- scale_colour_manual(name = "Name",values = myColors)
AMR_class_sum$Group <- factor(AMR_class_sum$Group,levels = c("Arrival", "Exit", "Shipment"))


fig1 <- AMR_class_sum %>% 
  ggplot(aes(x=Group,y=sum_class, fill=Name, label= percentage)) +
  geom_bar(stat='identity')+
  #facet_wrap(~Cow_ID) +
  #scale_fill_manual(name = "Name",values = myColors) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample group') +
  ylab('16S normalized AMR gene abundance\n') 
#ggtitle('            Total AMR Gene Abundance and proportion at the class level, by sample group\n')
fig1
ggsave("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/figures/Figure1_AMR_barplot.jpeg", width = 30, height = 20, units = "cm")
#dev.off()
### By average of AMR class counts

AMR_class_mean <- amr_melted_raw_analytic[Level_ID=="Class", .(mean_class= mean(Normalized_Count)),by=.(Group, Name)]
AMR_class_mean$Name = droplevels(AMR_class_mean$Name)
AMR_class_mean <- AMR_class_mean[with(AMR_class_mean, order(-mean_class)), ]
AMR_class_mean$Name = factor(AMR_class_mean$Name ,levels=c("Phenicol","Bacitracin","Cationic antimicrobial peptides", "Multi-drug resistance" ,"Aminoglycosides",   "betalactams" ,"MLS" ,"Tetracyclines"))
AMR_class_mean[,total:= sum(mean_class), by=.(Group)]
AMR_class_mean[,percentage:= round(mean_class/total, digits=2) ,by=.(Group, Name) ]
AMR_class_mean[percentage < .03,percentage:='']

fig_mean <- AMR_class_mean %>% 
  ggplot(aes(x=Group,y=mean_class,fill=Name,label= percentage)) +
  geom_bar(stat='identity')+
  #facet_wrap(~Cow_ID) +
  scale_fill_manual(name = "Name",values = myColors) +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank()
  ) +
  xlab('\n Sample group') +
  ylab('16S normalized AMR gene abundance\n') +
  ggtitle('Average AMR Gene Abundance and proportion at the class level, for each treatment group at Day 1 and day 11 \n')
fig_mean


####
####

# MICROBIOME

####
####


microbiome_phylum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.( Name)]
microbiome_phylum[,total:= sum(sum_phylum)]
microbiome_phylum[,study_proportion:= sum_phylum/total ,by=.(Name) ]
microbiome_phylum$Name = droplevels(microbiome_phylum$Name)
rare_phyla <- microbiome_phylum[study_proportion < .0001 , Name]
test_sum <- microbiome_phylum[,.(total = sum(study_proportion))]
length(rare_phyla)

## Changing name of 
microbiome_phylum_edited <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count),Normalized_Count),by=.( ID, Name,Packaging, Treatment)]
microbiome_phylum_edited[Name %in% rare_phyla ,Name := 'Low abundance phyla']
microbiome_phylum_edited$Name = droplevels(microbiome_phylum_edited$Name)
microbiome_phylum_edited[,total:= sum(Normalized_Count), by=.(ID)]
microbiome_phylum_edited[,percentage:= sum_phylum/total ,by=.(ID, Name,Packaging, Treatment) ]
microbiome_phylum_edited$Phylum <- microbiome_phylum_edited$Name

ggplot(microbiome_phylum_edited, aes(x = ID, y = percentage, fill = Phylum)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging,scales = "free") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample ID') +
  ylab('Relative abundance\n') 




myColors <- rev(brewer.pal(5,"Set3"))
names(myColors) <- levels(microbiome_phylum_mean$Name)
colScale <- scale_colour_manual(name = "Name",values = myColors)

microbiome_phylum_mean$Group <- factor(microbiome_phylum_mean$Group,levels = c("Day1-Treated", "Day1-Untreated", "Day11-Treated", "Day11-Untreated"))






rare_phyla <- Micro_phylum_sum[percentage < .05,Name]

## Changing name of 
Micro_phylum_sum[Name %in% rare_phyla ,Name := 'Low abundance phyla']

microbiome_phylum_edited[,total:= sum(Normalized_Count), by=.(Group)]





## By ID, faceted by packaging
Micro_phylum_sum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(ID, Name, Packaging, Treatment)][order(-Treatment )]
Micro_phylum_sum$Name = droplevels(Micro_phylum_sum$Name)
Micro_phylum_sum[,total:= sum(sum_phylum), by=.(ID)]
Micro_phylum_sum[,percentage:= sum_phylum/total ,by=.(ID, Name) ]
Micro_phylum_sum$Phylum <- Micro_phylum_sum$Name
#AMR_class_sum[,percentage:= round(sum_phylum/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(Micro_phylum_sum, aes(x = ID, y = percentage, fill = Phylum)) + 
  geom_bar(stat = "identity")+
  facet_wrap( ~ Packaging,scales = "free") +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample ID') +
  ylab('Relative abundance\n') 



### Figure 5 - Microbiome Phylum stacked graph ####

fig1 <- microbiome_melted_analytic %>% 
  filter(Level_ID=="Phylum") %>%
  ggplot(aes(x=reorder(ID))) +
  geom_bar(aes(fill=factor(Name), weight=Normalized_Count))+
  #facet_wrap(~Cow_ID) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample type') +
  ylab('Count of raw observations\n') +
  ggtitle('raw AMR Mechanism hits - by cow ID \n')
fig1


microbiome_phylum <- microbiome_melted_analytic[Level_ID=="Phylum", .(sum_phylum= sum(Normalized_Count)),by=.(Group, Name)]
microbiome_phylum[,total:= sum(sum_phylum), by=.(Group)]
microbiome_phylum[,proportion:= sum_phylum/total ,by=.(Group, Name) ]
rare_phyla <- microbiome_phylum[proportion < .03,Name]

## Changing name of 
microbiome_phylum_edited <- microbiome_melted_analytic[Level_ID=="Phylum", .(Normalized_Count,Group),by=.(Name,ID)]
microbiome_phylum_edited[Name %in% rare_phyla ,Name := 'Low abundance phyla']
microbiome_phylum_edited[,total:= sum(Normalized_Count), by=.(Group)]

microbiome_phylum_mean <- microbiome_phylum_edited[, .(mean_phylum = mean(Normalized_Count)) ,by=.(Group,Name) ]
microbiome_phylum_mean[, total := sum(mean_phylum) ,by=.(Group) ]
microbiome_phylum_mean[, percentage := round(mean_phylum/total,digits=2) ,by=.(Group, Name) ]
microbiome_phylum_mean[, percentage:= as.character(percentage)]
microbiome_phylum_mean[percentage=='0', percentage := '<0.01']

microbiome_phylum_mean$Name = droplevels(microbiome_phylum_mean$Name)
#microbiome_phylum_mean$Name = rev(reorder(microbiome_phylum_mean$Name, X=as.numeric(microbiome_phylum_mean$phylum), FUN=sum))
microbiome_phylum_mean$Name = factor(microbiome_phylum_mean$Name ,levels=c("Low abundance phyla","Actinobacteria","Proteobacteria","Bacteroidetes","Firmicutes"))

myColors <- rev(brewer.pal(5,"Set3"))
names(myColors) <- levels(microbiome_phylum_mean$Name)
colScale <- scale_colour_manual(name = "Name",values = myColors)

microbiome_phylum_mean$Group <- factor(microbiome_phylum_mean$Group,levels = c("Day1-Treated", "Day1-Untreated", "Day11-Treated", "Day11-Untreated"))


fig1 <- microbiome_phylum_mean %>% 
  ggplot(aes(x=Group,y=mean_phylum, fill=Name, label= percentage)) +
  geom_bar(stat='identity')+
  #facet_wrap(~Cow_ID) +
  scale_fill_manual(name = "Name",values = myColors) +
  geom_text(size = 5, position = position_stack(vjust = 0.5)) +
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=24),
    strip.text.y=element_text(size=24, angle=0),
    axis.text.x=element_text(size=16, angle=20, hjust=1),
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  xlab('\n Sample group') +
  ylab('CSS normalized counts\n') 
#ggtitle('Average phyla counts and proportion, by sample group \n')
fig1

ggsave("/home/enrique/Dropbox/WRITING/NCBA2_JAN2018/Microbiome_manuscript/FrontiersMicroSubmission/Review/Review_submission_files/figures/Figure5_Microbiome_barplot_update.jpeg", width = 40, height = 30, units = "cm")



