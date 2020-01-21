# Rebecca figures
library(tidyverse)
library(RColorBrewer)
library(ggthemes)

# Make sure to run the "staging_script.R" first. Then you can modify the main objects created automatically.
## Run the script that handles resistome data.
#source('staging_script.R')

# Set keys for both metadata files, add metadata
setkey(amr_melted_raw_analytic,ID) 
setkey(amr_melted_analytic,ID) 
setkey(metadata,ID)
amr_melted_raw_analytic <- amr_melted_raw_analytic[metadata]
amr_melted_analytic <- amr_melted_analytic[metadata]


################################################
## Figure 1 - Sample relative abundance plots ##
################################################

# Sum class abundance by sample
AMR_class_sum <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(ID, Name, Time, Group)][order(-Group )]
AMR_class_sum[,total:= sum(sum_class), by=.(ID)][,percentage:= sum_class/total ,by=.(ID, Name) ]


#Determine order of abundant classes
AMR_class_total <- amr_melted_analytic[Level_ID=="Class", .(sum_class= sum(Normalized_Count)),by=.(Name)]
AMR_class_total[,total:= sum(sum_class)]
AMR_class_total[,percentage:= sum_class/total ,by=.(Name) ][order(sum_class)]

# Drop factors and set order with least abundant first
AMR_class_sum$Name = droplevels(AMR_class_sum$Name)
AMR_class_sum$Name = factor(AMR_class_sum$Name ,levels=c("Tunicamycin","Rifampin","Fosfomycin","Glycopeptides","Fluoroquinolones","Trimethoprim",
                                                                   "Bacitracin","Cationic antimicrobial peptides","MLS",
                                                                   "Phenicol","Sulfonamides","betalactams","Aminoglycosides","Multi-drug resistance" ,"Tetracyclines"))

# Add all other metadata to data object
AMR_class_sum$Group = factor(AMR_class_sum$Group, levels = c("Farm 1 Pre","Farm 1 Post","Farm 2 Pre", "Farm 2 Post"))
setkey(AMR_class_sum, ID)
setkey(metadata,ID)
AMR_class_sum <- AMR_class_sum[metadata]

# Change name of AMR class metadata column
AMR_class_sum$Class <- AMR_class_sum$Name
hist(AMR_class_sum[Name == 'Tetracyclines',percentage])
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_class_sum, aes(x = reorder(ID, Order), y = percentage, fill = Class)) + 
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Group , scales="free_x", ncol=4, strip.position = "bottom") +
  #scale_fill_brewer(palette="Dark2") +
  theme(strip.placement = "outside") + # code for changeing location of facet_wrap labels
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=16, angle=0),
    strip.text.y=element_text(size=16, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    axis.title.x = element_text(size = 20),
    legend.position="bottom",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  #ggtitle("Core resistome composition by collection time (only classes > 0.1%)") +
  xlab('Sample ID') +
  ylab('Mean relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1) 

ggsave("Figure1-Sample_AMR_composition_by_group.jpeg", width = 60, height = 30, units = "cm")

####################################################################
## Figure 2 - mean group relative abundance plots (AMR class) ##
####################################################################

# Sum class abundance by sample

AMR_class_mean_group <- AMR_class_sum[, .(mean_percentage_by_group = mean(percentage)), by = .(Name, Group)]
AMR_class_mean_group[, .(total = sum(mean_percentage_by_group)), by = .(Group) ]

#AMR_class_median <- AMR_class_sum[, .(median_percentage_by_group = median(round(percentage,digits = 3))), by = .(Name, Group)]
#AMR_class_median[, .(total = sum(median_percentage_by_group)), by = .(Group) ]



# Change order of factors 
AMR_class_mean_group$Group = factor(AMR_class_mean_group$Group, levels = c("Farm 1 Pre","Farm 1 Post","Farm 2 Pre", "Farm 2 Post"))


# Drop factors and set order with least abundant first
AMR_class_mean_group$Name = droplevels(AMR_class_mean_group$Name)
AMR_class_mean_group$Name = factor(AMR_class_mean_group$Name ,levels=c("Tunicamycin","Rifampin","Fosfomycin","Glycopeptides","Fluoroquinolones","Trimethoprim",
                                                         "Bacitracin","Cationic antimicrobial peptides","MLS",
                                                         "Phenicol","Sulfonamides","betalactams","Aminoglycosides","Multi-drug resistance" ,"Tetracyclines"))


# Change name of AMR class metadata column
AMR_class_mean_group$Class <- AMR_class_mean_group$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_class_mean_group, aes(x = Group, y = mean_percentage_by_group, fill = Class)) + 
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Group , scales="free_x", ncol=4, strip.position = "bottom") +
  #scale_fill_brewer(palette="Dark2") +
  theme(strip.placement = "outside") + # code for changeing location of facet_wrap labels
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=16, angle=0),
    strip.text.y=element_text(size=16, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    axis.title.x = element_blank(),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  #xlab('Sample ID') +
  ylab('Mean relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1) 

ggsave("Figure2-AMRclass_composition_by_group.jpeg", width = 60, height = 30, units = "cm")

####################################################################
## Figure 3 - mean group relative abundance plots (AMR mechanism) ##
####################################################################
#AMR_mech_sum <- amr_melted_analytic[Level_ID=="Mechanism", .(sum_mech= sum(Normalized_Count)),by=.(ID, Name,Group)]
#AMR_mech_sum[,total:= sum(sum_mech), by=.(ID)][,percentage:= sum_mech/total ,by=.(ID, Name) ]



# Sum class abundance by sample
AMR_mech_mean_group <- amr_melted_analytic[Level_ID=="Mechanism", .(sum_mech= sum(Normalized_Count)),by=.(ID, Name,Group, Group_label)][order(-Group )]
#AMR_mech_mean_group[,total:= sum(sum_mech), by=.(ID)]


total_mech <- amr_melted_analytic[Level_ID=="Mechanism", .(sum_mech= sum(Normalized_Count)),by=.(Name)]
total_mech <- total_mech[,total:= sum(sum_mech)][,percentage:= sum_mech/total ,by=.(Name) ][order(percentage)]

AMR_mech_mean_group_CORE <- total_mech[percentage > .01][,1]
# Drop factors not in use
AMR_mech_mean_group_CORE$Name = droplevels(AMR_mech_mean_group_CORE$Name)
AMR_mech_mean_group_CORE_names <- c("Aminoglycoside N-acetyltransferases","23S rRNA methyltransferases","Class C betalactamases","Penicillin binding protein","Aminoglycoside efflux pumps","Phenicol efflux pumps",
                                    "Aminoglycoside O-nucleotidyltransferases","Tetracycline resistance major facilitator superfamily MFS efflux pumps",
                                    "MDR regulator", "Aminoglycoside O-phosphotransferases", "Sulfonamide-resistant dihydropteroate synthases",
                                    "Class A betalactamases","Multi-drug efflux pumps","Tetracycline resistance ribosomal protection proteins")

setkey(AMR_mech_mean_group_CORE,Name)
setkey(AMR_mech_mean_group,Name)
AMR_mech_mean_group_core <- AMR_mech_mean_group[Name %in% AMR_mech_mean_group_CORE_names]
AMR_mech_mean_group_core[,total:= sum(sum_mech), by=.(ID)]
AMR_mech_mean_group_core[,percentage:= sum_mech/total ,by=.(ID, Name) ]
AMR_mech_mean_group_core <- AMR_mech_mean_group_core[, .(mean_percentage_by_group = mean(percentage)), by=.(Name,Group)]


AMR_mech_mean_group_rare <- AMR_mech_mean_group[!(Name %in% AMR_mech_mean_group_CORE_names)]
AMR_mech_mean_group_rare[,total:= sum(sum_mech), by=.(Group)]
AMR_mech_mean_group_rare[,percentage:= sum_mech/total ,by=.(Group, Name) ]
AMR_mech_mean_group_rare[,mean_percentage_by_group := mean(percentage), by=.(Group)]

# Change order of factors 
AMR_mech_mean_group_core$Group = factor(AMR_mech_mean_group_core$Group, levels = c("Farm 1 Pre","Farm 1 Post","Farm 2 Pre", "Farm 2 Post"))

AMR_mech_mean_group_core$Name <- factor(AMR_mech_mean_group_core$Name ,levels= c("Aminoglycoside N-acetyltransferases","23S rRNA methyltransferases","Class C betalactamases","Penicillin binding protein","Aminoglycoside efflux pumps","Phenicol efflux pumps",
                                                                                 "Aminoglycoside O-nucleotidyltransferases","Tetracycline resistance major facilitator superfamily MFS efflux pumps",
                                                                                 "MDR regulator", "Aminoglycoside O-phosphotransferases", "Sulfonamide-resistant dihydropteroate synthases",
                                                                                 "Class A betalactamases","Multi-drug efflux pumps","Tetracycline resistance ribosomal protection proteins"))


# Change name of AMR class metadata column
AMR_mech_mean_group_core$Mechanism <- AMR_mech_mean_group_core$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_mech_mean_group_core, aes(x = Group, y = mean_percentage_by_group, fill = Mechanism)) + 
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Group, scales="free_x", ncol=4, strip.position = "bottom") +
  #scale_fill_brewer(palette="Dark2") +
  theme(strip.placement = "outside") + # code for changeing location of facet_wrap labels
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=16, angle=0),
    strip.text.y=element_text(size=16, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    axis.title.x = element_blank(),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("Core resistome composition by collection time (only mechanisms > 1% abundant)") +
  #xlab('Sample ID') +
  ylab('Mean relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1) 

ggsave("Figure3-AMRmech_composition_by_group.jpeg", width = 60, height = 30, units = "cm")



#####
# Figure 4
#####


ggplot(AMR_class_sum, aes(x = Class, y = sum_class, fill = Class)) + 
  geom_boxplot() +
  facet_wrap( ~ Group , scales="free_x", strip.position = "bottom" ) +
  #scale_fill_brewer(palette="Dark2") +
  theme(strip.placement = "bottom", 
        strip.background = element_blank()) + # code for changeing location of facet_wrap labels
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=16, angle=0),
    strip.text.y=element_text(size=16, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    axis.title.x = element_blank(),
    legend.position="bottom",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  #ggtitle("Core resistome composition by collection time (only classes > 0.1%)") +
  # xlab('Sample ID') +
  ylab('Normalized counts') +
  scale_fill_tableau("Tableau 20", direction = -1) 
ggsave("Figure4-AMR_Class_count_distribution.jpeg", width = 40, height = 20, units = "cm")


#####
# Figure 5 - AMR class Median count stacked bar plot
#####

AMR_class_median <- AMR_class_sum[, .(median_percentage_by_group = median(round(percentage,digits = 3))), by = .(Name, Group)]
AMR_class_median[, .(total = sum(median_percentage_by_group)), by = .(Group) ]

AMR_class_median$Class <- AMR_class_median$Name

ggplot(AMR_class_median, aes(x = Group, y = median_percentage_by_group, fill = Class)) + 
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Group , scales="free_x", ncol=4, strip.position = "bottom") +
  #scale_fill_brewer(palette="Dark2") +
  theme(strip.placement = "outside") + # code for changeing location of facet_wrap labels
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=16, angle=0),
    strip.text.y=element_text(size=16, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    axis.title.x = element_blank(),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  #xlab('Sample ID') +
  ylab('Median relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1) 

ggsave("Figure5-AMRClass_Median_composition_by_group.jpeg", width = 60, height = 30, units = "cm")


#####
# Figure 5 - AMR mechanism Median count stacked bar plot
#####
AMR_mech_median_group_core <- AMR_mech_mean_group[Name %in% AMR_mech_mean_group_CORE_names]
AMR_mech_median_group_core[,total:= sum(sum_mech), by=.(ID)]
AMR_mech_median_group_core[,percentage:= sum_mech/total ,by=.(ID, Name) ]
AMR_mech_median_group_core <- AMR_mech_median_group_core[, .(median_percentage_by_group = median(percentage)), by=.(Name,Group)]

# Change order of factors 
AMR_mech_median_group_core$Group = factor(AMR_mech_median_group_core$Group, levels = c("Farm 1 Pre","Farm 1 Post","Farm 2 Pre", "Farm 2 Post"))

AMR_mech_median_group_core$Name <- factor(AMR_mech_median_group_core$Name ,levels= c("Aminoglycoside N-acetyltransferases","23S rRNA methyltransferases","Class C betalactamases","Penicillin binding protein","Aminoglycoside efflux pumps","Phenicol efflux pumps",
                                                                                 "Aminoglycoside O-nucleotidyltransferases","Tetracycline resistance major facilitator superfamily MFS efflux pumps",
                                                                                 "MDR regulator", "Aminoglycoside O-phosphotransferases", "Sulfonamide-resistant dihydropteroate synthases",
                                                                                 "Class A betalactamases","Multi-drug efflux pumps","Tetracycline resistance ribosomal protection proteins"))


# Change name of AMR class metadata column
AMR_mech_median_group_core$Mechanism <- AMR_mech_median_group_core$Name
#AMR_class_sum[,percentage:= round(sum_class/total, digits=2) ,by=.(ID, Name) ] removes some with low proportions
ggplot(AMR_mech_median_group_core, aes(x = Group, y = median_percentage_by_group, fill = Mechanism)) + 
  geom_bar(stat = "identity",colour = "black")+
  facet_wrap( ~ Group, scales="free_x", ncol=4, strip.position = "bottom") +
  #scale_fill_brewer(palette="Dark2") +
  theme(strip.placement = "outside") + # code for changeing location of facet_wrap labels
  theme(
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=16, angle=0),
    strip.text.y=element_text(size=16, angle=0),
    axis.text.x=element_blank(), #element_text(size=16, angle=20, hjust=1)
    axis.text.y=element_text(size=22),
    axis.title=element_text(size=26),
    axis.title.x = element_blank(),
    legend.position="right",
    panel.spacing=unit(0.1, "lines"),
    plot.title=element_text(size=32, hjust=0.5),
    legend.text=element_text(size=15),
    legend.title=element_text(size=20),
    panel.background = element_rect(fill = "white")
  ) +
  ggtitle("Core resistome composition by collection time (only mechanisms > 1% abundant)") +
  #xlab('Sample ID') +
  ylab('Mean relative abundance') +
  scale_fill_tableau("Tableau 20", direction = -1) 

ggsave("Figure6-AMRmech_Median_composition_by_group.jpeg", width = 60, height = 30, units = "cm")
