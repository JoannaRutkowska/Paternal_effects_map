#Title:  "Mapping the past, present and future research landscape of paternal effects"
#Author: "Joanna Rutkowska, Malgorzata Lagisz, Russell Bonduriansky, Shinichi Nakagawa"
#Date:   2 May 2020

# prepare 

#sessionInfo()
#R version 3.6.0 (2019-04-26)
#Platform: x86_64-apple-darwin15.6.0 (64-bit)
#Running under: macOS Sierra 10.12.6

rm(list=ls()) #clean up

remove.packages("tibble")

library(tidyverse) #tidy family and related pacakges below
library(magrittr) #extending piping
library(DescTools) 
library(fulltext) #
library(readxl) #reading in Excel files
library(janitor) #data cleaning
library(bibliometrix) #bibliometric analyses
library(ggplot2) 
library(ggbeeswarm) #making bee-swarm plots possible
library(plotly)    #interactive plots using ggplot2
library(easyalluvial) # allovial plots


#  Main data files

empirical <- read_excel("2019-07-11-fulltext-empirical-included.xlsx")
nonempirical <- read_excel("2019-05-25-fulltext-codes-nonempirical-included.xlsx")

#change variable names in the non-empirical file
nonempirical<- nonempirical %>% rename(
  Doi = "copy & paste doi",
  Form = "What publication form does the record have?", 
  Type = "What type of paper is it or claims to be?", 
  Taxon = "What is its taxonomic scope?", 
  Primary_focus = "What is its primary focus?", 
  Secondary_focus = "What is its secondary focus?")

colnames(nonempirical)[12] <- "primary_other"
colnames(nonempirical)[14] <- "secondary_other"
names(nonempirical)


# Cluster analysis (from VosViewer), after Nettle 2019

#Procedure:
# - Read in map and cluster files created by VOSViewer.
# - Records included in our systematic review were searched forin Scopus database on 16th September 2019. 
# - Together with their citaions, they were uplodaed to VOSViewer. 
# - VOSViewer indietified that the records are grouped in 3 clusters and created two files (network and map) which are uplodaed here.

network <- read.delim(file = "2019-09-16-3cl-networktrial.txt", header = TRUE, sep = "\t", dec = ".")
#str(network)
network %>% rename(from.id = X1, to.id = X8 , links = X0.0227) -> network #rename columns

map <- read.delim("2019-09-16-3cl-maptrial.txt")
#str(map)
#map %>% rename(id = ď.żid) -> map #rename columns

#Make matrix of connections
clusters <- select(map, id, cluster)
str(clusters)

#Merge and reshape data frames
d <- merge(network, clusters, by.x="from.id", by.y="id")
d <- merge(d, clusters, by.x="to.id", by.y="id")
parta <- xtabs(~d$cluster.x + d$cluster.y) 
partb <- t(parta) #transpose
link.matrix <- parta + partb
xtabs(~map$cluster)

#Get numbers of papers per cluster
n <- as.vector(xtabs(~map$cluster))

#Calculate expected number of links within and between clusters
expected <- n%*%t(n)
diag(expected) <- diag(expected)-n
expected
#calculate ratio of actual numbers of links to expected number of links
f <- link.matrix/expected

#Make data frame of the connection indices
x <- c(rep("A1", times=3), rep("A2", times=3), rep("A3", times=3))
y <- rep(c("A1", "A2", "A3"), times=3)
z <- as.vector(f)
l <- data.frame(x, y, z)

#Assign names to cluster numbers
l$x <- as.character(l$x)
l$x[l$x == "A1"] <- "Med"
l$x[l$x == "A2"] <- "Tox"
l$x[l$x == "A3"] <- "EcoEvo"
l$y <- as.character(l$y)
l$y[l$y == "A1"] <- "Med"
l$y[l$y == "A2"] <- "Tox"
l$y[l$y == "A3"] <- "EcoEvo"

l$y <- factor(l$y, levels = c('Med','Tox','EcoEvo'))
l$x <- factor(l$x, levels = c('Med','Tox','EcoEvo'))

#Plot figure showing connections between the clusters
figure=ggplot(l, aes(x=x, y=z, fill=x)) + 
  geom_bar(stat="identity") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  facet_grid(~y) + 
  xlab("Cluster") + 
  ylab("Connection index") + 
  scale_fill_discrete(name="Cluster") + 
  coord_cartesian(ylim=c(0, 0.5)) +
  theme(aspect.ratio=0.7) +
  ggsave("Fig-3b.pdf")


# Merge data files with bibliometric clustering

#Preparing empirical data file for merging the records based on Doi
empirical$Doi <- tolower(empirical$Doi)
empirical$Doi <- gsub("/", ".", empirical$Doi)
empirical$Doi <- gsub("-", ".", empirical$Doi) 
empirical$Doi <- gsub("–", ".", empirical$Doi) 
empirical$Doi <- gsub("_", ".", empirical$Doi) 
empirical$Title <- tolower(empirical$Title) 

#Puting title in the place of Doi (if there is no Doi)  in a new columen called link
empirical <- mutate(empirical, link = Doi)
empirical <- mutate(empirical,link = ifelse(Doi == "..", Title, link))

#Preparing map data (from VOSViewer) for merging records based on Doi 
map$url <- gsub("https://doi.org/", "", map$url)

map %>% rename(Doi = url) -> map
map$Doi <- tolower(map$Doi)
map$Doi <- gsub("/", ".", map$Doi)
map$Doi <- gsub("-", ".", map$Doi) 
map$Doi <- gsub("–", ".", map$Doi) 
map$Doi <- gsub("_", ".", map$Doi) 
map$description <- gsub(".*</td></tr><tr><td>Title:</td><td>(.+)</td></tr><tr><td>Source:</td><td>.*", "\\1", map$description)
map %>% rename(Title = description) -> map

#Putting title in the place of Doi (if there is no Doi) in a new columen called link
map <- mutate(map, link = Doi)
map <- mutate(map,link = ifelse(Doi == "", Title, link))


#Adding cluster ID to the empirical layer
empirical_clusters <- left_join(empirical, map, by="link") #merging empirical data file with map file by doi 
write.csv(empirical_clusters, 'empirical_clusters.csv') #saving empirical_clusters to a .csv file
clusters_assigned <- data.frame(xtabs(~cluster, empirical_clusters)) #separating empirical_clusters that are assigned to one of the 3 clusters
nocluster <- filter(empirical_clusters, is.na(empirical_clusters$cluster)) #identifying records not assigned to clusters - needs fixing or removing
write.csv(nocluster, 'nocluster.csv') #saving unassigned records in a .csv file


# Bibilometric analysis aimimg at presenting number of citations 

# - Number of citaions is visualised for all paperes included in the systematic review, Thus, it is not restriced to the paperes assigend to clusters. 
# - Records included in our systematic review were searched for in Scopus database on 23rd July 2019. 
# - This was done separatly for empircial and non-empirical paperes. 
# - In both cases, the papers, together with their citations, were uplodaed here. 

#Getting the data from Scopus database for empirical layer
bibE <- convert2df(file = "2019-07-23-empiricalBib.bib", dbsource = "scopus", format = "bibtex")
names(bibE)
write.csv(bibE, "bibE_as_df.csv", row.names = FALSE)
#str(bibE)

#Getting the data from Scopus database for non-empirical layer
bibNE <- convert2df(file = "2019-07-23-nonempiricalBib.bib", dbsource = "scopus", format = "bibtex")
names(bibNE)
write.csv(bibNE, "bibNE_as_df.csv", row.names = FALSE)
#str(bibNE)

#For both layers, adding additional columns with layer label
bibE <-  mutate(bibE, map_layer = "empirical")
bibE <- mutate(bibE, BN = "NA")
bibNE <- mutate(bibNE, map_layer = "non_empirical")
names(bibE)
setdiff(colnames(bibE), colnames(bibNE))

#Dropping columns that are only in one of the data frames
bibE %>% select (-c(tradenames, manufacturers, molecular_seqnumbers)) -> bibE

#Joining data of both layers
bibE_NE <- union(bibE, bibNE)
#str(bibE_NE)

#Cleaning records that cause trouble
grep("\\(NEOTIGASON/SORIATANE?", bibE_NE$TI)
bibE_NE$TI <- gsub("\\(NEOTIGASON/SORIATANE?","(NEOTIGASON/SORIATANE)?", bibE_NE$TI)
grep("NEOTIGASON/SORIATANE", bibE_NE$CR)
bibE_NE$CR[37]
bibE_NE$CR <- gsub("\\(NEOTIGASON/SORIATANE","(NEOTIGASON/SORIATANE)", bibE_NE$CR)

write.csv(bibE_NE, "bibE_NE.csv", row.names = FALSE) #saving in a .csv file
#bibE_NE <- read.csv("bibE_NE.csv")


# Merging the file with cluster assigment, making figure}

#Preparing for joining the file based on Doi - cleaning DOI codes
bibE_NE$Doi <- gsub("https://doi.org/", "", bibE_NE$DI)
bibE_NE$Doi <- tolower(bibE_NE$Doi)
bibE_NE$Doi <- gsub("/", ".", bibE_NE$Doi)
bibE_NE$Doi <- gsub("-", ".", bibE_NE$Doi) 
bibE_NE$Doi <- gsub("_", ".", bibE_NE$Doi) 

#Joining based on Doi
joined <- left_join(bibE_NE, map, by = 'Doi')

#Preparing data for the figure
joined$cluster <- factor(joined$cluster)
joined$cluster <- as.character(joined$cluster)
joined$cluster[joined$cluster == "1"] <- "Med"
joined$cluster[joined$cluster == "2"] <- "Tox"
joined$cluster[joined$cluster == "3"] <- "EcoEvo"
joined$cluster <- factor(joined$cluster, levels = c('Med','Tox','EcoEvo'))


# Making figure showing total citation count (TC) according to publication year (PY) and map layer and cluster. 

#Paperes not assigned to one of three clusters are in grey.
ggplot(joined, aes(x=PY, y=TC, color = cluster, shape = map_layer))+
  geom_point(size=4, alpha=0.6)+
  scale_size(guide="none")+
  scale_alpha(guide="none")+
  scale_x_continuous(breaks=c(1970,1980,1990,2000, 2010, 2020))+
  theme_classic()+
  labs(x= "Publication year", y = "Total citations")+
  theme_classic()+ 
  theme(legend.position = "none") +
  theme(aspect.ratio=0.45)+
  ggsave("Fig-3c.pdf")  


# Visulasiation of the empirical layer of the map as an alluvial plot

#Preparing the data by separation multiply exposure into separate rows of data
empirical_clusters_expo <- empirical_clusters %>% separate_rows(Exposure_category, sep = "; ")

#Preparing the data by selecting the most important variables to be used in the plot
empirical_clusters_expo1 <- select(empirical_clusters_expo, Taxon, Source_of_population, Exposure_category, Maternal_exposure, Offspring_exposure, cluster)

#Changing duplicated names of categories
empirical_clusters_expo1$Taxon <- as.character(empirical_clusters_expo1$Taxon)
empirical_clusters_expo1$Taxon[empirical_clusters_expo1$Taxon == "human"] <- "men"
empirical_clusters_expo1$Maternal_exposure <- as.character(empirical_clusters_expo1$Maternal_exposure)
empirical_clusters_expo1$Maternal_exposure[empirical_clusters_expo1$Maternal_exposure == "no"] <- "none"

#Separating into clusters and changing columns into factors
clusterMed <- subset(empirical_clusters_expo1, empirical_clusters_expo1$cluster=="1")
clusterTox <- subset(empirical_clusters_expo1, empirical_clusters_expo1$cluster=="2")
clusterEE <- subset(empirical_clusters_expo1, empirical_clusters_expo1$cluster=="3")
clusterMed$Taxon <-as.factor(clusterMed$Taxon)
clusterMed$Source_of_population <- as.factor(clusterMed$Source_of_population)
clusterMed$Exposure_category <- as.factor(clusterMed$Exposure_category)
clusterMed$Maternal_exposure <- as.factor(clusterMed$Maternal_exposure)
clusterMed$Offspring_exposure <- as.factor(clusterMed$Offspring_exposure)
clusterMed$cluster <- as.factor(clusterMed$cluster)
clusterTox$Taxon <- as.factor(clusterTox$Taxon)
clusterTox$Source_of_population <- as.factor(clusterTox$Source_of_population)
clusterTox$Exposure_category <- as.factor(clusterTox$Exposure_category)
clusterTox$Maternal_exposure <- as.factor(clusterTox$Maternal_exposure)
clusterTox$Offspring_exposure <- as.factor(clusterTox$Offspring_exposure)
clusterTox$cluster <- as.factor(clusterTox$cluster)
clusterEE$Taxon <- as.factor(clusterEE$Taxon)
clusterEE$Source_of_population <- as.factor(clusterEE$Source_of_population)
clusterEE$Exposure_category <- as.factor(clusterEE$Exposure_category)
clusterEE$Maternal_exposure <- as.factor(clusterEE$Maternal_exposure)
clusterEE$Offspring_exposure <- as.factor(clusterEE$Offspring_exposure)
clusterEE$cluster <- as.factor(clusterEE$cluster)

#Setting color pallets for each cluster for Figure 4a (fllows and starata have separate pallets)
col_vectorMed3a = c('Abiotic habitat'="#6b452b",'Age'="#6fa6ab", 'Alcohol'="#b0b54e", 'Chemical substance'="#7081a1",  'Diet'="#edcb8a", 'Drug'="#70697d",'Physiological factor'="#a2a39b",'Psychological factor'="#1b1b1c" )
col_vectorMed3b = c('human'="#828585",'domesticated'="gray55",'captive'="gray10", 'wild'="gray60", 'men'="#828585",'non-human mammal'="gray5", 'fish'="gray15", 'other vertebrate'="gray25", 'other invertebrate'="gray35", 'plant'="gray45", 'Abiotic habitat'="#6b452b",'Age'="#6fa6ab",  'Alcohol'="#b0b54e", 'Chemical substance'= "#7081a1",  'Diet'="#edcb8a", 'Drug'="#70697d",'Physiological factor'="#a2a39b",'Psychological factor'="#1b1b1c")
col_vectorTox3a = c('Abiotic habitat'="#6b452b",'Chemical substance'="#7081a1",  'Diet'="#edcb8a", 'Drug'="#70697d", 'Physiological factor'="#a2a39b")
col_vectorTox3b = c('human'="#828585", 'captive'="gray10", 'wild brought into captivity'="gray50",'men'="#828585",'non-human mammal'="gray5", 'arthropod'="gray40", 'Abiotic habitat'="#6b452b",'Chemical substance'="#7081a1",  'Diet'="#edcb8a", 'Drug'="#70697d", 'Physiological factor'="#a2a39b")
col_vectorEE3a = c('Abiotic habitat'="#6b452b",'Age'="#6fa6ab", 'Chemical substance'="#7081a1", 'Diet'="#edcb8a", 'Past experience'="#3a3e3f", 'Physiological factor'="#a2a39b",'Psychological factor'="#1b1b1c")
col_vectorEE3b = c('captive'="gray10",  'wild brought into captivity'="gray50",'wild'="gray60", 'non-human mammal'="gray5", 'bird'="gray30" , 'fish'="gray15", 'other vertebrate'="gray25",'arthropod'="gray60" , 'other invertebrate'="gray35", 'plant'="gray45",'Abiotic habitat'="#6b452b",'Age'="#6fa6ab", 'Chemical substance'="#7081a1", 'Diet'="#edcb8a", 'Past experience'="#3a3e3f", 'Physiological factor'="#a2a39b",'Psychological factor'="#1b1b1c")

#Setting color palets for each cluster for figure 4b 
col_vectorMed4 = c('Abiotic habitat'="#6b452b",'Age'="#6fa6ab",  'Alcohol'="#b0b54e", 'Chemical substance'="#7081a1",  'Diet'="#edcb8a", 'Drug'="#70697d",'Physiological factor'="#a2a39b",'Psychological factor'="#1b1b1c", 'none'="black", 'independent'="#828585", 'diallel'="gray", 'NA'="white", 'no'="black", 'yes'="gray")
col_vectorTox4 = c('Abiotic habitat'="#6b452b",'Chemical substance'="#7081a1",  'Diet'="#edcb8a", 'Drug'="#70697d", 'Physiological factor'="#a2a39b", 'none'="black", 'independent'="#828585", 'diallel'="gray", 'no'="black")
col_vectorEE4 = c('Abiotic habitat'="#6b452b",'Age'="#6fa6ab", 'Chemical substance'="#7081a1", 'Diet'="#edcb8a", 'Past experience'="#3a3e3f", 'Physiological factor'="#a2a39b",'Psychological factor'="#1b1b1c", 'none'="black", 'independent'="#828585", 'diallel'="gray", 'NA'="white", 'no'="black", 'yes'="gray")


#Graphs for figure 4a
allMed3 <- alluvial_wide(select(clusterMed, Source_of_population, Taxon, Exposure_category ), fill_by = 'last_variable', order_levels = c('human','domesticated','captive', 'wild', 'men','non-human mammal', 'fish', 'other vertebrate', 'other invertebrate', 'plant','Abiotic habitat','Age',  'Alcohol', 'Chemical substance',  'Diet', 'Drug','Physiological factor','Psychological factor'),stratum_label_size = 7, 
                        col_vector_flow = col_vectorMed3a, col_vector_value = col_vectorMed3b) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0,0.2)) +
  labs(title= 'Med cluster', y = "Frequency") +
  theme(plot.title = element_text(size = 20, face="bold")) +
  ggsave("Fig-Med-4a.pdf", scale = 6, width = 10, units=c("cm") )

allTox3 <- alluvial_wide(select(clusterTox, Source_of_population, Taxon, Exposure_category ), fill_by = 'last_variable',   order_levels = c('human', 'captive', 'wild brought into captivity','men','non-human mammal', 'arthropod' ), stratum_label_size = 7,
                         col_vector_flow = col_vectorTox3a, col_vector_value = col_vectorTox3b) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0,0.2)) +
  labs(title= 'Tox cluster', y = "Frequency") +
  theme(plot.title = element_text(size = 20, face="bold")) +
  ggsave("Fig-Tox-4a.pdf", scale = 6, width = 10, units=c("cm"))

allEE3 <- alluvial_wide(select(clusterEE,Source_of_population, Taxon, Exposure_category ), fill_by = 'last_variable',  order_levels = c('captive',  'wild brought into captivity','wild','non-human mammal', 'bird', 'fish', 'other vertebrate','arthropod', 'other invertebrate', 'plant', 'Abiotic habitat','Age','Chemical substance',  'Diet', 'Past experience', 'Physiological factor','Psychological factor'), stratum_label_size = 7,
                        col_vector_flow = col_vectorEE3a, col_vector_value = col_vectorEE3b) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0,0.2)) +
  labs(title= 'EcoEvo cluster', y = "Frequency") +
  theme(plot.title = element_text(size = 20, face="bold")) +
  ggsave("Fig-EE-4a.pdf", scale = 6, width = 10, units=c("cm") )


# Graphs for figure 4b
allMed4 <- alluvial_wide(select(clusterMed, Exposure_category, Maternal_exposure, Offspring_exposure ), fill_by = 'first_variable', order_levels = c('Abiotic habitat','Age',  'Alcohol', 'Chemical substance',  'Diet', 'Drug','Physiological factor','Psychological factor', 'none' , 'independent' , 'diallel' , 'NA' , 'no' , 'yes' ), stratum_label_size = 7,
                        col_vector_flow = col_vectorMed4, col_vector_value = col_vectorMed4) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0,0.2)) +
  labs(title= 'Med cluster', y = "Frequency") +
  theme(plot.title = element_text(size = 20, face="bold")) +
  ggsave("Fig-Med-4b.pdf", scale = 6, width = 10, units=c("cm") )

allTox4 <- alluvial_wide(select(clusterTox, Exposure_category, Maternal_exposure, Offspring_exposure), fill_by = 'first_variable',   order_levels = c('Abiotic habitat', 'Chemical substance',  'Diet', 'Drug', 'Physiological factor', 'none', 'independent', 'diallel' , 'no' ), stratum_label_size = 7,
                         col_vector_flow = col_vectorTox4, col_vector_value = col_vectorTox4) + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0,0.2)) +
  labs(title= 'Tox cluster', y = "Frequency") +
  theme(plot.title = element_text(size = 20, face="bold")) +
  ggsave("Fig-Tox-4b.pdf", scale = 6, width = 10, units=c("cm"))

allEE4 <- alluvial_wide(select(clusterEE,Exposure_category, Maternal_exposure, Offspring_exposure), fill_by = 'first_variable',  order_levels = c('Abiotic habitat','Age', 'Chemical substance', 'Diet', 'Past experience', 'Physiological factor','Psychological factor', 'none' , 'independent' , 'diallel' , 'NA', 'no', 'yes'), stratum_label_size = 7,
                        col_vector_flow = col_vectorEE4, col_vector_value = col_vectorEE4) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  scale_x_discrete(expand = c(0,0.2)) +
  labs(title= 'EcoEvo cluster', y = "Frequency") +
  theme(plot.title = element_text(size = 20, face="bold")) +
  ggsave("Fig-EE-4b.pdf", scale = 6, width = 10, units=c("cm") )


 
# Timeline of emprical studies

#Preparing data for the timeline figure
freq3 <- data.frame(xtabs(~Year + Exposure_category, empirical_clusters_expo))
freq3$Year <- as.numeric(as.character(freq3$Year))

#Plotting the timeline figure
ggplot(freq3, aes(Year, Freq)) +
  geom_area(aes(fill = Exposure_category)) +
  scale_fill_manual(values = c("#6b452b", "#6fa6ab","#b0b54e", "#7081a1", "#edcb8a", "#70697d", "#3a3e3f", "#a2a39b",  "#1b1b1c","#a89985")) +
  scale_y_continuous(name="Frequency") +
  scale_x_continuous(breaks=c(1970, 1980, 1990, 2000, 2010, 2020)) +
  theme_classic() +
  theme(aspect.ratio=0.4) +
  ggsave("Fig-2a.pdf")  


# Visulasiation of the non-empirical layer of the map 

#Timeline for non-empirical studies
#Preparing data for timeline figure
freq2 <- data.frame(xtabs(~Year + Type, nonempirical)) 
freq2$Year <- as.numeric(as.character(freq2$Year))

#Removing values with zeros 
freq23 <- filter(freq2, Freq !=0)

#Adding zeros at the begining of time series
first_years <- group_by(freq2, Type) %>%
  summarise(Year = min(Year) - 1) %>%
  filter(Year > 1970) %>%
  mutate(Freq = 0)

freq24 <- bind_rows(freq2, first_years) %>% arrange(Year, Type)
freq24$Type <- factor(freq24$Type, levels = c('commentary-perspective', 'narrative review','systematic review family','theoretical paper'))

#Plotting the timeline of non-empirical papers
ggplot(freq24, aes(x=Year, y=Freq,fill = Type) ) +
  geom_area() +
  theme_classic() +
  scale_fill_manual(values = c('theoretical paper' = "#141414", 'systematic review family' = "#828585" ,'narrative review' = "gray", 'commentary-perspective' = "#4f5157" )) +
  expand_limits(x = c(1960, 2020)) +
  scale_x_continuous(breaks=c(1970,1980,1990,2000, 2010, 2020)) +
  scale_y_continuous(breaks=c(0,10,20), name="Frequency") +
  theme(aspect.ratio=0.2) +
  ggsave("Fig-2b.pdf")


# Visualising focus of non-empircal studies

#Preparing data for visualising focus of the studies
nonempirical_focus <- data.frame(xtabs(~Taxon + Primary_focus + Secondary_focus, nonempirical))
nonempirical_focus <- filter(nonempirical_focus, Freq !=0)
nonempirical_focus$Taxon <- factor(nonempirical_focus$Taxon, levels = c('animals and plants', 'animals', 'mammals', 'humans'))
nonempirical_focus <- nonempirical_focus %>% arrange(Taxon, Primary_focus, -Freq)
nonempirical_focus <- nonempirical_focus %>% mutate(proportion = Freq/sum(Freq)) 
nonempirical_focus <- nonempirical_focus %>% mutate(ymax = cumsum(proportion)) 
nonempirical_focus <- nonempirical_focus %>% mutate(ymin = ifelse(ymax!=0, lag(ymax), '0')) 
nonempirical_focus$ymin[is.na(nonempirical_focus$ymin)] <- 0

#Saving the transformed data
write.csv(nonempirical_focus, 'nonempirical_focus.csv')

#Making figure
ggplot(nonempirical_focus) + 
  geom_rect(aes(fill=Primary_focus, ymax=ymax, ymin=ymin, xmax=2.9, xmin=0.2)) +
  geom_rect(aes(fill=Secondary_focus, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
  xlim(c(0, 4)) + 
  theme(aspect.ratio=0.7) +
  coord_polar(theta="y")  +
  facet_wrap(~ Taxon, nrow = 1) +
  scale_fill_manual(values = c('age'="#6fa6ab",'alcohol'="#b0b54e", 'metabolic disorders'="#edcb8a", 'assisted reproduction techniques'= "#c03728",'effects of drugs/toxins'= "#7081a1", 'ecology & evolution'= "#e2d51d", 'general'= "#828585",  'none'= "#3a3e3f",'offspring cancer'="#e68c7c",'proximate mechanisms'="#c83e74",'other'="#daa03d" )) +
  theme_void() +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  ggsave("Fig-2c.pdf")


# Calculations of %s of empirical clusters for the main text

#Preparing empirical clusters
empirical_clusters$cluster <- as.character(empirical_clusters$cluster)
empirical_clusters$cluster[empirical_clusters$cluster == "1"] <- "Med"
empirical_clusters$cluster[empirical_clusters$cluster == "2"] <- "Tox"
empirical_clusters$cluster[empirical_clusters$cluster == "3"] <- "EcoEvo"

#Stage of male exposure
empirical_clusters_male_stage <- empirical_clusters %>% separate_rows(Exposure_stage, sep = ", ") #Separating males stages into rows
freq_male_stage <- data.frame(xtabs(~Exposure_stage + cluster, empirical_clusters_male_stage)) 
freq4 <- freq_male_stage %>% group_by(cluster) %>% mutate(prop=Freq/sum(Freq))
freq4

#Mate choice across the field 
freq_mate_choice <- data.frame(xtabs(~Mating_in_study, empirical_clusters)) 
freq_mate_choice <- freq_mate_choice %>% mutate(proportion = Freq/sum(Freq))

#Maternal effects per cluster
freq_maternal_effects <- data.frame(xtabs(~Maternal_effects + cluster, empirical_clusters)) 
freq6 <- freq_maternal_effects %>% group_by(cluster) %>% mutate(prop=Freq/sum(Freq))
freq6

#Maternal exposure per cluster
freq_maternal_exposure <- data.frame(xtabs(~Maternal_exposure + cluster, empirical_clusters)) 
freq7 <- freq_maternal_exposure %>% group_by(cluster) %>% mutate(prop=Freq/sum(Freq))
freq7 <- freq7 %>% arrange(Maternal_exposure, cluster)
freq7

#Offspring studied until adulthood per cluster
freq_offspring_last_stage <- data.frame(xtabs(~Latest_offspring_stage + cluster, empirical_clusters)) 
freq8 <- freq_offspring_last_stage %>% group_by(cluster) %>% mutate(prop=Freq/sum(Freq))
freq8 <- freq8 %>% arrange(Latest_offspring_stage, cluster)
freq8

#Accounting for offspring sex per cluster
freq_offspring_sex <- data.frame(xtabs(~Sex_specific_effects + cluster, empirical_clusters)) 
freq9 <- freq_offspring_sex %>% group_by(cluster) %>% mutate(prop=Freq/sum(Freq))
freq9 <- freq9 %>% arrange(cluster, Sex_specific_effects)
freq9

#Studying grand-offspring across the field 
freq_grand_offspring <- data.frame(xtabs(~Grand_offspring_traits, empirical_clusters)) 
freq_grand_offspring <- freq_grand_offspring %>% mutate(proportion = Freq/sum(Freq))
freq_grand_offspring

#Offspring exposure per cluster
freq_offspring_exposure <- data.frame(xtabs(~Offspring_exposure + cluster, empirical_clusters)) 
freq10 <- freq_offspring_exposure %>% group_by(cluster) %>% mutate(prop=Freq/sum(Freq))
freq10 <- freq10 %>% arrange(cluster, Offspring_exposure)
freq10


# Additional notes
# - The figure with visualisation of the clusters (Figure3a) was created in VosViewer
# - All figure panels were assembled ourside R environment and edited for clarity and presentation.


# Create a resaerch compendium with holepunch

#holepunch makes a binder instance from this R project

#install and load holepunch
devtools::install_github("karthik/holepunch")
library(holepunch)
write_compendium_description(package = "Paternal_effects_map", description = "Mapping the past present and future research landscape of paternal effects") #create compedium description file
write_dockerfile(maintainer = "ML") #Dockerfile will automatically pick the date of the last modified file, match it to that version of R and add it
generate_badge() #generate a badge for the readme file

# push the code to GitHub 
# click on the badge or use the function below to get the build ready ahead of time.
 


write_install() # Writes install.R with all your dependencies
write_runtime() # Writes the date your code was last modified. Can be overridden.
generate_badge() # Generates a badge you can add to your README. Clicking badge will launch the Binder.
# ----------------------------------------------
# push the code to GitHub
# ----------------------------------------------
# Then click the badge on your README or run
build_binder() # to kick off the build process


# #install.R
# write_install()
# #runtime.txt
# write_runtime()
# generate_badge() #add a link on the readme of the github readme to the binder instance 

