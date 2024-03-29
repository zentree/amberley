# Creating pedigree file for Amberley
# Luis A. Apiolaza
# School of Forestry
# University of Canterbury
#

# Changing default conversion of strings to factors
# (source of extreme pain when debugging)
options(stringsAsFactors = FALSE)

#### Setting Working directory
# This directory has to be changed to reflect
# the location of the files in your computer
setwd('~/Documents/Research/2012/amberley/')


# Reading pedigree of families used in Amberley
# and then matching with field codes
amb.codes <- read.table('ShafCodes.txt', header = TRUE)

amb.codes$Mum <- paste(amb.codes$FSERIES, amb.codes$FCLONE, sep = '-')
amb.codes$Dad <- paste(amb.codes$MSERIES, amb.codes$MCLONE, sep = '-')

# We only need the pedigree information from this dataset
ped .info <- data.frame(Family = amb.codes$Family, 
                       Mum = amb.codes$Mum, 
                       Dad = amb.codes$Dad)

# There is a couple of family codes in the trial that were not part of the 
# original group of seed offered and not contained in ShafCodes.txt 
# (Confirmed on 2 June 2009)
miss.codes <- data.frame(Family = c('06-520', '06-524'), 
                         Mum = c('268-262','875-066'), 
                         Dad = c('268-248','875-076'))

ped.info <- rbind(ped.info, miss.codes)

# Up to this point we have Family codes, Mum code and Dad code
# for each Family. We need to create records in the pedigree
# specifying unknown parents for each Mum and Dad (the base
# population)
base.ids <- unique(c(ped.info$Mum, ped.info$Dad))
pedigree <- data.frame(tree.id = base.ids, 
                       Mum = rep(NA, length(base.ids)), 
                       Dad = rep(NA, length(base.ids)))

# Once we have the pedigree for the base population we can move
# to the Amberley trial trees. This is an example with the
# 2009 assessment
amber <- read.csv('data20091016.csv', header = TRUE)
amber$Block <- factor(amber$Block)
summary(amber)

# Here we generate a unique tree identifier (to use in the genetic
# evaluation)
amber$tree.id <- 10001:(10000+dim(amber)[1])
head(amber)

# Joining trial data to family info
amber <- merge(amber, ped.info, by='Family', sort = FALSE, all.x = TRUE)
trial.ped <- data.frame(tree.id = amber$tree.id, 
                        Mum = amber$Mum, 
                        Dad = amber$Dad)

pedigree <- rbind(pedigree, trial.ped)

# Only if we wanted to save pedigree
# write.csv(pedigree, file = 'pedigree.csv', quote = FALSE, row.names = FALSE)
