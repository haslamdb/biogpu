setwd("C:/Users/dbhas/Downloads")
list.files()

genes <- read.table("fluoroquinolone_ressitance_table.tsv", 
                 header = TRUE, 
                 sep = "\t", 
                 comment.char = "#", 
                 fill = TRUE, 
                 quote = "")
names(genes)<-c("Organism", "Protein", "Project", "Isolate", "Contig", "start", "stop", "Strand", "Element_symbol", "Element_name", "Type", "Scope", "Subtype", "Class", "Subtype", "Method", "Coverage", "Coverageb" )
genes$Coverage<-as.numeric(genes$Coverage)
# genes<-subset(genes, genes$Coverage >=98)
genes<-subset(genes, ! genes$Protein == "")
genes <- genes[order(genes$Coverage, decreasing = TRUE), ]

genes<-genes[,c(1:16)]
genes$start <- gsub("X", "", genes$start)
genes$stop <- gsub("X", "", genes$stop)
genes$Strand<-ifelse(genes$start >= genes$stop, "-", "+")
genes$Org_Symbol<-paste(genes$Organism, genes$Element_symbol, sep = "-")
genes$Name_Symbol<-paste(genes$Element_name, genes$Element_symbol, sep = "-")
genes<-na.omit(genes)



unique_mutations_table<-subset(genes, ! duplicated(genes$Name_Symbol))
unique_mutations_table$gene <- sapply(strsplit(unique_mutations_table$Element_symbol, "_"), `[`, 1)
unique_mutations_table$mutation <- sapply(strsplit(unique_mutations_table$Element_symbol, "_"), `[`, 2)
unique_mutations_table$mutation<-gsub("-", "", unique_mutations_table$mutation)
unique_mutations_table$wt <- sub("^([A-Z])[0-9]+[A-Z]$", "\\1", unique_mutations_table$mutation)
unique_mutations_table$location <- as.integer(sub("^[A-Z]([0-9]+)[A-Z]$", "\\1", unique_mutations_table$mutation))
unique_mutations_table$mut <- sub("^[A-Z][0-9]+([A-Z])$", "\\1", unique_mutations_table$mutation)
unique_mutations_table$species <- sapply(strsplit(unique_mutations_table$Organism, " "), function(x) paste(x[1:2], collapse = " "))

mutations_table<-unique_mutations_table[,c(24, 19, 22, 21, 23, 2, 6, 7)]

write.csv(mutations_table, file = "Quinolone_resistance_mutation_table.csv")

allels<-read.table("fluoroquinolone_resistance_alleles_table.tsv", 
                   header = TRUE, 
                   sep = "\t", 
                   comment.char = "#", 
                   fill = TRUE, 
                   quote = "")
names(allels)<-c("Organism", "Protein", "Project", "Isolate", "Contig", "start", "stop", "Strand", "Element_symbol", "Element_name", "Type", "Scope", "Subtype", "Class", "Subtype", "Method", "Coverage", "Coverageb" )
allels$Coverage<-as.numeric(allels$Coverage)
# allels<-subset(allels, allels$Coverage >=98)
allels<-subset(allels, ! allels$Protein == "")
allels <- allels[order(allels$Coverage, decreasing = TRUE), ]


allels<-allels[,c(1:16)]
allels$start <- gsub("X", "", allels$start)
allels$stop <- gsub("X", "", allels$stop)
allels$Strand<-ifelse(allels$start >= allels$stop, "-", "+")
allels$Org_Symbol<-paste(allels$Organism, allels$Element_symbol, sep = "-")
allels$Name_Symbol<-paste(allels$Element_name, allels$Element_symbol, sep = "-")
allels<-na.omit(allels)

unique_allels_table<-subset(allels, ! duplicated(allels$Name_Symbol))
unique_allels_table$gene <- sapply(strsplit(unique_allels_table$Element_symbol, "_"), `[`, 1)
unique_allels_table$mutation <- sapply(strsplit(unique_allels_table$Element_symbol, "_"), `[`, 2)
unique_allels_table$mutation<-gsub("-", "", unique_allels_table$mutation)
unique_allels_table$wt <- sub("^([A-Z])[0-9]+[A-Z]$", "\\1", unique_allels_table$mutation)
unique_allels_table$location <- as.integer(sub("^[A-Z]([0-9]+)[A-Z]$", "\\1", unique_allels_table$mutation))
unique_allels_table$mut <- sub("^[A-Z][0-9]+([A-Z])$", "\\1", unique_allels_table$mutation)
unique_allels_table$species <- sapply(strsplit(unique_allels_table$Organism, " "), function(x) paste(x[1:2], collapse = " "))

allels_table<-unique_allels_table[,c(24, 19,  2, 6, 7)]

write.csv(allels_table, file = "Quinolone_resistance_gene_table.csv")


