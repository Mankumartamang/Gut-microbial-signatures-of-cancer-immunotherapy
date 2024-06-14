# Install MixOmics package.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")

#This R project involves the preprocessing and analysis of data from the four studies( 1. Baruch et al , 2. Lee et al, 3. Spencer et al, 4. Matson et al.; 
total num of samples= 328, num of shared ASV among four studies= 37).

#UPLOAD ABUNDANCE AND METADATA:
abundance1 <- read_excel("Abundancedata1.xlsx")# abundance1 is the integrated dataset containing ASV counts data from four studies.
metadata1 <- read_excel("metadata1.xlsx")# metadata1 contains metadata information of the four studies. Our variable of interest is "clinical response to immunotherapy".
class(abundance1)

#DATA PRE-PROCESSING:
#Since our microbiome data are compositional, there are three data preprocessing steps to be done on the raw count data:applying the offset,2. pre-filtering, 3. center log ratio(CLR) transformation.

#1.Applying the offset: Add an offset of 1 to the whole raw data matrix, as this helps handle zero counts and is recommended when using mixomics.
offsetdata <- abundance1 + 1 #when ran this code, it gave error report ("non-numeric argument to binary operator),the abundance1 data contains both character and numeric arguments.

str(abundance1)
#so, the above operation involves adding 1 to all count values and obviously abundance data frame contains samples on the first column(which is character vector),so can't do arithmetic operation(adding 1) on character vector
#so,below steps will only select numeric vectors for adding 1.


#extract the numeric count data from abundance1.The code below [, -1] will select all columns except the first.
countdata1 <- abundance1[, -1]

#add 1 to all count values.
countdata1_1added <- countdata1 +1

#convert the result back to a data frame and add the sample identifiers.
countdata1_1added <- as.data.frame(countdata1_1added)
abundance2 <- cbind(sample = abundance1[, 1], countdata1_1added)


#check the data frame correctly combined.
view(abundance2)

#2. Pre-filtering:
#removing all ASV for which sum of counts are below a threshold (0.01%). make sure this function requires samples to be on rows and counts on columns.

low.count.removal <- function(
  countdata1_1added, # OTU count df of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
) 
{
  keep.asv = which(colSums(countdata1_1added)*100/(sum(colSums(countdata1_1added))) > percent)
  data.filter = countdata1_1added[,keep.asv]
  return(list(data.filter = data.filter, keep.asv = keep.asv))
}

# call the function then apply on the offset data
result.filter <- low.count.removal(countdata1_1added, percent=0.01)
data.filter <- result.filter$data.filter

#add the sample identifiers to the data frame.
abundance2.filtered <- cbind(sample = abundance2[, 1], data.filter)

#An extra step that is recommended is to check how heterogeneous the library sizes are per sample. Here the sum of all counts per sample is calculated and represented in a barplot
#calculate feature counts. must select only feature count values from the data frame to generate the barplot.
feature_counts <- rowSums(abundance2.filtered[, -1])

#extract sample identifiers.
sample_identifiers <- abundance2.filtered$sample

#create barplot.
barplot(height = feature_counts, 
        names.arg = sample_identifiers, 
        las = 1,  # Rotate the labels to make them readable if they are long
        cex.names = 0.5,  # Adjust the label size
        col = "blue",  # Color of the bars
        main = "Feature count distribution by sample",
        xlab = "samples",
        ylab = "Feature Count")


#NOTE: Depending upon the context, a criteria to remove samples with extremely large library size may be applied, however, this project has not used that for this analysis. 

#3.centered log ratio (clr) transrformation.
#clr transformation may be used directly in some functions of mixomics (pca, plsda) or the data can be clr transformed first, then processed further.

data.clr <- logratio.transfo(as.matrix(abundance2.filtered[,-1]), 
                             logratio = 'CLR', offset = 0)

view(data.clr)

---------------------------
#DATA PROCESSING:

pca <- pca(data.clr) # undergo PCA on CLR transformed data

plotIndiv(pca,  # plot samples
          group = abundance2.filtered$sample, 
          title = 'Data1 PCA',
          legend = TRUE)

