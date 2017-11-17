# Internship at the Cordelier Research Center - Paris

Scripts written and used during a part-time bioinformatics internship at the CRC (1 day/week during 5 months)

## twoCellstoKM

This script need a data matrix from the R package MCP-counter as input and is used to study couples of immune cell populations.

MCP counter works on transciprtomics data and provides a matrix describing 
the infiltration of tumoral tissue by immune cell populations.

Functions in it compute variations of K-means clustering with K = 4 
(similar to a double cutoff, beetween high and low infiltrated, for each cell).
it then draw a Kaplan Meier estimator of the survival function.


## mergeClust

Working on clusters instead of values allows to merge different kind of data (for example chips and RNAseq) and draw a new
Kaplan-Meier estimator with a bigger cohort
