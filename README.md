## Generalized Dissimilarity Modeling Workshop 

In the summer of 2023, the Earth Observation and Remote Sensing Lab at UC Merced invited identified researchers interested in spatial ecology to explore the statistical technique of generalized dissimilarity modeling through an intensive week-long workshop.   

Modeling the compositional turnover of species, known as beta diversity, is often a challenging statistical process, due to the complexities of nonlinearity, as well as the violation of standard independence assumptions common to the statistical approaches we are used to. To overcome some of these challenges, we will explore generalized dissimilarity modeling (Ferrier et al. 2007), which has been developed to specifically account for these issues in modeling. 

### **Purpose:**   

To gain familiarity with the statistical technique of generalized dissimilarity modeling (GDM) for analyzing and predicting spatial patterns of compositional turnover (i.e., beta/zeta diversity).  

### **Objectives:**

Develop an understanding of the statistical framework underpinning GDMs (i.e., matrix regression, I-splines, link functions, assumptions, and limitations). 
Apply GDMs to community composition datasets (Delta, vernal pools, and riparian) to analyze and describe spatial patterns of compositional turnover using remote sensing and other raster-based datasets. 
Explore extensions of GDMs for various purposes, such as Sparse Generalized Dissimilarity Modeling for high dimensional datasets (i.e., hyperspectral) and Multi Site Generalized Dissimilarity Modeling for zeta diversity. 

### R Packages Required: 

gdm - https://mfitzpatrick.al.umces.edu/gdm/ 

sgdm – https://rdrr.io/github/steppebird/sparsegdm/ 

zetadiv - https://www.rdocumentation.org/packages/zetadiv/versions/1.2.1 

### Summary of GDMs (from Ferrier 2007 and Mokany 2022): 

Generalized dissimilarity modelling (GDM) is a statistical technique for analyzing and predicting spatial patterns of turnover in community composition (beta diversity) across large regions. The approach is an extension of matrix regression, designed specifically to accommodate two types of nonlinearity commonly encountered in large-scaled ecological data sets:  

- the curvilinear relationship between increasing ecological distance and observed compositional dissimilarity between sites  
- the variation in the rate of compositional turnover at different positions along environmental gradients.

GDM can be further adapted to accommodate special types of biological and environmental data including, for example, information on phylogenetic relationships between species and information on barriers to dispersal between geographical locations. The approach can be applied to a wide range of assessment activities including visualization of spatial patterns in community composition, constrained environmental classification, distributional modelling of species or community types, survey gap analysis, conservation assessment, and climate-change impact assessment. The modeled results provide meaningful information on compositional turnover and their environmental drivers, allowing for hypothesis testing and mapping of compositional change. 

### Suggested Readings: 

Ferrier, Simon, Glenn Manion, Jane Elith, and Karen Richardson. "Using generalized dissimilarity modelling to analyse and predict patterns of beta diversity in regional biodiversity assessment." Diversity and Distributions 13, no. 3 (2007): 252-264. 

Mokany, Karel, Chris Ware, Skipton NC Woolley, Simon Ferrier, and Matthew C. Fitzpatrick. "A working guide to harnessing generalized dissimilarity modelling for biodiversity analysis and conservation assessment." Global Ecology and Biogeography 31, no. 4 (2022): 802-821. 

Leitao, Pedro J., Marcel Schwieder, Stefan Suess, Inˆes Catry, Edward J. Milton, Francisco Moreira, Patrick E. Osborne, Manuel J. Pinto, Sebastian Van der Linden, and Patrick Hostert. "Mapping beta diversity from space: Sparse Generalised Dissimilarity Modelling (SGDM) for analysing high‐dimensional data." Methods in Ecology and Evolution 6, no. 7 (2015): 764-771. 

Leitão, Pedro J., Marcel Schwieder, and Cornelius Senf. "sgdm: An R package for performing Sparse Generalized Dissimilarity Modelling with tools for gdm." ISPRS International Journal of Geo-Information 6, no. 1 (2017): 23. 

Latombe, Guillaume, Cang Hui, and Melodie A. McGeoch. "Multi‐site generalised dissimilarity modelling: using zeta diversity to differentiate drivers of turnover in rare and widespread species." Methods in Ecology and Evolution 8, no. 4 (2017): 431-442. 
