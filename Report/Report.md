
# Case Studies in Bioinformatics, Report

## Module 1: Is the hourglass model for gene expression really supported by the data?


#### Antoine Gurtler, December 2018


### Introduction
During vertebrate ontogeny, different genes or gene sets are expressed. These genes drive the development through multiples stages. During the 19th, based on morphological observations, the first embryologists logically proposed a stepping-stone model where old embryogenic characteristics are conserved. This model leads to place a 'phylotypic' stage at the beginning of the development (funnel model). Since, it has been observed that a certain stage of development was morphologically common to all kind of vertebrates. But this stage is not at the beginning of the development. It is found in the middle of the development. It gave birth to the hourglass model. According to it, the mid-embryonic period is conserved among all vertebrates when embryos show, before and after this stage, divergences (1).  This morphological observation was first molecularly supported by Hox genes expression (2). This stage is supposed to be the most constraint and so, the most conserved. Some studies tried to highlight a relation between gene conservation and gene expression among developmental stages. Domazet-Lošo and Tautz (3) used phylostratigraphy and stage-specific gene expression data to compute an index representing the evolutionary age of the transcriptome for each stage. Phylostratigraphy gives to genes from a genome a phylogenetic rank. 
This report investigates the figure 1a of the Domazet-Lošo and Tautz publication (3). It explains the process used to product this figure, from the untreated data to the figure and discuss the latter. 


### Methods
We reproduced the figure 1a of the paper (3) on python 3 using the following packages: pandas, numpy, GEOparse and matplotlib.pyplot. We used the available expression data from the article and the dario_age_index which links probID, gene name and age of the gene together. 

First, we have to dowload and extract gene expression data and metadata from the GPL file. These meaningful metadata are the stage, time, sex and sample name. 


```python
import pandas as pd
import numpy as np
import GEOparse
import matplotlib.pyplot as plt
```


```python
############# Download the data
file_name = 'GSE24616'
gse = GEOparse.get_GEO(geo=file_name, destdir="./")
########### Read in age index data
age_index = pd.read_csv('danio_age_index.txt', sep='\t', header=None)
age_index.columns = ["GeneID","ProbeID","age"]
########### Set ProbeID as the index of dataframe
age_index.set_index('ProbeID',inplace=True)

```


```python
############### Extract GSE metadata
characteristics = {"stage":[],"time":[],"sex":[],"sample_name":[]}
for gsm_name, gsm in gse.gsms.items():
    characteristics["stage"].append(gsm.metadata['characteristics_ch1'][1].split(":")[1].rstrip())
    characteristics["time"].append(gsm.metadata['characteristics_ch1'][2].split(":")[1].rstrip())
    characteristics["sex"].append(gsm.metadata['characteristics_ch1'][3].split(":")[1].rstrip())
    characteristics["sample_name"].append(gsm_name)
char_df = pd.DataFrame(characteristics,index = characteristics["sample_name"])


# Exctract expression data and add ProbeID as the index of gene expression dataframe

data = gse.pivot_samples('VALUE') 
gpl = list(gse.gpls.values())[0]
data.set_index(gpl.table.SPOT_ID,inplace=True)
data.head()
```


Then we reduced the gene expression datafram to the probID found in the age_index index. This file also contains the gene associate with the prob and an integer linked to the gene age (stratigraphy index) (1 = old, 14 =young)


```python
matched_data =data.join(age_index,how='inner').groupby(level=0).last()
```

Some genes are reported by multiples probes, so we merged their expression together to have one line per gene. Authors of the original paper (3) did not do it (4). Then, we excluded male samples.

Then, we calculated the average gene expression for similar time points. There is 61 time points, separated in eleven developmental stages. We also associate a developmental stage to the corresponding time points. 


```python
# Average out multiple transcripts
unique_data = matched_data.groupby('GeneID').mean()

# Select only female and mixed sample.
char_df= char_df[char_df["sex"]!=' male']

########## find samples of the same time points
experiment_index=[]
time_stamps = char_df.time.unique()
for t in time_stamps:
    experiment_index.append(char_df[char_df['time']==t].index.tolist())

########### average the samples for similar time points
set_mean = {}
stages = []
timestamps=[]
for d in range(len(experiment_index)) :
    sample_list=experiment_index[d]
    set_mean[d]=unique_data[sample_list].mean(axis=1)
    stages.append(char_df[char_df.index.isin(sample_list)].stage[0])
    timestamps.append(char_df[char_df.index.isin(sample_list)].time[0])
mean_data = pd.DataFrame(set_mean)
mean_data.columns= stages
print(mean_data.shape)
mean_data.head()
```

The final table has 12892 rows and 61 columns, for 12892 genes and 61 time points. 

### Results. 

With this table containing the average gene expression per time point, we can calculate the transcriptome gene index (TAI) of each time point as decribed in the paper. 

$$TAI_s = \frac{\sum_{i=1}^n ps_i e_i}{\sum_{i=1}^n e_i}$$

* $ps_i$ is the phylostratum number of the gene i and $e_i$ its expression.






```python
########### Calculating TAI
age_indices =unique_data['age']
expression_data = mean_data.values
product = np.dot(expression_data.T,age_indices)
mean_expression = expression_data.T.sum(1)
TAI = np.divide(product,mean_expression)
```

We reproduced the figure 1a of the publication by plotting the 61 calculated TAIs. 


```python
#### define color map
unique_stages= list(set(stages))
color_list = plt.cm.tab10(np.linspace(0, 1, len(unique_stages))) 
color = {unique_stages[i]:color_list[i] for i in range(len(unique_stages))}
my_col= [color[i] for i in stages]

plt.figure(figsize=(20,10))
plt.scatter(range(len(TAI)),TAI,color=my_col,linestyle='-')
plt.xticks(range(len(TAI)),timestamps,rotation=45,size=12)
plt.ylabel("Transcriptome Age Index",size=30)
plt.xlabel("Developmental timing",size=30)
markers = [plt.Line2D([0,0],[0,0],color=c, marker='o', linestyle='') for c in color.values()]
plt.legend(markers,color.keys(),loc='upper center')
plt.show()
# write your code here 
```


![png](output_12_0.png)


*Figure 1: **Figure 1a reproduction** The x-axis represents the different time points, y-axis is the TAI. Color of the points corresponds to the developmental stage.*

According to this figure, we see that genes expressed during the mid-embryonic period have the lower TAI, so they are the oldest genes. 

However, if we look at the distribution of gene expression values, we see that most of the data have a really low expression. It means that only a few genes, ones with a strong expression, contribute to the pattern observed, which, surely by chance, correponds to a hourglass pattern. 


```python
plt.figure()
plt.hist(expression_data.reshape((expression_data.shape[0]*expression_data.shape[1],1)),bins=100)
plt.title('Distribution of the gene expression values')
plt.xlabel('Gene expression value')
plt.ylabel('Number of observations')
plt.show()
```


![png](output_15_0.png)


*Figure 2: **Distribution of the gene expression values**.*

If we apply a log-transformation to the data. They display a log-normal distribution which allows to compare data to each other. 


```python
plt.figure()
plt.hist(np.log(expression_data.reshape((expression_data.shape[0]*expression_data.shape[1],1))),bins=100)
plt.title('Distribution of the gene expression values with log-transformed data')
plt.xlabel('Gene expression value')
plt.ylabel('Number of observations')
plt.show()
```


![png](output_17_0.png)


*Figure 3: **Distribution of the log transformed gene expression**.*

Then, we can redo the plot with the normalized data.


```python
########### Calculating TAI


age_indices =unique_data['age']
expression_data = np.log(mean_data.values)
product = np.dot(expression_data.T,age_indices)
mean_expression = expression_data.T.sum(1)
TAI = np.divide(product,mean_expression)

#### define color map

unique_stages= list(set(stages))
color_list = plt.cm.tab10(np.linspace(0, 1, len(unique_stages))) 
color = {unique_stages[i]:color_list[i] for i in range(len(unique_stages))}
my_col= [color[i] for i in stages]

plt.figure(figsize=(20,10))
plt.scatter(range(len(TAI)),TAI,color=my_col,linestyle='-')
plt.xticks(range(len(TAI)),timestamps,rotation=45,size=12)
plt.ylabel("Transcriptome Age Index",size=30)
plt.xlabel("Developmental timing",size=30)
markers = [plt.Line2D([0,0],[0,0],color=c, marker='o', linestyle='') for c in color.values()]
plt.legend(markers,color.keys(),loc='upper center')
plt.savefig('TAI.png')
plt.show()
```


![png](output_19_0.png)


*Figure 4: **Transcriptome age index with log-transformed data**.* 

The normalized data changes the pattern of the graph. The lowest TAI are found during the first steps of the developpment. Trasncriptome age decreases over time. 

We can also test the sensibility to outliers by using another kind of normalization which is the z-score. This score indicates how many standard deviation a value is from the mean. 

Z-score formula is the following:

$$Z = \frac{X-\mu}{\sigma}$$

where X is the value of the element, $\mu$ the mean and $\sigma$ the standard deviation. 
source:https://stattrek.com/statistics/dictionary.aspx?definition=z_score

If we plot the data using this transformation, we observe an outlier on the fourth value of the segmentation. 



```python
########### Calculating TAI
age_indices =unique_data['age']
expression_data = (mean_data.apply(lambda x:(x-(sum(x)/len(x)))/np.std(x))).values
product = np.dot(expression_data.T,age_indices)
mean_expression = expression_data.T.sum(1)
TAI = np.divide(product,mean_expression)

#### define color map

unique_stages= list(set(stages))
color_list = plt.cm.tab10(np.linspace(0, 1, len(unique_stages))) 
color = {unique_stages[i]:color_list[i] for i in range(len(unique_stages))}
my_col= [color[i] for i in stages]

plt.figure(figsize=(20,10))
plt.scatter(range(len(TAI)),TAI,color=my_col,linestyle='-')
plt.xticks(range(len(TAI)),timestamps,rotation=45,size=12)
plt.ylabel("Transcriptome Age Index",size=30)
plt.xlabel("Developmental timing",size=30)
markers = [plt.Line2D([0,0],[0,0],color=c, marker='o', linestyle='') for c in color.values()]
plt.legend(markers,color.keys(),loc='upper center')
plt.show()


```


![png](output_21_0.png)


### Discussion:


Is the hourglass model for gene expression really supported by the data? We can answer this question now and the answer is no. This study shows how important is to normalize data. In fact, only outliers and strong expressed genes contribute to obtain a hourglass model. Data seems to support the funnel model, the one where the phylotypic stage is at the beginning of the developpement.  



### References:


1. 	Kalinka AT, Tomancak P. The evolution of early animal embryos: conservation or divergence? Trends Ecol Evol. 2012 Jul;27(7):385–93. 
2. 	Duboule D. Temporal colinearity and the phylotypic progression: a basis for the stability of a vertebrate Bauplan and the evolution of morphologies through heterochrony. Dev Camb Engl Suppl. 1994;135–42. 
3. 	Domazet-Lošo T, Tautz D. A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature. 2010 Dec;468(7325):815–8. 
4. 	Piasecka B, Lichocki P, Moretti S, Bergmann S, Robinson-Rechavi M. The Hourglass and the Early Conservation Models—Co-Existing Patterns of Developmental Constraints in Vertebrates. PLoS Genet [Internet]. 2013 Apr 25 [cited 2018 Dec 29];9(4). Available from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3636041/
 




```python

```
