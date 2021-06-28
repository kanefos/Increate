# INCREATE: Single-cell RNA Sequencing Workflow 

Workflow for the analysis of single-cell RNA sequencing data.

## Contents
1. Preprocessing and quality control.
2. Clustering and phenotyping.
3. Differential gene expression.
4. Compositional analysis.
5. Cell-cell communication analysis
6. single-cell TCRs

Other useful resources:
- https://github.com/crazyhottommy/getting-started-with-genomics-tools-and-resources
- https://github.com/mdozmorov/scRNA-seq_notes

## Setup
- Running on conda _'IncreateV2'_ environment. 
- Operate the Increate pipeline itself (i.e. scripts to process data and perform analyses) out of the appropriate directory
- However, data is physically _stored_ inside dataset-specific directories which are specified when a script (in R or Python) is initiated 
    - Not solely datasets, also subclusters, integrated data, etc
- Utilise standardise file nomenclature and formats, i.e. h5ad, functions to move from R to Python


## 1. Preprocessing and Quality Control

Processing the cell x gene matrix. Removing dead/dying cells, doublet events, and overt contaminants (i.e. RBCs).

Strategy for preprocessing: 
1. Firstly, run conventional QC (tripartite, acquire doublet scores, annotate TCR info etc)
2. Run full clustering pipeline on this object, with high-resolution clustering (leiden resolution = 3 & 4)
3. On this object, isolate stressed/dying cell and doublet clusters
4. Then, rerun to generate the first dataset clustered object. From here, can isolate lineage subsets etc (which may have additional QC steps performed on them)

### Sequence alignment
Standard 10X STAR alignment
- Issues are reported with CellRanger:
  - https://github.com/10XGenomics/baofastq/issues/23
  - Solution? See here: https://github.com/Zha0rong/Bam_Parser_10X

Also interested in in _kallistobus_: 
- https://github.com/Acribbs/scflow
- https://github.com/pachterlab/kb_python 

### 1.1 DataImport.py
- Processing from uploaded data or raw feature-barcode-matrice files (i.e. own experiments).
- Externally, process data to either a standardise .h5ad file (raw data), or an assembly of X, obs and var
- Perform metadata standardisation as seen below, ensure for integration when similarities occur accross datasets shared labels are appropriately named (i.e. 'Lung' versus 'lung')

#### Cells/obs
Essential
- Study name or data source: 'study'
- Sample name or batch: 'sample'
Additional:
-  - Tissue status i.e. disease/disease stage: 'tissue_status'
Other:
- Usually study-specific, no standardisation less important, i.e. cell sort: 'sort' (identical to 'sample')

#### Genes/var
Essential:
- Ensemble ID: 'ENS'
- Gene name: 'gene' (also for 'var_names')
Additional metadata:
- Will include chromosome (for X/Y sex removal) 
- This also annotates genes localised to the mitochondrial genome. Use *biomart* first.

#### Ambient expression
- I experimented with SoupX - while it seems impressive, rough experiments with the b16sc data were unsatisfying. I think if perfected and really paid attention too, would be worthwhile, but sufficient analysis possible without it.

### 1.2 PreProcessing
Misc  notes
- QC covariates should be considered _jointly_ when univariate thresholding decisions are made
- Go easy on number of genes - T-cells can have a low number of genes per cell
- Need to really closely massage this data, on a sample-wise, then experiment-wise, then dataset-wise way.
- Biggest QC for ample prep readout is dead cells.

#### TCR Data
- If 10X vdj acquired, can be applied early (to the `umi.h5ad` object) in a reduced form (`has_ir`, `chain_pairing`, `clone_size`) etc, as can be informative about cellular composition.

#### Gene filtering
- Minimum 3 cells expressed in 

#### Conventional QC
Count depth
 - >1300-1500
Number of genes
- >300-500
Percentage mitochondrial counts 
- Human <10%, mouse <5%. 
Percentage ribosomal counts
- Use sparingly. Had success in removing RBC from haematological samples, but there is too much variation in baseline cellular % Ribo expression to rely on this metric for downstream work. 

#### Doublets
- _doubletdetection_ tool, default settings
- _scrublet_ tool, default settings
- Do on raw matrix, sample-wise (note, post-QC)
- By using 2 orthogonal approaches, better interrogate likelihood of double encapsulation events
- Useful tools, but marker genes are usually key here also.

#### Misc Steps
MALAT1
- Remove MALAT1 - general QC? 
- https://kb.10xgenomics.com/hc/en-us/articles/360004729092-Why-do-I-see-high-levels-of-Malat1-in-my-gene-expression-data-
Sex Chromosomes
- Mark Y-linked (male exclusive)
- Can view women by their expression of XIST (expressed solely on inactive X chromosome)
Ig Genes
- Mark via gene string start IGH/IGL/IGK
- Dotplot suggests not as much an issue but worth doing

#### Contaminant Cells
- Coarse cell type identification
    - Very low-level: i.e. red blood cells, tumour cells

#### Assessing outputs
See output_df, preprocessing_output_df.csv (note: saved to graphs dir)
- Self-explanatory but info can be further acquired in-script
]

### 1.3 PreNorm.py
#### Sample-wise Coarse Clustering
- Leiden, scanpy standard workflow
- Export dataset and coarse leiden clustering
    - Scran documentation: "Only a rough clustering is required here, as computeSumFactors is robust to a moderate level of DE within each cluster."
- Export a .h5ad file of the dataset with sample-wise clustering md? i.e. all in one anndata object, but can pull out individuals 
    - For visualising issues around scran negative size factors 

#### Subset reanalysis
- This is the point in the workflow where re-normalisation, correction, and clustering etc of specific subsets would occur, i.e. T-cells, Myeloid cells etc
- There is a specific section dedicated to creating data like this

### 1.4 ScranNorm.R
Normalisation
- scran size factor normlisation, coarse Leiden clustering as input clusters
- May need to cut any samples/coarse clusters giving rise to negative size factors (in PreNorm.py)
- Apply to dataset, log2-transform
- 
### 1.4.1 sctransform.R
Normalisation
- Perform the suposedly count-depth robust sctransformation and HVG selection tool within Seurat wrapper
- After using on Bailur 2019 whole-dataset, did seem to identify some residual RBCs but also (on a macro scale) seem to cause some CD4/CD8 mixing not evident using scran
    - So: notably different to scran at this level
- See clustering/misc dir - may partially mask T-cell subsets at whole-dataset level
- Warning! `sctransform_matrix.npz` is a vast file - about 7GB of my 500GB memory Mac!

### 1.5 PostNorm.py
#### Apply scran normalisation
- Numpy vectorized computing critical here
#### Total-count normalisation
- i.e. scanpy default. Also valid. Had success with low-depth, non-batch corrected b16sc data.
#### Gene Scaling 
- Done later. Only needs performing on HVGs used to calculate PCA.

### 1.x Integration.py
- Integration is to be done, when possible, purely for visualisation purposes (at least when using Harmony, a batch-corrected expression matrix can always be used later, i.e. DEGs)
    - i.e. project, DEGs (& pathways), pseudotime, all done on individual datasets but then projects onto same low dim graph 
- Perhaps as a test - remove ½ of a sample/cell type/cluster and then run entire pipeline again - see if we get the same results?
#### Integration workflow
- Utilises Increate_V1 legacy scripts: Preintegration.py, DataIntegrationV4_Oetjen.R, Postintegration.py
Workflow
0. Uses QC'd, sample-wise coarse Leiden scran-normalised
1. Study-wise scran:modleGeneVar() (may require block= on sample), union of outputs (recommended threshold 0.1/0.5)
2. batchelor:multiBatchPCA, batch=*sample* (not *study*, but could be tried also)
3. Input harmony, default settings with iteration of theta values on study/sample (final employed no cosine normalisation)
4. harmony-adjusted PCs feed into scanpy clustering pipeline
Adam email:
- "One of my postdocs has a pipeline that she finished where she tests the SCT and harmony. This is for a separate endometrial tissue project and im waiting for results from her. A PhD student in the lab will use the pipeline on myeloma patient samples and I will have the data soon, hopefully."
### 1.x Integration.R
- Notes for above



## 2. Clustering-Phenotyping

### 2.1 Clustering
Takes normalised data with a majority of contaminant populations removed.

HVG Selection
- Scanpy defaults
Batch Correction with Harmony
- Harmony with default settings on `samples` 
- Annotate PCA to matrix
- Tried without, usually looks messy (both on UMAP and with clustering)
Prerequisites
- HVG selection: have used dataset-wise scanpy defaults, but still slightly skeptical - idea of combining HVHs across samples still appealing
- Gene scaling: still no consensus so for now can run both
    - May need to do anyway just for my data to interface with other resources (which were scaled)
    - Note: could I use non-scaled clustering / visualisations with scaled raw data? sounds a bit iffy


#### Clustering pipeline

##### How many PCs to use? 
- Reduce noise noise
- Choose the number of PCs such that we can explain 90% of the initial data dispersion

##### Neighbors

n_neighbors
- controls balance between local and global structure, by constraining size of local neighborhood
- the lower the value, the more local the concentration
- description is pretty accurate, when tested a range of values on 
    - medium (bailur all) nonscaled data (scaled scanpy HVGs) seemed to change a lot below 20-30 after which stablised
    - small (bialur T cells) nonscaled data (scaled scanpy HVGs) roughly same as above in clustering but umap fairly constant throughout

metric
- A known metric’s name or a callable that returns a distance.
- Ehsan recommend cosine
    
- Switching between these two does change a lot - but again it's slightly hard to interpret
    - On an intermediate size dataset, seemes to make a little messier? i.e. over-preserving local structure?
- May be better to just outright go with cosine - ehsan recommended and can tune downstream anyway


##### UMAP

###### min_dist and spread

min_dist
- default 0.1, recommend initially trying between 0.1 and 1)
- above 3, broke (on bailur T cells)
spread (default 1)
- below 0.3 and above 5 didn't work

spread doesn't move clusters much, but min_dist does


# - a, b (harder to explain, start at both =1 and mess around a bit)



##### Leiden and Louvain algos
- These require network graphs - hence extra info on the anndata object
- n_neighbours can be altered! critical
- Think I should properly just keep this as anndata-centric
Within scanpy:
- The procedure of clustering on a Graph can be generalized as 3 main steps:
    1) Build a kNN graph from the data
    2) Prune spurious connections from kNN graph (optional step). This is a SNN graph.
    3) Find groups of cells that maximizes the connections within the group compared other groups.
- If you recall from the integration, we already constructed a knn graph before running UMAP. Hence we do not need to do it again, and can run the community detection right away.
    - Basically: done when we do `sc.pp.neighbors()`
- The modularity optimization algoritm in Scanpy are Leiden and Louvain.

- Cluster granularity - most important on i.e. just T-cells
    - Graph here might also be interesting?

- Iterate resolution! Then pull out some cluster ML metadata

- Leiden/Louvain and clustree?
    - Would be a nice way to support idea of transitional clusters that breakaway early over resolution
    - Does clustree support resolution data input thusly? yes!
    - Also is clustree rigerous against changing metric (k, resolutionm, etc) that doesn't result in a change of cluster numbers / membership?


#### Returns
- A clustered dataset, be it a full dataset or a sub-cluster (i.e. on cell type)
- Dataset composed of solely cells of interest
- Several Leiden resolutions for different n_clusters

- After clustering, may ned to remove contaminant populations, and re-run.
- Should also run pipeline from raw data again on: subclusters of interest, i.e. T-cells, TIMS, etc
    - Remember: sub-clustering can be done as many times as on likes (at least until number of cells limitng)


### 2.x Clustering_extra

This is an additional script / resource of materials for the downstream deeper analysis of specific clusters, after having done clustering/phenotyping - also some visual material i.e. clustree

- Various other clustering approaches, i.e. Leiden, Louvain, Heirarchal etc
- Cluster metadata, i.e. QC and numerical info, i.e. depth, number of genes etc
- Clustering metadata:
    - Cluster membership as a member of resolution (perhaps a graph?)
    - Entropy, see ROGUE

#### Heirarchal Clustering
- Only clustering method to include the hierarchical information about subclusters, may be useful for subcluster/subtype heterogeneity
- While Scanpy implements some decent dendrograms, not true clustering and instead based on final label itself

##### Agglomerative clustering
- AKA bottom-up approach or hierarchical agglomerative clustering (HAC), is more informative than the unstructured set of clusters returned by flat clustering. 
- This clustering algorithm does not require us to prespecify the number of clusters. 
- Bottom-up algorithms treat each data as a singleton cluster at the outset and then successively agglomerates pairs of clusters until all clusters have been merged into a single cluster that contains all data.
- The results can be visualized as a dendogram to help recognize the moment the algorithm should be stopped to get optimal results.

##### Divisive clustering
- AKA top-down approach. 
- This algorithm also does not require to prespecify the number of clusters. 
- Top-down clustering requires a method for splitting a cluster that contains the whole data and proceeds by splitting clusters recursively until individual data have been splitted into singleton cluster.
- Versus agglomerative cluster, divisive algorithm is  more accurate.
    - Agglomerative clustering makes decisions by considering the local patterns or neighbor points without initially taking into account the global distribution of data. 
    - These early decisions cannot be undone. whereas divisive clustering takes into consideration the global distribution of data when making top-level partitioning decisions.
    -  _TooManyCells_ is a divisive hierarchical clustering workflow to catch and visualize the cell clades. 
        - It partitions the cells into two groups with the spectral clustering algorithm in an iterative manner. 
        - However, TooManyCells does not give the full hierarchy and makes it hard to set as specific cluster number.

#### Accuracy Metrics
- Unlike classification, an accuracy metric cannot depend on the labels but only on the goodness of split.
- There are _internal_ and _external_ goodness metrics: external metrics use the information about the known true split while internal metrics do not use any external information and assess the goodness of clusters based only on the initial data. 
- The optimal number of clusters is usually defined with respect to some internal metrics.

##### Label-independent:
Silhouette
- Lets us estimate the quality of the clustering using only the initial, unlabeled sample and the clustering result
- The silhouette distance shows to which extent the distance between the objects of the same class differ from the mean distance between the objects from different clusters
- Ascribes a number to 'cluster tightness’
- Closer to -1 suggests incorrect clustering, while closer to +1 shows that each cluster is very dense
`metrics.silhouette_score(scaled_feature_data, cluster_labels)`
Calinski Harabaz Index
- The Calinski Harabaz Index is the ratio of the variance of a datapoint compared to points in other clusters, against the variance compared to points_within_its cluster. 
- Since we want this first part to be high, and the second part to be low, a high CH index is desirable. Unlike other metrics we have seen, this score is not bounded.
`metrics.calinski_harabasz_score(scaled_feature_data, cluster_labels)`


#### Gueguen-2021 clustering metadata
- clustree cluster membership as a function of cluster resolution
    - A split that occurs early between 2 clusters indicates the importance and robustness of the differences between the two
    - 'The early separation (low resolution of 0.2) between the resident and circulating part of CD8+ indicates that these 2 subparts of the datasets are the strongest difference within the data'
- Also termed *cluster granularity* 
- Note: this actually required fairly "poor" clustering, i.e. didn't data peak a single good silhouette score for their intermediate pop

#### ROGUE
Cluster entropy with ROGUE (Liu et al 2020)
- ROGUE proposed the concept of cluster purity and introduced a conceptually novel statistic, named ROGUE, to examine whether a given cluster is a pure cell population.
- https://github.com/PaulingLiu/ROGUE
- Used to infer relative homogeneity/heterogeneity within specific clusters
- Similar to other tools like heircharcal clustering / membership as a function of resolution

Interpretation of rogue value
-  D cell population with no significant ds (maximal rogue reduction, it's a formula) for all genes will receive a ROGUE value of 1, indicating it is a completely pure subtype or state. 
-  In contrast, a population with maximum summarization of significant ds will yield a purity score of ~0
Zheng 2020:
- Using a selection of FACS-sorted single-cells (2g), CD4/CD8 naïve T cells and CD4 memory T cells (known to be highly homogeneous populations) had high detected  ROGUE values (0.94–1), as expected. 
- In contrast, both CD14 monocytes and CD34+ cells are mixtures of diverse subtypes and received relatively low ROGUE values (~0.8; Fig. 2h), thus confirming their heterogeneity.
- Tumor-infiltrating DC2 cells have been proven to be highly heterogeneous populations, in Fig 2J shown to have a reduced ROGUE value
- ROGUE-guided analysis not only discovered novel cell subtypes, but also enabled the detection of the true signals in specific pure subpopulations.
- Identified poorly-clustered cells? i.e. ROGUE values broadly increase

ROGUE might be interesting to play around with, esp with lots of labels of different cells types & samples & disease etc?


#### clustree

Deciding what resolution to use can be a difficult question when approaching a clustering analysis. One way to approach this problem is to look at how samples move as the number of clusters increases. This package allows you to produce clustering trees, a visualisation for interrogating clusterings as resolution increases.

To build a clustering tree we need to look at how cells move as the clustering resolution is increased. Each cluster forms a node in the tree and edges are constructed by considering the cells in a cluster at a lower resolution (say) that end up in a cluster at the next highest resolution (say). By connecting clusters in this way we can see how clusters are related to each other, which are clearly distinct and which are unstable. Extra information about the cells in each node can also be overlaid in order to help make the decision about which resolution to use.



### 2.1-2: Phenotyping
- starting to suspect seperating into two not the best idea
- this section is not as linear as the QC / clustering

However: broad approach (current ideas:
1. Projection tools (Azimuth, ProjecTILs etc)
2. DEGs with scran
3. fgsea for dedicated signatures / pathways / contrasts



### 2.2 Coarse Phenotyping

#### Phenotyping ideas Ideas

- Do NOT want to treat single-cell RNA seq like flow data, optimise entire spectrum
- Gene signatures: surely the TIGIT paper bulkRNAseq could be used as a general infiltrating score?
- Other bulk scores? Or other tools? 
- Think bio: i.e.
    - How to seperate NK-T from CD8s?
    - What are the different NK and B-cell subsets? Seen Cheng-2021-style review!

Residual contaminant population removal


#### Cycling cluster annotation
- Annotate with Swedish Scanpy tutorial methodology including regev lab cell cycle genes.
    - Importantly, this IS done post-normalisation/scaling etc
    - However, does difer a lot between different approaches, i.e. scaling versus non-scaling
- At the basic level, there is some connection between count depth and these cell scores
    - Do not be distracted by the total genes! depth is of interest here
- Does not need to be wildly restrictive, can re-slice additional cycling clusters detected later out

- Slice out bona fide cycling clusters for downstream DEG testing
    - Can always re-import clustered object and rename for whole-dataset cluster ID umap

#### Seurat Reference Projection

##### Pediction scores: 
- Cell prediction scores range from 0 to 1 and reflect the confidence associated with each annotation. 
- Cells with high-confidence annotations (for example, prediction scores > 0.75) reflect predictions that are supported by mulitple consistent anchors.

##### Mapping scores: 
- This value from 0 to 1 reflects confidence that this cell is well represented by the reference. T
- The mapping score reflects how well the unique structure of a cell’s local neighborhood is preserved during reference mapping.

#### Azimuth Reference Projection

Leverages a 'reference-based mapping' pipeline that inputs a counts matrix of gene expression in single cells, and performs normalization, visualization, cell annotation, and differential expression (biomarker discovery)

Impessive companion tool: useful but can't be relied on.
- Marker fishing: worked welll for HSPCs, not so great for ILCs

Should have faith in my clustering: recall the stain-glass window analogy.
- An individual 10X is really a flawed insight into the transcriptome, we need the centroid transcriptome of whole-clusters for sufficient "depth"
- See the Azimuth predictions for i.e. the B-cell clusters: get a couple of CD4 TCM and Tregs as predicted (which I doubt)
- Also: if I do a dotplot of markers etc and see the cells labelled as NK cells having recued CD3D etc, does that say much?
    - What's to say they simply didn't get marked as T cells due to this dropout of CD3D / other markers, but during whole clustering their total transcriptome did group them? have faith in my approaches!
- HOWEVER! that in mind, might be interesting to try some DEGs or TI approaches without these cells - see if things look cleaner/better

How useful for T-cells? / other subsets
- perhaps use downloaded Gueguen and map onto?
    - tried Gueguen, didn't work but pretty sure I did it wrong :/
- both of these ideas remain succeptible to having sufficient cell numbers
- with contribution of my own markers / tools (i.e. gsea) etc

Functionally may not need to totally abide by, i.e. ASDC versus pDC - much point differentiating here?

##### Another reference atlas tool:
- https://abseqapp.shiny.embl.de/
- From Triana 2021


#### Coarse: DEGS 
- Can use basic scanpy functionality, can simply describe as "scanpy defaults"
- Should be done early, leaf with statistically-supported DEGs
- Save coarse DEGs to object with correct label for leiden resolution
    - With Bailur, had good results with Leiden resolution 0.2 and 0.6
- I am retaining mito/ribo ect at this stage - can make is slightly difficult to identify some subsets i.e. active T-cells, but overall not much of an issue until later 
    - Note: these were remove for clustering / DEGs (I think?)

#### Coarse: Phenotyping
- Have a relatevely sound array of marker genes
- Well-informed by Azimuth


#### Coarse: Labelling and Subsetting
- Save objects out, move from here

#### Remove doublets from pops?
- i.e. if ONE T-cell has x2 many counts, likely doublet



### 2.3 Finer Phenotyping

#### Notes
- For now, will work with cells identified via scaled UMIs - but will re-normalise this data in several ways
- If one way seems to perform best on subset of data like this that isn't scaled UMIs, may return to this step and re-implement this approach for subset selection solely for completeness
- For plotting (via scanpy), utilise nonscaled normalised data ('nonscaled' layer) which is then z-normalised gene-wise within plot (it's functionally identical and aids plotting)



#### Finer DEGs (as opposed to coarse DEGs)
- This is what I envisage being used to find the differences between, say, two clusters of interest in the name of making violin plots etc
    - So for finer phenotyping, i.e. different T-cell effector subsets
- This is also where I'm keeping info about stats/DEGs in single-cell analysis
- See above: will need to regress out / remove non-testable genes i.e. ribo / mito etc (I assume, may be less of an issue on lineage subsets)
- make an anndata.layer for stats testing which has all genes mentioned above removed?

- Slice on gene of interest (!= ribo, mt, etc)
- How to control for cluster size? What is the actual version of my mental image of downsampling?
- How make sample = batch?
- Need to control for count depth?
- New idea: couple leiden resolution and stats testing - so do testing on macro clusters of i.e. T cells, B cells, etc, then once established do testing within clusters for T-cell subsets etc
    - This suggests should do iterative increased resolution reclustering of macro clusters versus set resolution high to begin with?
    - Does raise a question about cells that move upon re-cluster, i.e. on the periphary
    - Also, how does scanpy / other methods store stats info? a dict would be nice...
- Could very easily make a DEG testing layer

- Other contrasts we could do: i.e. normalised expression of TCF7 in T-cells in hd v mgus v mm (controlling for size), etc
- Think on other options, conditon- and cluster-wise
- Be more downstream-focused, i.e. volcanoes plots between clusters etc
- What about observed / expected stats test (if expected random?)
    - Pearson's chi-squared test is used to determine whether there is a statistically significant difference between the expected frequencies and the observed frequencies in one or more categories of a contingency table.
    - http://archive.bio.ed.ac.uk/jdeacon/statistics/tress9.html

Tommy Tang: DEGs for Marker Genes
- Use original data (versus corrected/integrated), blocking on batch (scran does this)

Statistical power
- somewhat related: https://satijalab.org/howmanycells/
- See large mdozmorov notes .md file: several other power tools
- Adhoc testing?

DEGs on batch-corrected counts
- Will get email about batch-corr from Adam
- Should always be on normalised data though
- Perhaps do step-wise on n clusters (i.e. from early time point - establish broad lineager marker genes i.e. CD3D which may be less significant when testing between multiple sub-clusters which include similar cell types)

#### Marker Gene Discovery with Scran

Bioconductor OSCA: http://bioconductor.org/books/release/OSCA/marker-detection.html#pairwise-tests-between-clusters

General:
- Do cluster pairwise contrasts, then combine into a single rank versus one cluster versus all other as this is sensitive to population composition and cell type abundances
- Pairwises comparisons also have more info: FCs
- Welch t-test: easy choice as it is quickly computed and has good statistical properties for large numbers of cells
- Upregulated genes are prefered for marker discovery (not to say totally invalid in sc stats-testing)

Scran approach: `findMarkers()`
- Scran's philosophy of findMarkers() is to identify a combination of marker genes that together uniquely define one cluster against the rest, indicated by the Top column
- Collects top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster
- `summary.logFC` is defined as the log-fold change from the comparison with the lowest p-value

What works:
- Welch's t-test: 
    - `pval=some` for a comprimise of strigency and generousity 
    - Only investigate upregulated genes.
    - `findMarkers(sce, groups=sce$leiden, block=sce$sample, pval.type="some", direction="up")`
- Binomial test
    - For lowly-expressed genes, also reduce's effect of mt/ribo
    - Also upregulated genes only
    - `findMarkers(sce, groups=sce$leiden, block=sce$sample, test="binom", direction="up")`
- All corrected for batch
- in terms of volcanoes, both generate fairly readable plots, although for both p-value does cap

Stats notes:
- Custom DE Methods such as edgeR and DESeq2 not routinely used for marker detection
    - Many of these methods rely on empirical Bayes shrinkage to share information across genes in the presence of limited replication which is unnecessary when there are large numbers of “replicate” cells in each group, and does nothing to solve the fundamental n=1 problem in these comparisons. 
    - These methods also make stronger assumptions about the data (e.g., equal variances for linear models, the distribution of variances during empirical Bayes) that are more likely to be violated in noisy scRNA-seq contexts. 
    - From a practical perspective, they require more work to set up and take more time to run.
- An alternative to `block=` is to define a design matrix containing the batch of origin as the sole factor, which is functional, can increase power, and can handle situations where multiple batches contain unique clusters - however use of linear models makes some strong assumptions (if batch effect is not consistent across clusters, the variance will be inflated and the log-fold change estimates will be distorted, variances are also assumed to be equal across groups, which is not true in general). 
    - Thus, we generally recommend the use of `block=` where possible
- data snooping/data dredging/fishing

Question:
- if I have 2 clusters and want to find DEGs between them, how do I do this with scran?
    - i.e. what's my "normal" stats testing approach
    - `pairwiseBinom` or  `pairwiseTTests` using `restrict=` to specific groups of interest
    - this works for factors passed to `groups=`, so to say do T cells in MGUS/MM would need to pass only T-cells in the sce object then do testing on `groups=disease` to properly model (and `exclude` HD, in this example)
    - so can be done within same testing framework
- Also a pesudobulk approach: see Liu 2021 systems immunology C-19 paper (requires a cutoff of cell numbers of total UMIs per sample-cluster pseudobulk contrast - can do in Limma voom)




#### DEGs: working with mt and ribo genes
- When testing, included mito/ribo/HSP (MRH) genes
- This does yield an abundance of such genes, particuarly for subsets like naive CD4 T cells
- Excluding them yielded similar results to simply sorting these genes out of the results file - to properly remove them for testing I think would need to fully detract these genes from the transcriptome prior to clustering - a step I am not willing to make
- Therefore, will simply work around at the testing level (just need to scroll through a few MRH markers)
- I did try FindAllMarkers() on Gueguen dataset: no ribo genes as far I can see, some HSPs but it was biology discussed/recognised in original paper
- Not to say I won't isolate these genes / gene groups later 

Working within sctransform ideas framework: housekeeping genes (i.e. RPL11):
`# - count depth for ribo? i.e. high num = higher counts (by chance)
sc.pl.umap(T_dataZ, color=['leiden', 'total_counts'], legend_loc='on data', legend_fontoutline=1)
sc.pl.scatter(T_dataZ, x='total_counts', y='RPL11') # RPL11 = high score as DEF
sc.pl.scatter(T_dataZ, x='total_counts', y='RPL11', color='leiden') #hard to interpet but useful
# - hmmm maybe`

#### DEGs into FGSEA
- Had best results with `findMarkers(sce, groups=sce$leiden, pval.type="some", direction="up")`
- Tried some FC cutoffs, >1.5 looked best (want around 1000 genes in? it's signal:noise ratio tradeoff)
- Incorporating downregulated genes (negative FC) confused the signal a bit - so perhaps only look for enrichment (versus depletion) score

!!!!
big realisation about Guo signatures
- the CD4 DEGs are between the two sets of T cells, not collectively among CD4s
- so (if we had more Tregs) could be done on Tregs alone, but are NOT part of the CD4 testing pannel
- will neaten stuff up slightly
- suggests removing Tregs akin to removing MAITs? idk
- this makes the final graph for CD4 fgsea pretty okay (suggested naive-like, but perhaps more too it? see GATA3 dotplot - hence further investigation)
!!!!





#### Should re-do it all, just ensure works
- while here, ruminate on if I am totally happy with clustering
    - i.e. didn't get HVGs on scaled data
    - also no MT/ribo exclusion for clustering?
    - should probably run through all again!
- it's nice to think this way, especially when ensuring I know what I did, but can be sat here optimising forever - if confident, run with

(lots of above can be cut)


#### True use of statistics
what needs testing?
- i.e. cluster/ID-wise, could do HD v MM v MGUS etc (determine cellular signatures specific to diseased BM single-cells), then old-school GSEA?
    - sample size/replicates, power, etc
Validation
    - this is the name of the game: should try to see if any changes in cellular abundance are reflected in healthy BM or MM datasets

Starting to realise... the power (pun intended) of gene signatures isn't in having eleborate sources of signatures that look flahsy, but instead to have suitable statistics across conditions of interest - it's generating the DEGs / stats that really matters!!
- As a rule:
- Fundemental:
    - GO BP
    - GSEA hallmarks collection: http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=H
    - KEGG
    - Reactome
- Immune:
    - Blood txnomic modules: https://github.com/shuzhao-li/BTM
    - Extensive but http://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=IMMUNESIGDB

"Level of testing": at T-cell level, at PBMC-esc level - think about where a contrast of DEGs are needed (esp in reference to the signature/dataset being contrasted with - think of Guo story!)
- BM immune Atlases might be nice here



#### More critical stats reading!
Effect sizes
- https://www.leeds.ac.uk/educol/documents/00002182.htm
- https://www.simplypsychology.org/effect-size.html

#### CD4 heterogenity: ROGUE-clustree
How do the T cell clusters look with ROGUE?
- is there any implied heterogeneity?
- might be an idea to do this soon-ish before doing loads of DEG interpretation on CD4s (ie SHOULD I be expecting any sort of interesting signal? Or are they dimply all naive CD4s?)

overly mean (or mean plus/minus sd?) shouette or ROGUU onto clustree - better visualisse cluster stability  across resolution
- all this inter-R file movement may need some thinking about

*IF* it works with Leiden resolution instead of kmeans, use clustree to probe CD4 heterogeneity? As a support of rogue

May acc be simpler to just use kmeans? But then again Leiden resolutions (perhaps 0.1 or 0.05 increments versus 0.2?) should work fine

There may NOT be any underlying structure to the CD4s
- And even if there is who cares? Unless there’s some enrichment in MM v MGUS etc != matter

Lit requirement: NEED DER LEUN 2020 FOR CD4s - perhaps NOT TILs as MGUS?? More conventional CD4 subsets

aside:
Anything interesting about this ribo gene expression particularly?
https://www.sciencedirect.com/science/article/pii/S1074761303002681


#### stats testing: rRNA fraction
- Ribosomal RNA read fraction, Ribosomal protein read fraction
- Possible that degradation of RNA leads to more templating of rRNA-fragments.
- Proportion ribosomal proteins may be an artifact from handling of samples for cells of the same celltype
- But EEIFs etc suggest biological? not sure about this idea

#### what about simply filtering from stats test output?
- becomes more readable, but underlying issues remain


#### scHaystack (in Advanced mode)
- With cassette of DEGs? i.e. pre-cluster DEG for downstream cluster ID (avoid issue of subclusters of i.e. T cells not having CD3D a DEG all other clusters etc)
- Outputs (example graph)
    - Density: colored by the relative density of cells expressing each gene
	- Level: 
- scHaystack can take high-dim coordinates, which could be batch-corr/integrated PCA output by Harmony
- Could this apporach be used to assess intra-cluster heterogenity (i.e. within 'T cells') from a macrostandpoint - be fun to show more heterogeneity within T cells than other cell types 
    - Run on just T-cells, see how it does!
- What about using DEGs output from scHaystack as HVGs?
	- If we then cluster these DEGs, could even form an input for some modelling? i.e. block on different DEG-clusters? Idk
- The clustered of DEGs is a strange idea - presumebly it would align with regular clustering? Therefore you overlay the two to reveal DEGs in support of your clustering, could pull out biologically genes (or even gene groups if could be bothered) from the DEG list/clusters?
- Re: low sample depth - scHaystack could be very well-suited to reduced depth 10X data, as it would make data binary

Idea:
- Look at very final dataframe from advanced mode - has raw statistical info (alongside clusters, ignore for this precise point) so could simply pull out genes of interest and use this info?
- Not cell-cluster driven, so how fit? idk - but could wrangle to make any intelligent/novel use of stats

Also: course of novel genes? i.e. do discovery / volcanos/heatmaps etc here?
- The top haystack DEGs do make nice marker genes, very good expression on dots/umap
    - However - note sometimes NOT uniform across cluster as indicated by sc haystack
- be nice to see how this works with more clusters too! (may partially eleviate very above point)
- *discovery*: return to bio!

Can have many clusters! To make leiden resolution? Or iterate k here too?

If do combine this expression data to cluster, may need to think hoew matrix-creating-leiden clusters interface with these - i.e. same expression format



#### gseapy

#####  Hallmarks
- Pathway and TCGA data analysis. 
- To characterize and detect the pathway signals in specific fibroblast subtypes, we performed pathway analyses using hallmark pathways from the molecular signature database54 with GSVA55. 
- Copy from Gueguen for T-cells!

#### Gene Signatures
- See covid_cycling.R for most up-to-date version of Seurat's AddModuleScore()
- Beware signatures over unsuperivsed clustering! 
	- See Cheng 2021 macrophage section, do clusters/cluster integrity/cluster metadata first, then assign biology
- Gueguen (2021) utilised many useful gene signatures (from other datasets too) - see Table S3
- Sade-Feldman et al. (2018): signatures from melanoma patients responding (good response) versus not responding (bad response) to ICB. 
	- Nice - but need an even more specific 'tumour-reactive signature'? - does it exist? surely just active/cyotoxic?
- Metabolic signatures
- 	Definetely worth thinking about! See Gueguen (2021), can just be a mention but certainly useful
	- 	Very exhaustive (pun intended) signatures which they use to really interrogae cluster phenotype - imagine if they're used a better visualisation algo and had a more competant graph!
- Signatures are nice, but remember to be exploratory
	- i.e. "derive my own signatures", not simply test/try others
- See each dataset's .md/.txt file for full info
- On z-norm? don't want a gene w high expression to offset (although unlikely as 10X)
- Perhaps also make reduced Guo sigs, i.e. first 20 genes (MOST DEG'd) - think of 10X depth etc

#### Gene Set Variation Analysis
- https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-7

#### Classifier:
- Adam mentioned using datasets i.e. zavidij to train classifiers? did he mean that in a literal machine learning sense - if get reply about harmony counts, could ask in reply for advice (papers, best classifiers, etc)
- Fancy, but really, really do not see the point of something like this
- Adam reply:
	- "Anna (PhD student), has written an annotation pipeline and is working on classifying the cells using the human cell atlas and microarray data from PBMCs, then the next step will involve building her own classifier. The first stage would be curating public data that has been well annotated (Although I have a feeling the public data isn’t very good). Otherwise we will develop a classifier using our patient bone marrow samples, manually annotate them, then build a classifier from that."

#### Cell scores with UCell
- Suposedly more powerful than conventional Seurat/Scanpy approach


#### Dir of signature genes

- Published
- Curated databases
- Other tools

Make one long CSV, something like:
FirstAuthor / Cell_ID / CellType / Genes
i.e.
Guo18 / CD8-LAYN / ExhaustedCD8 / CD111,CD3,IFNH,QWERTO69,etc

Then have some tool to extract a list from the df entry of comma-sep values


#### Later doublet-detection…
could run a quick loop for cell-wise co-expression of distinct cell signatures in high-gene cells - i.e. see if one cells with many genes scores highly in both T-cells and DCs...
- ofc would then need to re-run clustering pipeline...


#### Gueguen has some nice marker genes too
- fgsea? general genes of interest? idk


#### Bayesian belief networks:
https://machinelearningmastery.com/introduction-to-bayesian-belief-networks/


#### SCENIC
- Computational method to infer Gene Regulatory Networks and cell types from single-cell RNA-seq data.
- https://github.com/aertslab/SCENIC
- This is a little more complex than the reglative simplicity of CellChat/Cellphone etc - employs actual RNA transcription factor pathways, etc 
    - 'REGULONS'
- Would actually suggest reading full paper if I want to understand, but then again I think if we treated this as sort of analagous to gene scores and simply use to infer pathway activity etc between conditions might be fine?
    - would need to go back to some primary lit, connect specific TFs and pathways to observed phenotypes, ideally also link to Chen's flow?
- First encountered in Ren 2021 (Zemin Zhang sc COVID19 paper) - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7857060/
    - could serve as a good resource for the sort of TFs etc we expect to see in single-cel data
    - however, would rely on immune subset reviews for specific TFs connected to biology!

Scenic workflow steps
0. To be ran on a high-performance system (i.e. the cluster)
1. Download a list of correct transcription factors, ranking/motif databases
2. Prep data: SCENIC ran on the uncorrected counts matrix of cells of interest
    - ensure cells have enough expressed/detected genes for ranking (define treshold manually), see Cellular Enrichment
3. Run pyscenic in shell using pythonic objects
    - outputs an adjascenices matrix for TF and targets
5. Regulons predictions
    - outputs a list of adjascenies between a TF and it's targets stored in MOTIFS_FNAME
6. Cellular enrichment 
    - outputs final scenic output
7. Visualisation of SCENIC AUC metrix output
    - a step in itself, output format complex

Runnning:
1. conda activate scenic_protocol




#### Cycling cells/clusters
Working with 2.2_Phenotypic isolated bone fide cycling clusters
- If I slice out a cycling cluster, can re-projection
- May be cycling-like cells / clusters within individual cell lineage sub-datasets - not a problem (will cluster as distinct) and may meld nicely with projected cells

- Cell cycle–related genes are expressed at high levels, and thus, cycling cells cluster independently in scRNA-seq analyses. 
- These clusters, however, may include cells that also bear underlying signatures from other clusters, which are in some way masked by highly expressed cell cycle genes. 
- Can investigate whether cells in cycling cluster are related to other clusters thusly:
	- Label transfer
	- Flow cytometry confirmation, i.e. Ki67, and co-expression with other lineage markers i.e. Trm
	- TCR clonotypes / TCR/CDR3 sharing

Approach also applicable to IFN signalling etc?
- May need a GEO IFN-stimulated gene list etc
- Annoyingly - may become relevant again upon re-examination... hope now




#### Immuno-Navigator - A database for gene coexpression in the immune system
- What it says on the tin
- https://genomics.virus.kyoto-u.ac.jp/immuno-navigator/
- Unsure how to use with large datasets

#### General
- HPSCs (in Zavidij for example)
- This is ~equiv to stroma in other tissue (i.e. local pops)
- Kwee thinks might be interesting, also wants Fibro I'm but unsure we'll get
- Mono -> Osteoclasts?





### Compositional analysis

#### MiloR
- can't do the pretty tripartide graphs in my head: outputs fc
- when assessing fold-change between multiple groups (i.e. hd, mgus, mm), how to then decide what is signif?
- repurpose the approches from scran DEGs? raises Q then of how to visualise etcs
- yes! so having run through there's def some thinking needed here - i.e. how to test between 3 groups (design matrix etc)
- re: FC, I do like the idea above of strategies from the scran testing framework - just a shame will lose beeswarm plots etc
- see here: https://r4ds.had.co.nz/model-basics.html#formulas-and-model-families
- miloR actually would work slightly better with large clusters, whose composition/phenotype can be investigated for condition-specific effects
    - to do this properly, need to be semi-pure AKA not containing contaminant pops (although suppose can visually ID and remove after)

- Cell-number ratio of clusters (across states, samples, etc)
- This was a specific context, but Cole Trapnell presented the idea of using an abundance matrix of cell types across samples into pesudotime for 'psuedostage' analysis


#### scCODA
Buttner 2020
- Total number of cells per sample is restricted in most single-cell technologies, implying that cell type counts are proportional in nature
- In turn, leads to a negative bias in cell type correlation estimation
    - For example, if only a specific cell type is depleted after perturbation, the relative frequency of others will rise. 
    - If taken at face value, this would lead to an inflation of differential cell types. 
    - Therefore, standard univariate statistical models that test compositional changes of each cell type independently may falsely deem certain population shifts as real effects, even though they were solely induced by the inherent negative correlations of the cell type proportions
- Draws inspiration from microbioma analysis, uses a Bayesian approach for cell type composition differential abundance analysis
    - models cell type counts with a hierarchical Dirichlet-Multinomial distribution that accounts for the uncertainty in cell type proportions and the negative correlative bias via joint modeling of all measured cell type proportions instead of individual ones
- Since compositional analysis always requires a reference to be able to identify compositional changes, scCODA currently relies on the manual specification of a reference cell type




#### Abundance II: ideas about Bailur analysis !!!! 30/04/21
If took the scran approach, could simply test neighbourhood composition for each disease stage against other two in turn

Then for each, have a UMAP of FC (up/down with Scale bar) for each disease 

Ie highlight what composition of each indiv stage looks like W supporting stats 

It could work?... idk 

As a good approach (working around miloR FC Q), how about start by introducing HD versus MM (support with lit and data i.e. Nakamura) - so general sense of TME
- Then move into MGUS

This idea actually has some interesting support: if we can properly support / characterise the healthy BM (given their flow work, similar cells in similar proportions to Oetjen would be sufficient, could even do a cheeky integration) then use this healthy as a basis to contrast with MGUS then MM
- I think showing similar phenotypes (projection/integration) and abundances  (overall and patient-wise) of Bailur HD with Oetjen would be super strong (esp as Oetjen has that nice spporting flow) - also that other BM Atlas I saw recently
- Returns to the idea of MM patients moving through an abnormal MGUS / being destined to progress
- Whereas MGUS patients have a different premalignant controlled disease - i.e. think on it let as a trajectory and more of a contrast of two disease states
- Pairwise contrasts here are suitable
- Built up from original idea of early-MM / SMM being this distinct entity to MGUS which has more in common with early malignancy then premalignancy
- "while true that an functionally all MM patients pass through MGUS and a temping hypothesis to consider a step-wise movement through discrete or continuous process of oncogenesis, recent work has demonstrated the differences between MGUS and MM at the level of progression suggestion patients (particularly smm patients...) may represent distinct disease states" (aka can be studied independently) "therefore contrasted HD to MGUS and MM independently, then contrasted the two disease conditions"
    - work around tools! Learning exercise ofc
- HD V MGUS
    - Less BM-like? Aka less myeloid cells etc? (Aka imagine oetjen v MGUS) - perhaps no! Simply a sampling thing, but then biology ofc is more CD4 cells


Thought - will miloR need all cells to have a label? i.e would only do on MM/HD clustered hoods etc

Scran's findMarkers()'s sumarry.logFC is defined as the log-fold change from the comparison with the lowest p-value
Combining p-values in accordance with Simes 1986???

post-CD4 paper idea: perhaps, yes, show a continuum of effectorness across the CD4 landscape in all cells
- but from there, anythink changed as a function of disease state? anything unique to MM (presumebly) and/or MGUS? shared/unique? could be fun



### Time-course Analysis

- Worth mentioning in single-cell systems cells are unsynchronized, which enables the capture of different instantaneous time points along an entire trajectory. Thus perhaps treating discrete time points as individuals is unproductive.

#### Compositional
- Poissan abundance testing?
    - edgeR on miloR? i.e. up the neighbourhoods stage then do testing? 
- Some interesting ResearchGate Q/As
- new scCODA ideas:
    - https://github.com/theislab/scCODA/issues/29

#### DEGs
- https://github.com/lengning/SCPattern
- edgeR tutorial:
    - but for bulk not single-cell? i.e. mean technical replicates



## 3. Ontogeny

# good overview/contrast !
https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#pseudotime-analysis
# !!!!!

key question: do I need to bother?
- how will looking at relateness better my hypothesis?
- is there a biological / abudnance / disease-driven / independent dataset-suggested REASON to look at this
- a case could be made that highly plastic subsets, i.e. CD4 T cells, naturally adopt a continuum of cell states and cannot be effectively phenotyped as discreate catagories
- if one possesses temporal data, i.e. time-course experiments, it may be worth looking at how different time points fall along a trajectory


### NBISweden: trajectory tutorials
- Have a dedicated set of trajectory tooold workshop:
    - https://github.com/NBISweden/workshop-scRNAseq/tree/master/labs/trajectory
- Have Seurat/Scanpy-specific tutorials too

### PAGA
- See 2.x_old_cluster.py 
- Also see here: https://github.com/theislab/paga
- blood/dahlin18, murine hematopoiesis (10x genomics) (https://nbviewer.jupyter.org/github/theislab/paga/blob/master/blood/dahlin18/dahlin18.ipynb)
    - `In [15]:` 'Increase the resolution in clusters 4 and 13 a bit more' - ???

### PHATE
- Get working and employ!

#### Ontogeny / relatedness 

i.e.
- Trajectory Analysis (& RNAv etc)
- PAGA
- TCR sharing between clusters
- Clustering as a function of resulution

I think working on sub-subclusters i.e. CD8 effector-memory-exhausted cells is totally fine, but it'll be cool to also to have something like ALL T-cells to serve as internal control
- RNAv on whole-dataset / subsets partially corrects this


Data prep
- Luecken 2020: "Regress out biological covariates only for trajectory inference and if other biological processes of interest are not masked by the regressed out biological covariate.

 Pseudotime Analysis / Trajectory Inference
- Notes:
	- Slingshot (Street et al, 2018), relatevely reliable tool: tutorial - https://bustools.github.io/BUS_notebooks_R/slingshot.html
	- TI with covariates: PhenoPath, Cambell (2018) - a downstream TI tool, once a rough pseudotime idea is established
	- See DownstreamAnalsis.md for lots of notes
	- Combine with (perfectly legal) subclusters, perhaps have transcriptomic info on different positions within a trajectory (however, expression may be sufficient)
	- Perhaps with myeloid only? More RNA
	- PAGA can confirm/corroberate pseudotime analysis: see Gueguen-2021 Fig 3
- Idea: scHaystack DEG clusters with pesudotime? 
- Broadly looking at genes that change as a function of pseudotime with clustering, be it DEGs of schaystack or not, is interesting
- If we have multiple pseudotime trajectories, can overlay them on one graph with expression of marker gene of interest overlayed (i.e. all T-cell dev pathways but then only have GZMK increasing for the trajectory of CD8 effectors, etc)


Single linear trajectory: SCORPIUS
- https://github.com/rcannood/SCORPIUS
- i.e. on only T cells
- however, should not solely expect single trajectories… may be multiple sources of Tex for example

!!!!
What is a trajectory tool that can interpret ANY trajectory, or a bifurcating one?
!!!!

- Older opinion:
    - IMO, unless doing like SMARTseq2 depth, TI may play the single of elucidating metadata-specific unique subpops (i.e., a terminally exhausted T cell only present in disease state)
        - However, with additional support i.e. RNA velocity, no reason couldn't also employ higher-evel analysis like PAGA, gene/cell networks
    - I now think that 10X depth TI may hold up, if biologically relevant, etc
    - Think about “meta-cells” - which are clusters of cells to reduce the effect of noise /effect of dropout in each single cell
    - Perhaps THIS could be combined with MiloR's ideas of many small clusters?



Pseudotime expression
- Developmental signatures? 
- Or simply use the reduction/gain of the cell-type signatures at the start and end of the trajectory?
- See 


Trajectories with covariates: Phenopath
- Campbell-2018
- Only covariate so far is disease stage - so perhaps if we get development NOT along a clear HD->MGUS->MM pesudotime axis (i.e. stages T-cell dev, excluding exhaustion, in every stage) we could use stage as a covar to investigate why
- I think this is a nice tool, but definetely requires an established trajectory first which a covariate is THEN tested on, with this initial trajectory being supported by bio, good pseudotime, and cytotrace/RNAv etc


PAGA
- Statistical analysis of connectivity in the k-nearest neighbor graph of cellcell expression similarities [partition-based graph abstraction (PAGA)
- See Wolf 2019

Silhouette Scores
- If the specific clusters represents a transitional state, then they should present lower clustering robustness, because cells from different clusters would enter this “state” and then differentiate into other states.
- Can quantify robustness of clustering using the silhouette score, which validates consistency by evaluating
how close each cell inside a cluster is to its neighboring cells within the same cluster (high score), compared with its neighboring cells in other clusters (low score).
- Transitional clusters would have a low silhouette score than their counterparts

Label Transfer
- A method introduced for the integration of multiple single-cell datasets
- Gueguen-2021: supplementary methods: 'Label transfer using a reference'
- T. Stuart, A. Butler, P. Hoffman, C. Hafemeister, E. Papalexi, W. M. Mauck III, Y. Hao, M. Stoeckius, P. Smibert, R. Satija, Comprehensive integration ofsingle-cell data. Cell 177, 1888–1902.e21 (2019).
- Cells from transitional clusters would show projection to the corresponding neighboring clusters
- **Would it work with scanpy's scanpy.tl.ingest?**
See swedish lab 06 celltype.ipynb - they define their own label transfer function! crikey

Cytotrace
- Novel: computational method that predicts the differentiation state of cells from single-cell RNA-sequencing data. 
- CytoTRACE leverages a simple, yet robust, determinant of developmental potential: the number of detectably expressed genes per cell, or gene counts. 
- https://cytotrace.stanford.edu/
- Ideally used on single cluster of cells assumed to lie on a trajectory (i.e. post-Slingshot/Monocle - for example Teff/Tex, etc)

Direct stratification
- Can trace a cell type's origin by using signature genes generated from each developmental subset to score all cells and the stratification of cells would give a reliable prediction for their origins
- See Cheng 2021 Fig S5B (works nicely!)
	- Defined signature genes by an absolute log2FC and adjusted Pvalue (variable depending on number of genes)
	- Calculated scores as the fraction of RNA in a cell belonging to signature genes
		- Employ of formulate, see paper
	- Plotted log10-transformed, with +ve and -ve axis to visualise positive/negative development


### Dana Pe'er Pseudotime Talk
- Can very easily get a pesudotime trajectory, but getting it right is difficult
- Requires deep and accurate psuedotime to properly understand biologu
- Introduces *Palantir*: complicated tool but very impressive, scanpy-wrapped
    - Git: https://github.com/dpeerlab/Palantir
    - Example: https://nbviewer.jupyter.org/github/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb
- Two large assumptuons of pesudotime
    - There is a starting point / cell
    - Cells move to more differntiated states (orientation)
- In vitro > vivo for all psuedotime analysis
- 

### RNA Velocity

#### Notes
- Zhang21 (CellPath) imply that one should perform RNAv first, so as to be able to get a broad look at the overall trends of a dataset's developmental trajectories
- Reveal areas of the data / some potential cell populations which might possess trajectories, agnostic to specific topologies (reducing relience on selecting specific TI tools for specific topologies)
- THEN, when have a general idea about structure of dataset, can select specific TI tools to address regions/pops etc (which correct topologies)
- Example large plot with velo info: DentateGyrus.ipynb in Velocyto example workflows

Marc input:
- velocity will tell you which state comes first in a tajectory 

Also see Zhang-2019-Landscape-and-dynamics-of-single-im.pdf
- Particularly Supp Fig 3E: RNAv for T cell exhaustion (ideal)
- It's just a couple of cells which fall between the two who DO have big RNA velocity to them that can infer about two states temporal order etc

NOTE: velocyto and scVelo are not equivalent
- Velocyto plot an arrow following each cell pointing to its future state,
- scVelo plot a streamline for a group of cells showing the dynamic trend in this group.

Dana Pe'er
- Vecolity: a smoothed tSNE/UMAP can somethings mislead the underlying *consistency* of velocity
    - need to really zoom in and see how consistent it looks
- Velocity is also only really reliable for the timeframe of splicing: around 8 hrs
    - defo NOT the same thing, but maybe look at (or think in framework of) Roger Geiger's preparedness fig 6/7? idk
    - Must have been a big T cell pseudotime/velocity development paper, no?
- Velocity may also simply not work in some situations: i.e. see haematopeisis
- Introduces *CellRank*, a velocity-integrated Palantir approach
    - Git: https://github.com/theislab/cellrank
    - Example: https://cellrank.readthedocs.io/en/latest/pancreas_basic.html
    - Notably has an internal consistency metric
- Can plot spliced expression of specific regulator (activ/repr) genes over pseudotime to reveal temporal dynamics



#### Velocyto

https://velocyto.org/velocyto.py/tutorial/index.html#
Four steps:
1. Installation: https://velocyto.org/velocyto.py/install/index.html#install
2. Running the CLI:https://velocyto.org/velocyto.py/tutorial/index.html#running-the-cli "Permalink to this headline")
    - This guide explains how to run the counting pipeline. It starts from .bam/.samfiles to yield a .loom file with counts divided in spliced/unspliced/ambiguous.
3. Analysis pipeline: https://velocyto.org/velocyto.py/tutorial/analysis.html
    - It's all here, and explained, but lacks visualisation / view of data, hence...
5. Exemplar notebooks: https://github.com/velocyto-team/velocyto-notebooks/tree/master/python
    - DentateGyrus.ipynb has one mega amazing graph

#### scVelo
- Alternative to velocyto

#### CellPath
Zhang 2021
https://github.com/PeterZZQ/CellPath (includes some example scripts)
- Very nice exploratory tool - to be applied to entire dataset post-RNAv analysis to acquire large and small trajectories on a dataset of different cells types alongside sub-trajectories, (smaller sub-flows within larger ones)
    - Specifically like the ability to detect very small trajectories, i.e. for small pops of interesting cells like NK cells and B cells etc
- Envisage this being applied early post-RNAv, then perhaps use specific tools like slingshot to support specific trajectories of interest
- However - while it may be useful to infer whereabouts subclustering could be done, would still do canonical clustering and assessment of cluster granularity / overall structure FIRST to avoid deliberating looking for subclusters via reclustering to specifically identify clusters which _may_ correspond to thus implied by this _single_ tool

#### Visualisation
veloViz: Atta 2021 (see folder in Coding&Bioinformatics)



#### Mitochondrial lineage tracing techniques
Dana Peer: "most exciting recent development"
- I think you just need to align reads to a dedicated mitochondrial genome?
- https://www.cell.com/cell/pdf/S0092-8674(19)30055-8.pdf
- https://elifesciences.org/articles/45105

MQuad Tool
- thread: https://twitter.com/yuanhuahuang/status/1376778067205443594?s=21
- preprint: https://www.biorxiv.org/content/10.1101/2021.03.27.437331v1
- "10x Genomics, as it had attracted a lot of attention recently owing to its accuracy of identifying CNVs and  clonal  structure  despite  its  low  coverage"
- doubts about MQuad! see https://twitter.com/caleblareau/status/1377081916738793474?s=21











## 4. Interactions

### Notes
- Javier: in Chris Tape's work, need all pops or communication won't make sense, i.e. can't miss any
- Nakamura et al 2020 review: real groundwork for the concilitation 

- At some point, here or in graphing, will need to re-import all finer-phenotypic/clustering labels etc and make a UMAP with all info

### CellPhoneDB
- co-expression filter (i.e. R-L pairs with co-expression, pull out) - for bulk (i.e. Nakamura) and sc?
- simple enough: see Bassez-2021 for a nice implmentation (see notes)

### NicheNet
- Predicts ligands driving the transcriptomic changes of target cells, **takes DEGs** (i.e. of clusters)
	- So perhaps could use a more specific/targeted follow-up (as not as broad)
- https://github.com/saeyslab/nichenetr
- Cheng 2021 employ it see their methods
- As takes DEGs, could perhaps select genes identified as differentially enriched by a more unusual temporal gene expression stats test, something better-suited/designed for time points over multiple sample
    - EdgeR example, something off of i.e. ResearchGate
- Has some extensive, nice-looking vignettes (very clear)

### CellChat
- cell-to-cell interaction tool
- https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
- good as general, perhaps confirm other results
- in pursuit of what goal though? no experimental condition from which to derive DEGs (apart from contrast to other datasets, but a bit whooly) - can it work on a purely observational-exploratory level?
- may be critical in solving the non-T cell fraction question
- Perhaps using simpler trajectory analysis (i.e. integrated our cells between HD/MGUS and MM cell and establish directionality) to establish only our cells in a postion (i.e. established in the TME, put them on a pedestool) then go in lots of detail about what they're doing i.e. cell signalling
    - Would be good if we could establish some nice T cell and myeloid pops 
    - If, in an unbiased way, we could establish directionality with the end-points having increased/decreased signalling, would be super nice
- (see above): The ultimate observation study? as in properly frame everything to support OUR data as important, then interrogate the hell out of it with tools like this
- This is no joke - need a serious understand of transcription factors and signalling networks to full optimise this tool. Can invisage prefacing this analysis with some big reviews (and diagrams) and then cementing results with more specific primary literature.
    - However... in context of our project, could be an oppurtunity to act as the Conciliator the corpus of signalling/cytokine/milueue literature in MM  (however, not SMM versus MM - a better MM sc dataset would be a good resource for contrast)
    - MM cytokinal mileue (20, Musolino 2017): https://pubmed.ncbi.nlm.nih.gov/29089667/
    - See also: Leblay 2020

#### CellChat run-through
- CellChat requires two user inputs: one is the gene expression data of cells, and the other is either user assigned cell labels 
- CellChat database(s) extensive, can use all or focus on a subset of biology (i.e. cell-cell secretory signalling) with one db
- No matter the faith in my hazy cell ID labels, at some point need to establish populations and work from there
    - can return to the umap later and look at expression of signalling molecules/genes etc
    - ALSO: one a communication network / concept is established, will then have to validate expression of these genes in another dataset
    - If we want to discuss signalling potential of mPCs, this will be essential
- `Preprocessing the expression data for cell-cell communication analysis`
    - this was a bit of a process - re-read later
- `Inference of cell-cell communication network`
    - Some options with calculating statisticals of interactions
    - "When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation."
        - Not problematic but definetely worth considering
- CellChat provides a function subsetCommunication to easily access the inferred cell-cell communications of interest.


### Other tools:
### CSOmap
Reconstruction of cell spatial organization from single-cell RNA sequencing data based on ligand-receptor mediated self-assembly
- https://www.nature.com/articles/s41422-020-0353-2
- Not downloaded PDF yet
### SoptSC
- https://academic.oup.com/nar/article/47/11/e66/5421812?login=true
- Suposedly prioritises intracellular L:R targets etc to verify
### Omnipath
- https://workflows.omnipathdb.org/
- I think it's a database vs a tool?
- See below...
### Squidpy R:L calling (leveraging Omnipath)
- https://squidpy.readthedocs.io/en/latest/auto_examples/graph/compute_ligrec.html
- Seems simple - could quickly run through




## 5. scTCRseq

TCR nomenclature reading
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2706371/


### Gueguen (2021)
- Some patient's scRNAseq libraries were seuqenced using the 10X Genomics Chromium Single Cell 3′ Solution, while for others RNA and TCR libraries were synthesised using the 10X Chromium Single Cell 5' V(D)J Enrichment Kit
	- They then integrated to correct for the impact 3' versus 5' variation
	- Depending on our libraries, may need to do?
- Clone sizes within clusters: 'expansion index'?
- Can calculate the conditional probability that if a particular TCR is observed in one cluster, then it will also be observed in the same or other clusters
	- i.e. quantify cluster relatedness from TCRs
	- Can then do unsupervised clustering of these results (Fig 4D) - identify shared clusters
	- Furthermore (it could be argued) high-intracluster TCR sharing highlights the relevance/correctness of our txnomic clustering strategy
	- NOTE: these kind of insights only hold up if are very confident in our txnomic clusters - i.e. haven't broken an in vivo irl populations in two computational subclusters and treating them as two. Lots of validation for both clustering & identification! Another reason to do full clustering QC i.e. silhouette score / function of resolution etc
- Computed the transition index using the single T-cell analysis by Rna-seq and Tcr TRACking (STARTRAC) method (33), as the likeliness of a cluster to share clones (fig. S4B)
	- Supervised representation of the computed probabilities of shared TCRs among clusters (Fig. 4F) 
	- Had a nice bit of text associated with it too - interpretation!
- Also did 'equilibrated sharing'? See FigS4D.
- When assessing preferential TCR sharing of different clusters, control for larger size of these specific clones (Fig S5F)

Remember to do all patient-wise as well!

### Tommy Tang notes
- https://github.com/crazyhottommy/TCR-BCR-seq-analysis
- Several tools and tutorials

Such as…

### Scirpy
- https://github.com/icbi-lab/scirpy
- Scanpy wrapped for TCRseq analysis, several tutorials:
- https://icbi-lab.github.io/scirpy-paper/wu2020.html
- https://icbi-lab.github.io/scirpy/tutorials/tutorial_3k_tcr.html

### TESSA
https://github.com/jcao89757/tessa
- Zhang 2021, TCR to expression signature tool
- Interested to see if this could be used with bulk data to get a better look at the reporteroire's phenotype



### Marc TCR Analysis Chat

#### TCR data
- Overall, a little easier than it looks
- Decombinator outputs a sample-wise list of TCR sequences with coresponding UMI counts
- You then have the option to aggregate by *translated* CDR3 sequences
    - Step which Marc always performs
    - This reduces the size of the dataset
    - No real advantage to retaining the entire TCR sequence of precide VDJ genes etc
    - This outputs a similar df of CDR3 seqs with UMI counts

#### Diverity measures
- One of the first steps you can always take is to generate diversity measures for each sample, Gini, Renyi etc
- Importantly: the size of the input will effect these (and other downstream) measures A LOT
    - High counts = harder to observe clonality
    - With samples with very high or very low counts, may need to totally isolate them (likely a sample prep/biological reason for standing out)
- For samples with very high counts, can take random sub-samples (weighted by size of overall sample) to try and replicate other samples: but lose indo and generally smooth slopes

#### Diversity and Clonality
- The two metrics are essentially opposite:
    - clonality/reduced diversity is indicate of Ag exposure & clonal responses

#### TCRs Across Samples
- i.e. sharing

#### Clustering
- Can also essentially always be done
- But also hard to derive a metric as equally impact by sample sizes
- Ideally need a good question to ask, ideally supported by multiple samples (i.e. elegant experimental design, pre-/post-relapse etc)
- Simple kernel script that does clustering will be sent / found on Chain RDSS - tools like iGraph etc simply for visualisation

#### Interface of single and bulk TCR-seq
- No real gold standard here, not a lot of work done here (so room to build!)
- TESSA is a novel concept but again totally new
- Can simply transfer labels
    - i.e. view bulk clustering, overlay phenotypes (perhaps of both a/b?)
    - A gradient of phenotypes? Expanded/clustered phenotypes versus others? 
    - With luck, will be able to get some insight into clonality / expansion of every single scTCR - so can classify levels of expansion/diversity in indivdidual clusters/regions
- Might need some more fundemental stuff, i.e. abundance/clonality across samples using both data sources
- Also fundemental - what of the TCR a/b sequences? How to use sc to inform bulk?
- Also worth mentioning all metrics discuss to far are purely on bulk - scTCR by 10X insufficient numbers for these approaches 
- Bulk alone can yield a ton of info for hypothesising, so perhaps simply need to build on that?



## 6. Regulatory Network Analysis

### Ken Olive thread
- https://twitter.com/KenOliveLab/status/1321451471515254784

the print:
- https://www.biorxiv.org/content/10.1101/2020.10.27.357269v1

What is “regulatory network analysis”? At it's heart, it's a way of extracting more useful information from expression profiles. A fundamental flaw of differential gene expression (DGE) analysis is the assumption that each gene is independent. Biologist KNOW this is not true!

DGE treats all 25K detectable genes as SEPARATE variables, performs 25K T-tests, and then slaps on a multiple hypothesis correction to make the statisticians less dyspeptic. There is no consideration of the relationships between genes! We KNOW genes are co-regulated in SETS by transcription factors and other REGULATORY FACTORS
- Other limitations of expression analysis: 
1. high variance for gene expression measurements
2. high dimensionality (25K!) means a big N in the denominator
3. Diff. expression does not suggest causality

RNA expression does not tell you anything about protein activity… or does it? You need an activity assay, right? What if I told you I could quantitatively measure enzyme activity using RNA expression data
- Transcription factors are simply enzymes that alter the abundance of RNA transcripts for specific target genes. 
- So RNAseq can quantify the products of the transcriptional enzymatic reaction, and therefore it can be used as an ASSAY to measure their PROTEIN ACTIVITY
- you need to know the targets for each regulatory factor. And that’s the trick. Because if you use databases or services that collate this info from literature, you mix together info from lots of cell types and are limited to what others have already found

The ARACNe algorithm, which uses information theory principles to computationally deduct target genes for regulatory factors, DE NOVO, using only gene expression profiles as input. The result is a global list of transcriptional relationships- a “regulatory network”.
- With a regulatory network as a scaffold, one can transform an expression profile of 25K genes to an ACTIVITY profile for perhaps 1500 regulatory factors. 
- Benefits:
1. Each activity measurement incorporates the expression of 100s of targets – way more precise/lower variance.
2. Notable dimensionality reduction
3. Incorporates relationships between genes
4. There is built-in mechanistic info. If you observe a difference in TF activity between two groups, they could be causal and you know their potential target gene effectors
- here? : https://github.com/califano-lab/ARACNe-AP

Critical piece here is the context specificity. The set of genes regulated by a TF in one cell type is hugely different than in another cell type

### Takeaways
Only thing scaring me off here is the fact there's so much literature already built up around concept of DEGs (at times supported with solid wet lab work)
- To count:two sections: 'crude' initial DEGs connected to existing sc lit, then more detailed Reg Network Analysis connected directly to wet lab and biology

Even if don't buy into the tool, criticism of DEGs noteworthy











## x. Graphics

Like the idea of having my own classes/data structures here :D

An issue here is that EVERY fucking package and module has it's own graphical functions

Will need to essentially export the raw data from everything
- However! there is a no need to really have it outside of a .h5ad or .RDS file, as long as it's stored
- So for this script/pipeline/stage - only really need the import data

Improving tSNE with initialisation parameters
- https://github.com/berenslab/rna-seq-tsne
- When one runs UMAP on just T-cells, seems a bit too blobby and stays as one large visual cluster, so tSNE may be clearer

Good notes:
https://twitter.com/aaronquinlan/status/1366073521642692608

Using PCA
- Think on PCA: this is a good way to visualise the "actual" data (even tho batch-corr)
- That is, a PCA with high % of variation explained by PCs can give an insight into underlying structure
- Maybe on a reduced number of cells - would be useful on just T cells etc

### Scanpy plotting
If I want to take advantage of scanpy plotting API, can manually set the adata.uns color palettes for different variables to. my own colors:
- matplotlib palettes: https://matplotlib.org/stable/tutorials/colors/colormaps.html
- `sc.set_figure_params(color_map=p)`
- See examples - need something which starts dark (for outline of umap) but then changes color quickly, but retains a recognisable very high color
    - `hot` and especially `afmhot` do this very well, but maybe starting at darkest black a bit much? A grey/bley might be nicer
- Visualise UMAP with normalsied data? Purely to allow visualisation
    - sc.pl.umap(dataZ, color=['CD8B'], use_raw=False, legend_loc='on data', legend_fontoutline=2)

### GIFs in matplot lib
https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30
