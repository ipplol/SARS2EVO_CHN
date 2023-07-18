# SARS-CoV-2 evolution in China during 20220701 and 20230531
This repository provides data and scripts used for the research [研究]()
The workflow was designed to run on a [WSL-enabled](https://learn.microsoft.com/en-us/windows/wsl/) Windows platform.
As most scripts were written in C# and one needs to adjust the file path manually, the [Microsoft VisualStudio](https://visualstudio.microsoft.com/zh-hans/) or other compatible IDE is needed.

These two files in the [google drive](https://drive.google.com/drive/folders/1W7KXdd2Q1crVCrd1eTMjZw41G0HiVvt_)
need to be unzipped and saved in _./China220701_230531/Data_ before processing:
_global_assignments.pb.gz_
_public-latestCHNadded.metadata.tsv.gz_
This can be done in a WSL environment with _"tar -zxvf [filename]"_
## Downloading and deduplication of sequences collected in China
Chinese SARS-COV-2 sequences collected between 1 July 2022 and 31 May 2023 were retrieved from the [GISAID](https://gisaid.org/) and [NGDC](https://ngdc.cncb.ac.cn/ncov/?lang=en) databases as well as their metadata file. 

To deduplicate the same sequences, a comprehensive comparison of sequence name, genome, collection date, and submit laboratory was conducted using:
***preprocessing_sequence.py***

Daily case data were downloaded from the [OWID](https://github.com/owid/covid-19-data) website.

## Calculation of lineage distribution over time
The ***CalLineageDistribution*** was used to calculate the lineage distribution in China from the metadata file. Both input and output files are stored in _./China220701_230531/CHNSubtree/CHNLineage_

## Finding evolutionary clade in China
The globally public SARS-COV-2 mutation-annotated tree (MAT) and corresponding metadata file were downloaded from the [McBroome et al.](http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/) website.

The metadata file of the Chinese sequences downloaded was combined into the metadata file of the MAT  manually.

The Chinese sequences obtained were placed on the UShER tree download using the UShER package (Linux)

The result [protobuf](https://protobuf.dev/) file _global_assignments.pb_  was first transferred to a C# readable JSON file using [matUtils](https://usher-wiki.readthedocs.io/en/latest/QuickStart.html#matutils) (Linux) from the [UshER toolkit](https://usher-wiki.readthedocs.io/en/latest/).

_“matUtils extract -i ./global_assignments.pb -d [_YOUR DIRECTORY_] -j global_assignments.json”_

Then we search the MAT tree for all mutation events and Chinese evolutionary clades using ***FindCHNSubtreeFromMAT***. It will generate several files:

*global_assignments.json.line*: A reformated JSON file.
*global_assignments.json.mutevent*: All mutation events between nodes.
*global_assignments.json.nodeinfo*: The metadata of each node on the MAT.
*global_assignments.json.CHNClade*: All Chinese evolutionary clades.

The three main clades were extracted using [matUtils](https://usher-wiki.readthedocs.io/en/latest/QuickStart.html#matutils) (Linux) :
_"matUtils extract -i ./China220701_230531/Data/global_assignments.pb -I node_1084141 -d ./China220701_230531/CHNSubtree -o node_1084141.pb"_

And then visualized with [usher_to_taxonium and Taxonium](https://usher-wiki.readthedocs.io/en/latest/tutorials.html#using-taxonium-to-visualize-phylogenies).

## Mutation events identification and incidence estimation
Meanwhile, the mutation incidence files have also been generated and saved at: _./China220701_230531/CHNSubtree/MutIncidence_

We then combine the incidence of clade BA.5.2.48 and BA.5.2.49 manually.

The incidence files were annotated to RBD amino acid mutations using ***Mutincidence2RBDincidence***.

The mutation incidence files were combined into a matrix using ***BuildingIncidenceMatrix***, and be saved at:
_./China220701_230531/ChinaVSAbroad/incidenceCorr_

Cosine similarity was calculated based on the R package.

The mutation event genome distribution was calculated using ***MutincidenceGenomeDistribution***.

The cumulative event fraction (P50) was calculated manually through [Microsoft Excel](https://www.microsoft.com/zh-cn/microsoft-365/excel).

## Calculation of the ACE2 binding affinity and the Escape score
The raw ACE2 binding and RBD expression score of BA.2 from [Bloom et al.](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1010951) was added together and saved in the _./China220701_230531/Data/bindExpr.txt_ as the ACE2 binding affinity changes in the paper.

The antibody spectrum, neutralizing activity, antibody epitope group, and raw mutation escape score were obtained from [Yunlong et al.](https://www.biorxiv.org/content/10.1101/2023.05.01.538516v2)'s studies. We combined the raw mutation escape score with antibody-neutralizing activity using ***CalMutEscapeScore*** and saved at _./China220701_230531/Data_, which contains several files:

_EscapeScore_PKU_NEW_BA.5.single.txt_: The mutation effect to BA.5 antibodies.
_EscapeScore_PKU_NEW_BA.5.12.txt_: The mutation effect to BA.5 antibodies of each epitope.
_EscapeScore_PKU_NEW_BA.5.12EscOnly.txt_: The mutation effect to BA.5 antibodies of each epitope without antibody-neutralizing activity.

The escape mutations of each epitope and their distribution were calculated using ***CalRBDmutEpitopeDistribution***.
Which output a _./China220701_230531/Data/MutationEsc12_NEW_BA5.txt_ files that record if the mutation is an escape mutation to the 12 epitopes.
The distribution is saved in _./China220701_230531/Group12/CladeRBDMut12Group.txt_

Then, the median Escape score and ACE2 binding changes of each lineage was calculated using ***SummaryMutincidenceACE2ESCMedian*** and saved in the _./China220701_230531/ChinaVSAbroad/BA5SubLineage/ACE2ESCMedian.txt_

To fit a normal distribution of the escape mutation proportion, the R package [fitdistrplus](https://cran.r-project.org/web/packages/fitdistrplus/index.html) was used.
_"FIT <- fitdist(x, "norm")"_

## Finding all possible one-step RBD mutations

Finding of all possible one-step RBD mutations was performed using ***SpikeOneStepMut***, it uses the spike gene sequence in the file _./China220701_230531/Group12/ACE2/BA5.Spike.fa_ as the reference, and output all possible mutations into file _./China220701_230531/Group12/ACE2/BA5.Spike.AvailableAA_

## Calculation of immune evasion of different strains
The ***ImmuneEscapeCalculator*** was used to estimate the immune evasion of different strains, which was based on the same formula as Bloom et al.'s  [escape calculator for the SARS-CoV-2 RBD](https://github.com/jbloomlab/SARS2-RBD-escape-calc/).  Input files were the sequences collected in China and abroad, extracted from _./China220701_230531/data/global_assignments.json.nodeinfo, and stored in the _**.MutDate.txt_

## Calculation of mutation R403K distribution over time
The ***Cal403KDistribution*** was used to calculate the lineage distribution of RBD mutation R403K. Results were saved into file _./China220701_230531/ChinaVSAbroad/403KDistribution/Global_History_AllLineage_403K_Distribution.tsv_