Functional Annotation
This package of scripts is designed to conduct functional annotation of variants identified in genome wide association studies of human disease.
The annotation is tailored for specific tissues and is currently set up for respiratory disease (specifically asthma).
The scripts are designed to be run in the following sequence:EQTLannotation.py,3DSNP_Api.py,MetaAnalysis.py
Users are required to download the database files referenced in the code from the various databases.The filenames are unchanged to aid identification of the correct files.
The output of the program will be a summmary of the genes identified as associated with the variants and the nature of the interaction identifed.A StringDB protein interaction network file and GO enrichment analysis will also be produced. 
