# Proteomics analysis pipeline - LSP
The proteomics analysis pipeline consists of a suite of tools that support the design and analysis of mass-spec based proteomics and phosphoproteomic measurements. In addition to custom scripts, the suite also incorporates pre-existing tools such as GSEA, VIPER and Hotnet in order to make network level inference of differential pathway activity.

# Using the pipeline
In the first example below, a TMT-labeled mass spec dataset is merged with user-provided metadata.  In addition, all outdated identifiers (Uniprot Id, Gene name) are identified and updated. The metadata file should be a csv file containing a column 'TMT_label' that correspond with the columns (tmt tags) in the dataset, and a second column caled 'Sample' which is a user provided name for the sample. The metadata file may optionally also contain additional columns that describe other categorical features of the samples.
```
import pandas as pd
from msda import preprocessing
meta_df = pd.read_csv(metadata_file.csv)
df = preprocessing.pd_import(data_file, meta_df)
```
If the experiment samples are split across two or more 10-plex experiments and a bridge sample is present in each 10-plex, then the datasets can be normalized based on the bridge sample and merged into a single dataset before subsequent analysis
```
df = preprocssing.pd_import([data_file1, data_file2], meta_df)
```
The next example shows how principal component analysis or hierarchichal clustering of the dataset can be done. 
``` 
from msda import clustering
clustering.pca(df, meta_df, num_components=2, label=None, output_file='figure.png')
clustering.hierarchical_clustering(df, output_file='figure.png')
```
If the metadata file contains additional columns with categorical information about the samples, then the PCA plots can also depict the categories using different colors. In order to do so, simply pass the name of the column in the metadata file to the argument 'label'.
``` 
clustering.pca(df, meta_df, num_components=2, label='tumor_subtype', output_file)
```

Mass-spec based proteomics measurements provide relative intensity measurements for each protein i.e intensities for a given protein can be compared between samples, but intensities between proteins cannot be comapred. Intensity based absolute quantifcation (iBAQ) is a  normaliation methods that enables comparison across proteins. The example below show how to map the dataset onto iBAQ space.
```
from msda import ups2_calibration as uc
df_ibaq = uc.compute_ibaq_dataset(df, organism='human', meta_df)
```

If UPS2-based calibration has been performed for one of the samples in the dataset (for eg:'Bridge'), then the intensity measurements can be mapped to concentrations (or copies/cell).
```
df_ups2 = pd.read_csv(ups2_dataset_file)
df_abundance = uc.compute_abundance(df, df_ups2, sample='Bridge')
```
GSEA

VIPER

