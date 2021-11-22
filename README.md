# DEL-Analysis
DNA encoded library analysis.  This is companion software to [NGS-Barcode-Count](https://github.com/Roco-scientist/NGS-Barcode-Count) for outputing analysis and graphs.

## Table of Contents
- [Installation](#installation)
- [Files Needed](#files-needed)
- [Run](#run)
- [Methods](#methods)

## Installation

Anaconda python required for the instructions below

### Create a del environment and activate

```
conda create -n del python=3.9
conda activate del
```

### Install

From pypl:  

```
pip install delanalysis
```

From source:  

```
git clone https://github.com/Roco-scientist/DEL-Analysis.git
cd DEL-Analysis
pip install --use-feature=in-tree-build .
```

## Files Needed
Output files from NGS-Barcode-Count

## Run

### Start

```
conda activate del
python
```

### Working with merged data output

All code below is within python  
  

```
import delanalysis

# Import merged data output from NGS-Barcode-Count.  This creates a DelDataMerged object
merged_data = delanalysis.read_merged("test_counts.all.csv")

# zscore, then quantile_normalize, then subtract background which is 'test_1'
merged_data_transformed = merged_data.binomial_zscore().subtract_background(background_name="test_1")

# Create a 2d comparison graph between 'test_2' and 'test_3' in the current directory and with a low end cutoff of 4
merged_data_transformed.comparison_graph(x_sample="test_2", y_sample="test_3", out_dir="./", min_score=4)

# Creates a DelDataSample object from a single sample from the merged object
test_2_data_transformed = merged_data_transformed.sample_data(sample_name="test_2")

# Create a 3d graph with each axis being a barcode within the current directory and a low end cutoff of 4
test_2_data_transformed.graph_3d(out_dir="./", min_score=4)

# Create a 2d graph within the current directory and a low end cutoff of 4
test_2_data_transformed.graph_2d(out_dir="./", min_score=4)

# Can all be done in one line
delanalysis.read_merged("test_counts.all.csv").binomial_zscore().subtract_background("test_1").sample_data("test_2").graph_3d("./", 4)


# Create a comparison graph for tri, di, and mono synthons
full = read_merged("../../test_del/test.all.csv")
double = read_merged("../../test_del/test.all.Double.csv")
single = read_merged("../../test_del/test.all.Single.csv")
full_double = full.concat(double)
full_double_single = full_double.concat(single)
full_double_single_zscore = full_double_single.binomial_zscore_sample_normalized()
full_double_single_zscore.subtract_background("test_1", inplace=True)
full_double_single_zscore.comparison_graph("test_2", "test_3", "../../test_del/", 0.002)
```

### Working with sample data output

All code below is within python  
  

```
import delanalysis

# Import sample data output from NGS-Barcode-Count.  This creates a DelDataSample object
sample_data = delanalysis.read_sample("test_1.csv")

# zscore
sample_data_zscore = sample_data.binomial_zscore()

# Create a 3d graph with each axis being a barcode within the current directory and a low end cutoff of 4
sample_data_zscore.graph_3d(out_dir="./", min_score=4)

# Create a 2d graph within the current directory and a low end cutoff of 4
sample_data_zscore.graph_2d(out_dir="./", min_score=4)
```

### Resulting graphs

The actual graphs will be interactive HTML graphs with hover data etc.  
  
From comparison_graph()  

![ "delanalysis.comparison_graph()" ](./comparison_graph.png)  

From graph_2d()  

![ "delanalysis.graph_2d()" ](./2d_graph.png)  

From graph_3d()  

![ "delanalysis.graph_3d()" ](./3d_graph.png)  

## Methods

### delanalysis methods to import data

|Method |Description|
|-------|-----------|
|read_merged(file_path) |Creates a DelDataMerged object which can use the methods below|
|read_sample(file_path) |Creates a DelDataSample object which can use the methods below|

### Common to merged data and sample data

Used with either delanalysis.read_merged() or delanalysis.read_sample() objects  

|Method |Description|
|-------|-----------|
|building_block_columns() |returns all column names which contain building block info|
|data_columns() |returns all column names which contain data|
|data_descriptor() |Returns data_type with underscores for file output|
|data_type |The data type of the DelData|
|to_csv(out_file) |Writes the DelData object to the out_file in csv format|
|zscore(inplace=False) |z-scores the data|
|binomial_zscore(del_library_size, inplace=False) |z-scores the data using the binomial distribution standard deviation|
|binomial_zscore_sample_normalized(del_library_size, inplace=False) |z-scores the data using the binomial distribution standard deviation and normalizes by sqrt(n). See: [Quantitative Comparison of Enrichment...](https://pubs.acs.org/doi/10.1021/acscombsci.8b00116)|
|enrichment(del_library_size, inplace=False) |count * library_size/ total_counts|
|update_synthon_numbers(unique_synthons_per_barcode: List[int]) |The number of unique synthons is inferred by the total uniques found in the data.  These numbers can be updated with this function|

### Merged data

Used with delanalysis.read_merged() which creates a DelDataMerged object

|Method |Description|
|-------|-----------|
|quantile_normalize(inplace=False) |quantile normalizes the data|
|sample_enrichment(inplace=False) |(sample_count/total_sample_count)/(non_sample_count/total_non_sample_count).  Still experimental as if the count only happens in one sample, a div 0 error occurs|
|subtract_background(background_name, inplace=False) |subtracts the background_name sample from all other samples|
|reduce(min_score, inplace=False) |Removes all rows from the data where no samples have a score above the min_score|
|merge(deldata, inplace=False) |Merges DelDataMerged data into the current DelDataMerged object|
|sample_data(sample_name) |Returns a DelDataSample object from the DelDataMerged object.  This is needed for the 2d and 3d graph|
|select_samples(sample_names: List, inplace=False) |Reduces the data to the listed sample names|
|comparison_graph(x_sample, y_sample, out_dir, min_score=0) |Outputs a comparison graph of x_sample vs y_sample names.|

### Sample data

Used with delanalysis.read_sample() which creates a DelDataSample object  

|Method |Description|
|-------|-----------|
|reduce(min_score, inplace=False) |reduces the data to only data greater than the min_score|
|max_score() |Returns the maximum score within the data|
|data_column() |Returns the data column name|
|graph_2d(out_dir, min_score=0) |Produces two subplot 2d graphs for the different barcodes of a DelDataSample.|
|graph_3d(out_dir, min_score=0) |Produces 3d graphs for the different barcodes of a DelDataSample.|
