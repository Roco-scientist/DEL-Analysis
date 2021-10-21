# DEL-Analysis
DNA encoded library analysis.  This is companion software to <a href=https://github.com/Roco-scientist/NGS-Barcode-Count>NGS-Barcode-Count</a> for outputing analysis and graphs.

## Table of Contents
<ul>
<li><a href=#installation>Installation</a></li>
<li><a href=#files-needed>Files Needed</a></li>
<li><a href=#run>Run</a></li>
<li><a href=#methods>Methods</a></li>
</ul>

## Installation

Anaconda python required for the instructions below

### Create a del environment and activate

```
conda create -n del python=3.9
conda activate del
```

### Install

From pypl:<br>

```
pip install delanalysis
```

From source:<br>

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

All code below is within python<br><br>

```
import delanalysis

# Import merged data output from NGS-Barcode-Count.  This creates a DelDataMerged object
merged_data = delanalysis.read_merged("test_counts.all.csv")

# zscore, then quantile_normalize, then subtract background which is 'test_1'
merged_data_transformed = merged_data.binomial_zscore().subtract_background("test_1")

# Create a 2d comparison graph between 'test_2' and 'test_3' in the current directory and with a low end cutoff of 4
merged_data_transformed.comparison_graph("test_2", "test_3", "./", 4)

# Creates a DelDataSample object from a single sample from the merged object
test_2_data_transformed = merged_data_transformed.sample_data("test_2")

# Create a 3d graph with each axis being a barcode within the current directory and a low end cutoff of 4
test_2_data_transformed.graph_3d("./", 4)

# Create a 2d graph within the current directory and a low end cutoff of 4
test_2_data_transformed.graph_2d("./", 4)

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
full_double_single_zscore.comparison_graph("test_2", "test_3", "../../test_del/", 0.002
```

### Working with sample data output

All code below is within python<br><br>

```
import delanalysis

# Import sample data output from NGS-Barcode-Count.  This creates a DelDataSample object
sample_data = delanalysis.read_sample("test_1.csv")

# zscore
sample_data_zscore = sample_data.binomial_zscore()

# Create a 3d graph with each axis being a barcode within the current directory and a low end cutoff of 4
sample_data_zscore.graph_3d("./", 4)

# Create a 2d graph within the current directory and a low end cutoff of 4
sample_data_zscore.graph_2d("./", 4)
```

### Resulting graphs

The actual graphs will be interactive HTML graphs with hover data etc. <br><br>

From comparison_graph()<br>

![ "delanalysis.comparison_graph()" ](./comparison_graph.png)<br>

From graph_2d()<br>

![ "delanalysis.graph_2d()" ](./2d_graph.png)<br>

From graph_3d()<br>

![ "delanalysis.graph_3d()" ](./3d_graph.png)<br>

## Methods

### delanalysis methods to import data

<table>
<tr>
<th>Method</th> <th>Description</th>
</tr>
<tr>
<td>read_merged(file_path)</td> <td>Creates a DelDataMerged object which can use the methods below</td>
</tr>
<tr>
<td>read_sample(file_path)</td> <td>Creates a DelDataSample object which can use the methods below</td>
</tr>
<tr>
<td></td> <td></td>
</tr>
</table>

### Common to merged data and sample data

Used with either delanalysis.read_merged() or delanalysis.read_sample() objects

<table>
<tr>
<th>Method</th> <th>Description</th>
</tr>
<tr>
<td>building_block_columns()</td> <td>returns all column names which contain building block info</td>
</tr>
<tr>
<td>data_columns()</td> <td>returns all column names which contain data</td>
</tr>
<tr>
<td>data_descriptor()</td> <td>Returns data_type with underscores for file output</td>
</tr>
<tr>
<td>data_type</td> <td>The data type of the DelData</td>
</tr>
<tr>
<td>to_csv(out_file)</td> <td>Writes the DelData object to the out_file in csv format</td>
</tr>
<tr>
<td>zscore(inplace=False)</td> <td>z-scores the data</td>
</tr>
<tr>
<td>binomial_zscore(del_library_size, inplace=False)</td> <td>z-scores the data using the binomial distribution standard deviation</td>
</tr>
<tr>
<td>binomial_zscore_sample_normalized(del_library_size, inplace=False)</td> <td>z-scores the data using the binomial distribution standard deviation and normalizes by sqrt(n). See: <a href=https://pubs.acs.org/doi/10.1021/acscombsci.8b00116>Quantitative Comparison of Enrichment...</a></td>
</tr>
<tr>
<td>enrichment(del_library_size, inplace=False)</td> <td>count * library_size/ total_counts</td>
</tr>
<tr>
<td>update_synthon_numbers(unique_synthons_per_barcode: List[int])</td> <td>The number of unique synthons is inferred by the total uniques found in the data.  These numbers can be updated with this function</td>
</tr>
</table>

### Merged data

Used with delanalysis.read_merged() which creates a DelDataMerged object

<table>
<tr>
<th>Method</th> <th>Description</th>
</tr>
<tr>
<td>quantile_normalize(inplace=False)</td> <td>quantile normalizes the data</td>
</tr>
<tr>
<td>sample_enrichment(inplace=False)</td> <td>(sample_count/total_sample_count)/(non_sample_count/total_non_sample_count).  Still experimental as if the count only happens in one sample, a div 0 error occurs</td>
</tr>
<tr>
<td>subtract_background(background_name, inplace=False)</td> <td>subtracts the background_name sample from all other samples</td>
</tr>
<tr>
<td>reduce(min_score, inplace=False)</td> <td>Removes all rows from the data where no samples have a score above the min_score</td>
</tr>
<tr>
<td>merge(deldata, inplace=False)</td> <td>Merges DelDataMerged data into the current DelDataMerged object</td>
</tr>
<tr>
<td>sample_data(sample_name)</td> <td>Returns a DelDataSample object from the DelDataMerged object.  This is needed for the 2d and 3d graph</td>
</tr>
<tr>
<td>select_samples(sample_names: List, inplace=False)</td> <td>Reduces the data to the listed sample names</td>
</tr>
<tr>
<td>comparison_graph(x_sample, y_sample, out_dir, min_score=0)</td> <td>Outputs a comparison graph of x_sample vs y_sample names.</td>
</tr>
</table>

### Sample data

Used with delanalysis.read_sample() which creates a DelDataSample object

<table>
<tr>
<th>Method</th> <th>Description</th>
</tr>
<tr>
<td>reduce(min_score, inplace=False)</td> <td>reduces the data to only data greater than the min_score</td>
</tr>
<tr>
<td>max_score()</td> <td>Returns the maximum score within the data</td>
</tr>
<tr>
<td>data_column()</td> <td>Returns the data column name</td>
</tr>
<tr>
<td>graph_2d(out_dir, min_score=0)</td> <td>Produces two subplot 2d graphs for the different barcodes of a DelDataSample.</td>
</tr>
<tr>
<td>graph_3d(out_dir, min_score=0)</td> <td>Produces 3d graphs for the different barcodes of a DelDataSample.</td>
</tr>
</table>
