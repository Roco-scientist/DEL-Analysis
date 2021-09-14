# DEL-Analysis
DNA encoded library analysis.  This is companion software to DEL-Decode for outputing analysis and graphs.

## Table of Contents
<ul>
<li><a href=#installation>Installation</a></li>
<li><a href=#files-needed>Files Needed</a></li>
<li><a href=#methods>Methods</a></li>
<li><a href=#run>Run</a></li>
</ul>

## Installation
Work in progress

## Files Needed
Output files from DEL-Decode

## Methods
Work in progress

### Merged data

<table>
<tr>
<th>Method</th>
<th>Description</th>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>

### Sample data

<table>
<tr>
<th>Method</th>
<th>Description</th>
</tr>
<tr>
<td></td>
<td></td>
</tr>
</table>

## Run
All code below is within python<br><br>

### Working with merged data output

```
import delanalysis

# Import merged data output from DEL-Decode.  This creates a DelDataMerged object
merged_data = delanalysis.read_merged("test_counts.all.csv")

# zscore, then quantile_normalize, then subtract background which is 'test_1'
merged_data_transformed = merged_data.zscore().quantile_normalize().subtract_background("test_1")

# Create a 2d comparison graph between 'test_2' and 'test_3' in the current directory and with a low end cutoff of 4
delanalysis.comparison_graph(merged_data_transformed, "test_2", "test_3", "./", 4)
```

### Working with sample data output

```
import delanalysis

# Import sample data output from DEL-Decode.  This creates a DelDataSample object
sample_data = delanalysis.read_sample("test_1.csv")

# zscore
sample_data_zscore = sample_data.zscore()

# Create a 3d graph with each axis being a barcode within the current directory and a low end cutoff of 4
delanalysis.graph_3d(sample_data_zscore, "./", 4)

# Create a 2d graph within the current directory and a low end cutoff of 4
delanalysis.graph_2d(sample_data_zscore, "./", 4)
```

