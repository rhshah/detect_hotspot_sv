# detect_hotspot_sv.py

Look for potential hotspot structural variants by applying filters to all SVs called by DELLY.

Some code adapted from [Ronak Shah's check_cDNA_contamination.py](https://github.com/rhshah/Miscellaneous/blob/master/check_cDNA_contamination.py).


### Usage:
```
usage: python detect_hotspot_sv.py [-h] [-s sample_name] [-t n]
                                [--iAnnotateSV iAnnotateSV.py]
                                [--genome {hg18,hg19,hg38}]
                                structural_variants.vcf
                                /path/to/output_directory/

positional arguments:
  structural_variants.vcf
                        VCF with deletion event calls
  /path/to/output_directory/
                        Directory to put output in

optional arguments:
  -h, --help            show this help message and exit
  -s sample_name, --sample-name sample_name
                        Prefix [default: filename of in_vcf]
  -t n, --threshold n   Output pseudogenes must have at least this many
                        qualifying deletion events. [default: 2]
  --iAnnotateSV iAnnotateSV.py
                        path to iAnnotateSV.py script [default: /home/shahr2/g
                        it/iAnnotateSV/iAnnotateSV/iAnnotateSV.py]
  --genome {hg18,hg19,hg38}
                        Reference genome version to use during annotation
                        [default: hg19]
```


### Example usage:
```
$ python detect_hotspot_sv.py example.vcf ./ -s my_sample
```


### Output:
- `<output directory>/<sample name>.hotspot_sv.txt`: a one-line, tab-delimited list of sample name and its processed pseudogenes
- Intermediate files in `<output directory>/scratch/`


### Requires:
- Python 2
- pandas
- gridmap
- [iAnnotateSV](https://github.com/rhshah/iAnnotateSV) and its required modules
- [iCallSV](https://github.com/rhshah/iCallSV)
- [PyVCF](https://github.com/jamescasbon/PyVCF)

## Script process

1. Filter the input VCF for PRECISE deletions with paired read count > 10
2. Convert the filtered VCF to tab-delimited input format for iAnnotateSV, using iCallSV.dellyVcf2Tab
3. Annotate using iAnnotateSV
4. For each gene, count the number of events that aren't in exons, "in frame" or "out of frame".If a gene has more than the threshold number of events, then conclude it's a processed pseudogene.

# run_detect_hotspot_sv_analysis.py 

Wrapper to run detect_hotspot_sv.py 

### Usage:
```
python run_detect_hotspot_sv_analysis.py  -h

usage: run_detect_hotspot_sv_analysis.py  [options]

Run processed-psuedogene_analysis on selected pools/samples using MSK data

optional arguments:
  -h, --help            show this help message and exit
  -mif folders.fof, --metaInformationFile folders.fof
                        Full path folders and sample names of files.
  -qc /some/path/qcLocation, --qcLocation /some/path/qcLocation
                        Full path qc files.
  -P /somepath/python, --python /somepath/python
                        Full path Pyhton executables.
  -dhs /somepath/detect_hotspot_sv.py, --detecthotspotsv /somepath/detect_hotspot_sv.py
                        Full path processpsuedogene.py executables.
  -ias /somepath/iAnnotateSV.py, --iAnnotateSV /somepath/iAnnotateSV.py
                        Full path iAnnotate.py executables.
  -q all.q or clin.q, --queue all.q or clin.q
                        Name of the SGE queue
  -qsub /somepath/qsub, --qsubPath /somepath/qsub
                        Full Path to the qsub executables of SGE.
  -t 5, --threads 5     Number of Threads to be used to run process-psuedogene
                        analysis on
  -v, --verbose         make lots of noise [default]
  -o /somepath/output, --outDir /somepath/output
                        Full Path to the output dir.
```

### Example usage:
```
$ python run_detect_hotspot_sv_analysis.py  -mif /location/to/meatdata.txt -qc /dmp/qc/location/path -ppg /path/to/processed-pseudogenes/detect_hotspot_sv.py -ias /path/to/iAnnotateSV.py -q sge_queue_name -qsub path/to/qsub -t 1 -o /path/to/output/directory -v 
```

### Meta information file format

Tab separated file with header:

| Run  | MRN  | Mnumber  | LIMS_ID  | DMP_ASSAY_ID | 12245 |
|---|---|---|---|---|---|
|  PoolName | SomeID  | SomeID  |  SomeID |  SomeID | Flag(0,1)  |

From here:
The **Run** and **LIMS_ID** are required.