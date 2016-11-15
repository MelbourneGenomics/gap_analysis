Master Gap Analysis
=======================

This script analyses the Coverage output, i.e. *cov.gz* files from a bioinformatics analysis.
It discovers systematic gaps across multiple coverage files through calculating the difference
between the average mean and standard deviation at a specific base and comparing it to a 
coverage depth threshold. This identifies whether a base or several bases are considered to be 
a systematic gap and also whether it should be added to the "Master Gap" BED file outputted.

Usage
=======================
```
usage: master_gaps.py [-h] --target_bed TARGET_BED --covfile COVFILE
                      [--threshold THRESHOLD] [-s SAVE_DATA_FOLDER]
                      [-l LOAD_DATA_FOLDER] [-o OUTPUT] [--skip]
```

This will output the master gap bed file to stdout 
```
python master_gaps.py --target_bed EXOME_CAPTURE.bed --covfile coverage_files.txt  
```

This will save output master gap bed file to an output location
```
python master_gaps.py --target_bed EXOME_CAPTURE.bed --covfile coverage_files.txt --output /var/tmp/master_gaps.bed
```

This will save the python and numpy data to a location to be loaded for further analysis later
```
python master_gaps.py --target_bed EXOME_CAPTURE.bed --covfile coverage_files.txt --save_data_folder /var/tmp/covdata --output /var/tmp/master_gaps.bed
```

This will load the python and numpy data and append information to the existing saved analysis then save the python and numpy data again
```
python master_gaps.py --target_bed EXOME_CAPTURE.bed --covfile new_coverage_files.txt --save_data_folder /var/tmp/covdata --load_data_folder /var/tmp/covdata --output /var/tmp/master_gaps.bed
```

coverage_files.txt
```
/var/tmp/coverage_file1.cov.gz
/var/tmp/coverage_file2.cov.gz
```

new_coverage_files.txt
```
/var/tmp/coverage_file3.cov.gz
/var/tmp/coverage_file3.cov.gz
```



