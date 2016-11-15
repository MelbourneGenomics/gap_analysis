Master Gap Analysis
=======================

This script analyses the Coverage output, i.e. *cov.gz* files from a bioinformatics analysis.
It discovers systematic gaps across multiple coverage files through calculating the difference
between the average mean and standard deviation at a specific base and comparing it to a 
coverage depth threshold. This identifies whether a base or several bases are considered to be 
a systematic gap and also whether it should be added to the "Master Gap" BED file outputted.


**Please Note**: The coverage files analysed should pertain to the specified TARGET_BED file 

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

Output
=======================
The output consists of a BED file with the headings:
chromosome, start position, end position, gene name, mean

Only the systematic gaps will be included in the bed file. 
If there are more than one adjacent bases considered systematic, they are treated as one
interval and the average mean for that interval is outputted in the master gap bed file.

Example output: 
```
chr1   325105  325132  LOC100133331    76.7407407407
chr1    721405  721436  Intergenic_2    22.2204301075
chr1    721441  721447  Intergenic_2    29.1111111111
chr1    762080  762093  LINC00115   24.5
chr1    762553  762571  LINC00115   22.2685185185
chr1    865702  865718  SAMD11  20.78125
chr1    874653  874662  SAMD11  26.6296296296
chr1    874829  874842  SAMD11  25.5384615385
chr1    876653  876688  SAMD11  20.6428571429
chr1    877473  877683  SAMD11  13.4055555556
chr1    877727  877762  SAMD11  18.4857142857
..
..
..
..
chrY    59343005    59343080    IL9R    12.32
```

