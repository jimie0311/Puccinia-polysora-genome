# The scripts for calculating HiC links

```
Conda install HiCexploer
#Follow the pipeline in Figure-1B to obtain hicpro_output, Note: the reference genome should be hapA+hapB
hicConvertFormat --matrices ./hicpro_output/hic_results/matrix/hapAB_all_150000_iced.matrix --inputFormat hicpro --outputFormat h5 --outFilename test.150000.matrix --bedFileHicpro ./hicpro_output/hic_results/matrix/hapAB/raw/150000/hapAB_all_150000_abs.bed
hicConvertFormat --matrices test.150000.matrix.h5 --inputFormat h5 -o 150000.tsv --outputFormat ginteractions
```

The informative data are in output file 150000.tsv.tsv. I list few lines. The last column represts hic contact. The %hic link of chr01A = sum(chr01Avschr01A)/sum[(chr01Avschr01A)+(chr01AvsChr01B)]*100%
```
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10050000        10200000        531.318121
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10200000        10350000        397.396175
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10350000        10500000        113.659477
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10500000        10650000        93.582654
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10650000        10800000        35.045091
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10800000        10950000        31.488322
Ppz_chr01A      10050000        10200000        Ppz_chr01A      10950000        11100000        34.151366
Ppz_chr01A      10050000        10200000        Ppz_chr01A      11100000        11250000        30.007855
Ppz_chr01A      10050000        10200000        Ppz_chr01A      11250000        11400000        22.104703
Ppz_chr01A      10050000        10200000        Ppz_chr01A      11400000        11550000        21.227431
Ppz_chr01A      10050000        10200000        Ppz_chr01A      11550000        11700000        23.082917
Ppz_chr01A      10050000        10200000        Ppz_chr01A      11700000        11850000        10.297713
Ppz_chr01A      10050000        10200000        Ppz_chr01A      11850000        12000000        6.888999
Ppz_chr01A      10050000        10200000        Ppz_chr01A      12000000        12150000        6.095324
Ppz_chr01A      10050000        10200000        Ppz_chr01A      12150000        12300000        11.334333
Ppz_chr01A      10050000        10200000        Ppz_chr01A      12300000        12450000        10.570572
Ppz_chr01A      10050000        10200000        Ppz_chr01A      12450000        12600000        11.983396
Ppz_chr01A      10050000        10200000        Ppz_chr01A      12600000        12750000        19.274579
```