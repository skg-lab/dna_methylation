# tag-count within bins for bismark file

山田作

```
time python bismark_tag_count.py test/test_1000line_CpG_naive_pTreg.sort.bismark.cov test/test_1000_chrominfo_mm10.bed test/test.out
```

# methyl_diffのラッパー

## install

https://github.com/EmanueleRaineri/methyl_diff

> Raineri, E., Dabad, M. & Heath, S. A Note on Exact Differences between Beta Distributions in Genomic (Methylation) Studies. PLoS One 9, e97349 (2014).

```bash
$ git clone https://github.com/EmanueleRaineri/methyl_diff.git
$ cd methyl_diff
$ make
```

makeがmacでは通らなかった。(おそらくpyenvのせい。pyenvは使わないほうがやはりいい。)そのかわり以下でコンパイルすることも出来る。

```bash
gcc methyl_diff.c -o methyl_diff
```

```bash
./methyl_diff < diff.in > out.txt
```

## usage

```
python hemi_methyl.py test/posi.txt test/nega.txt test/hemi.out.txt -t 8
```

Bonferroni補正のデフォルトalpha値は0.01
