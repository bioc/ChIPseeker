# ChIPseeker 1.23.1

+ update GEO data (51079/762820 GSM) (2019-12-20, Fri)

# ChIPseeker 1.22.0

+ Bioconductor 3.10 release


# ChIPseeker 1.21.1

+ new implementation of `upsetplot` (2019-08-29, Thu)
  - use `ggupset`, `ggimage` and `ggplotify`
+ `subset` method for `csAnno` object (2019-08-27, Tue)

# ChIPseeker 1.20.0

+ Bioconductor 3.9 release

# ChIPseeker 1.19.1

+ add `origin_label = "TSS"` parameter to `plotAvgProf` (2018-12-12, Wed)
  - <https://github.com/GuangchuangYu/ChIPseeker/issues/91>
  
# ChIPseeker 1.18.0

+ Bioconductor 3.8 release

# ChIPseeker 1.17.2

+ add `flip_minor_strand` parameter in `getTagMatrix` (2018-08-10, Fri)
  - should set to FALSE if windows if not symetric
  
# ChIPseeker 1.17.1

+ fixed issue of `vennpie` by adding pseudo-count +1 (2018-07-21, Sat)
  - <https://www.biostars.org/p/326456/>

# ChIPseeker 1.16.0

+ Bioconductor 3.7 release

# ChIPseeker 1.15.4

+ If the required input is a named list and user input a list without name,
  set the name automatically and throw warning msg instead of error <2018-03-14,
  Wed>
    - <https://support.bioconductor.org/p/106903/#106936>
+ change `plotAvgProf`'s default y label <2018-03-14, Wed>
    - <https://github.com/GuangchuangYu/ChIPseeker/issues/76>
+ plotAnnoBar now visualize barplot according to the order of input list
  (y-axis) (2018-02-27, Tue)
    - <https://github.com/GuangchuangYu/ChIPseeker/issues/73>
+ follow renaming of RangesList class -> IntegerRangesList in IRanges v2.13.12
    - <https://github.com/GuangchuangYu/ChIPseeker/commit/b62d7922fb61e58620bbb685e4def4fb863c8e81>

# ChIPseeker 1.15.3

+ options to ignore '1st exon', '1st intron', 'downstream' and promoter
  subcategory when summarizing result and visualization (2018-01-09, Tue)
    - <https://support.bioconductor.org/p/104676/#104689>
+ throw msg of 'file not found and skip' when requested url is not available
  when downloading BED file from GEO (2017-12-28, Thu)
    - <https://support.bioconductor.org/p/104491/#104507>
+ bug fixed of getGene (2017-12-27, Wed)
