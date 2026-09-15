library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicRanges)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peaks <- GRanges(
    c("chr1", "chr2", "chr20"),
    IRanges(start = c(1000000, 3000000, 3000000), width = 1000)
)

test_that("geneChr / geneStrand are characters, not factor codes", {
    ## issues #233, #247:
    ## as.data.frame() turns 'seqnames' / 'strand' into factors and assigning a
    ## factor into mcols() dropped the class, so geneChr / geneStrand came out
    ## as integers (e.g. 1/2 instead of chr1/chr2 and +/-)
    pa <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                       level = "transcript", verbose = FALSE)
    m <- mcols(as.GRanges(pa))

    expect_true(is.character(m$geneChr))
    expect_true(is.character(m$geneStrand))
    expect_true(all(m$geneStrand %in% c("+", "-")))
    expect_true(all(grepl("^chr", m$geneChr)))

    ## numeric columns must stay numeric
    expect_true(is.integer(m$geneStart))
    expect_true(is.integer(m$geneEnd))
})

test_that("peaks without any feature in TxDb are reported, not silently dropped", {
    ## issue #251: peaks on contigs/scaffolds that carry no gene used to be
    ## dropped without any message at all
    peaks2 <- GRanges(
        c("chr1", "chr1_gl000191_random", "chrM"),
        IRanges(start = c(1000000, 5000, 5000), width = 1000)
    )

    expect_warning(
        annotatePeak(peaks2, TxDb = txdb, tssRegion = c(-3000, 3000),
                     level = "transcript", verbose = FALSE),
        "peaks were dropped"
    )
})
