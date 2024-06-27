test_that("bed12 works", {
  expect_equal(read.table(system.file("extdata",'test_bed12.txt',package="CIRCUS"),
                          header = T, stringsAsFactors = FALSE, sep = '\t'),
               CIRCUS(ref_gpf_path = system.file("extdata",'transcript2gene.txt',package="CIRCUS"),
                                 ref_path = system.file("extdata",'annotation.bed',package="CIRCUS"),
                                 bed_path = system.file("extdata",'bed12_coords.bed',package="CIRCUS"),
                                 bed6=FALSE))
})
test_that("bed6 works", {
  expect_equal(read.table(system.file("extdata",'test_bed6.txt',package="CIRCUS"),
                          header = T, stringsAsFactors = FALSE, sep = '\t'),
               CIRCUS(ref_gpf_path = system.file("extdata",'transcript2gene.txt',package="CIRCUS"),
                                ref_path = system.file("extdata",'annotation.bed',package="CIRCUS"),
                                bed_path = system.file("extdata",'bed6_coords.bed',package="CIRCUS"),
                                bed6=TRUE))
})
