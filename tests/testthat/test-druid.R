test_that("bed12 works", {
  expect_equal(read.table(system.file("extdata",'test_bed12.txt',package="DRUID"),
                          header = T, stringsAsFactors = FALSE, sep = '\t',row.names = NULL),
               DRUID(ref_gpf_path = system.file("extdata",'transcript2gene.txt',package="DRUID"),
                                 ref_path = system.file("extdata",'annotation.bed',package="DRUID"),
                                 bed_path = system.file("extdata",'bed12_coords.bed',package="DRUID"),
                                 bed6=FALSE))
})
test_that("bed6 works", {
  expect_equal(read.table(system.file("extdata",'test_bed6.txt',package="DRUID"),
                          header = T, stringsAsFactors = FALSE, sep = '\t',row.names = NULL),
               DRUID(ref_gpf_path = system.file("extdata",'transcript2gene.txt',package="DRUID"),
                                ref_path = system.file("extdata",'annotation.bed',package="DRUID"),
                                bed_path = system.file("extdata",'bed6_coords.bed',package="DRUID"),
                                bed6=TRUE))
})
