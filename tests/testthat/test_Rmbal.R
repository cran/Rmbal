test_that("error if one p and Wp have the same length", {

   expect_equal(
      length(p <- c(3700, 3650, 3400, 3100, 2800, 2500, 2200, 1900, 1600, 1300, 1000, 700,
      600)),
      length(Wp <- rep(0, length.out = (length(p))))
   )

})

