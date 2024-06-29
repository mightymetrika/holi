# test_that("rstar_lr works", {
#   iris_bin <- iris
#   iris_bin$setosa <- ifelse(iris_bin$Species == "setosa", 1, 0)
#
#   rs <- rstar_lr(setosa ~ Sepal.Length + Petal.Length, # Sepal.Width + Petal.Width,
#                  .data = iris_bin)
#   data(Seatbelts)
#   rs <-rstar_lr(law ~ DriversKilled + VanKilled + drivers + kms, .data = Seatbelts)
#   expect_equal(2 * 2, 4)
# })

