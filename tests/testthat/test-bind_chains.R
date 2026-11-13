
skip_on_cran()

net <- set_agd_arm(smoking, studyn, trtn, r = r, n = n)

fit1 <- suppressWarnings(nma(net, iter = 50))
fit2 <- suppressWarnings(nma(net, iter = 50))
fit3 <- suppressWarnings(nma(net, trt_effects = "random", iter = 50))
fit4 <- suppressWarnings(nma(net, iter = 10))
fit5 <- suppressWarnings(nma(net, iter = 50, warmup = 40))
fit6 <- suppressWarnings(nma(net, iter = 50, chains = 2))

test_that("bind_chains error checks", {
  expect_error(bind_chains(fit1, "a"), "Can only combine stan_nma objects")
  expect_error(bind_chains(fit1, net), "Can only combine stan_nma objects")
  expect_error(bind_chains(fit1, fit3), "Models are not identical")
  expect_error(bind_chains(fit1, fit4), "Number of iterations not equal")
  expect_error(bind_chains(fit1, fit5), "Number of warmup iterations not equal")
})

test_that("bind_chains output", {
  bind1 <- bind_chains(fit1)
  bind12 <- bind_chains(fit1, fit2)
  bind126 <- bind_chains(fit1, fit2, fit6)

  expect_s3_class(bind12, "stan_nma")

  expect_equal(bind1$stanfit@sim$chains, 4)
  expect_equal(bind12$stanfit@sim$chains, 8)
  expect_equal(bind126$stanfit@sim$chains, 10)

  expect_equal(dim(as.array(bind1, pars = "d")), c(25, 4, 3))
  expect_equal(dim(as.array(bind12, pars = "d")), c(25, 8, 3))
  expect_equal(dim(as.array(bind126, pars = "d")), c(25, 10, 3))

  expect_equal(as.array(bind1), as.array(fit1))
  expect_equal(as.array(bind12)[,1:4,], unclass(as.array(fit1)))
  expect_equal(as.array(bind12)[,5:8,], unclass(as.array(fit2)))
  expect_equal(as.array(bind126)[,1:4,], unclass(as.array(fit1)))
  expect_equal(as.array(bind126)[,5:8,], unclass(as.array(fit2)))
  expect_equal(as.array(bind126)[,9:10,], unclass(as.array(fit6)))
})

test_that("bind_chains with list input", {
  expect_error(bind_chains(list(fit1), fit2), "Can only combine stan_nma objects")

  bind1 <- bind_chains(list(fit1))
  bind12 <- bind_chains(list(fit1, fit2))
  bind126 <- bind_chains(list(fit1, fit2, fit6))

  expect_s3_class(bind12, "stan_nma")

  expect_equal(bind1$stanfit@sim$chains, 4)
  expect_equal(bind12$stanfit@sim$chains, 8)
  expect_equal(bind126$stanfit@sim$chains, 10)

  expect_equal(dim(as.array(bind1, pars = "d")), c(25, 4, 3))
  expect_equal(dim(as.array(bind12, pars = "d")), c(25, 8, 3))
  expect_equal(dim(as.array(bind126, pars = "d")), c(25, 10, 3))

  expect_equal(as.array(bind1), as.array(fit1))
  expect_equal(as.array(bind12)[,1:4,], unclass(as.array(fit1)))
  expect_equal(as.array(bind12)[,5:8,], unclass(as.array(fit2)))
  expect_equal(as.array(bind126)[,1:4,], unclass(as.array(fit1)))
  expect_equal(as.array(bind126)[,5:8,], unclass(as.array(fit2)))
  expect_equal(as.array(bind126)[,9:10,], unclass(as.array(fit6)))
})

test_that("cbind alias", {
  expect_equal(cbind(fit1, fit2), bind_chains(fit1, fit2))
})


a1 <- as.array(fit1)
a2 <- as.array(fit2)
a3 <- as.array(fit3)
a4 <- as.array(fit4)
a5 <- as.array(fit5)
a6 <- as.array(fit6)

test_that("cbind.mcmc_array error checking", {
  expect_error(cbind(a1, "a"), "Can only combine mcmc_array objects")
  expect_error(cbind(a1, net), "Can only combine mcmc_array objects")
  expect_error(cbind(a1, a3), "Number of parameters not equal")
  expect_error(cbind(a1, a4), "Number of iterations not equal")
  expect_error(cbind(a1, a5), "Number of iterations not equal")
  expect_error(cbind(as.array(fit1, pars = c("d[1]", "d[2]")),
                     as.array(fit1, pars = c("d[2]", "d[4]"))),
               "Parameter names do not match")
})

test_that("cbind.mcmc_array output", {
  bind1 <- cbind(a1)
  bind12 <- cbind(a1, a2)
  bind126 <- cbind(a1, a2, a6)

  expect_s3_class(bind12, "mcmc_array")

  expect_equal(dim(bind1), c(25, 4, 178))
  expect_equal(dim(bind12), c(25, 8, 178))
  expect_equal(dim(bind126), c(25, 10, 178))

  expect_equal(bind1, a1)
  expect_equal(bind12[,1:4,], unclass(a1))
  expect_equal(bind12[,5:8,], unclass(a2))
  expect_equal(bind126[,1:4,], unclass(a1))
  expect_equal(bind126[,5:8,], unclass(a2))
  expect_equal(bind126[,9:10,], unclass(a6))
})
