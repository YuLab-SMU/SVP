data(hpda_spe_cell_dec)


test_that("runENCODE test",{
  res <- hpda_spe_cell_dec |> runENCODE(group.by='cluster_domain')
  expect_true(inherits(res, "SVPExperiment"))

  res2 <- res |> 
      gsvaExp() |> 
      assay() |> 
      as.matrix() |> 
      t() |>
      apply(2, function(x)names(x)[x>0])

  keep.nm <- unique(hpda_spe_cell_dec$cluster_domain)

  res3 <- lapply(keep.nm, 
                 function(x)colnames(res)[res$cluster_domain == x]) |>
      setNames(keep.nm)
  
  expect_equal(res2, res3)
})


