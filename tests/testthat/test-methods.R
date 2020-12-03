test_that("lp assignment method works", {

  load("data.rda")

  lp_adding_results =
    expect_type({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 2),
        return_details = TRUE
      )
    }, "list")

  expect_s3_class(plot_lp_loglik(lp_adding_results), "ggplot")
  expect_s3_class(plot_lp_loglik_edge(lp_adding_results), "ggplot")
  expect_s3_class(plot_sample_data(lp_adding_results), "ggplot")

})
