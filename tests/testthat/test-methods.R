test_that("lp assignment method works", {

  load("data.rda")

  lp_adding_results =
    expect_type({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data, 2),
          return_details = TRUE
        )
      })
    }, "list")

  lp_adding_results_plus_loss =
    expect_type({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data, 2),
          return_details = TRUE,
          max_loss_frac = 0.5
        )
      })
    }, "list")

  lp_adding_results_no_optim =
    expect_type({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data, 2),
          return_details = TRUE,
          optimize_values = FALSE
        )
      })
    }, "list")

  lp_adding_results_no_details =
    expect_s3_class({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data, 2),
          return_details = FALSE
        )
      })
    }, "phylo")

  lp_adding_results_bootstrap =
    expect_type({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data, 2),
          return_details = TRUE,
          n_bootstraps = 10
        )
      })
    }, "list")


  # alternative column names and various alternative inputs
  lp_data_alt = lp_data
  colnames(lp_data_alt[[1]]) =
    c("clone", "alt", "dp", "id", "ccf", "vaf", "cn_total")
  colnames(lp_data_alt[[2]]) =
    c("clone", "alt", "dp", "id", "ccf", "vaf", "copy_number")

  # a sample that should be skipped
  lp_data_alt[[3]] =
    lp_data_alt[[3]] %>%
    head(0)

  # a sample that should be renamed
  purity_vals[tree$tip.label[1]] = purity_vals[names(lp_data_alt)[1]]
  names(lp_data_alt)[1] = tree$tip.label[1]

  # various names for mm estimates
  lp_data_alt[[1]]$mm = 1
  lp_data_alt[[2]]$cn_mutated = 1

  lp_adding_results_alt_names =
    expect_type({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data_alt, 3),
          return_details = TRUE,
          max_vaf_bkgr = 0.05,
          purity_estimates = c(
            EPICC_C518_A1_G4_D1 = 0.5,
            EPICC_C518_A1_G10_L1 = 0.5,
            EPICC_C518_A1_G2_L1 = 0.5
          )
        )
      })
    }, "list")

  # test that duplicated sample was renamed:
  expect_match(
      lp_adding_results_alt_names$tree$tip.label,
      paste0(names(lp_data_alt)[1], "_"),
      all = FALSE,
      fixed = TRUE
    )

  # test that sample with no data was skipped:
  expect_no_match(
    lp_adding_results_alt_names$tree$tip.label,
    names(lp_data_alt)[3],
    all = TRUE,
    fixed = TRUE
  )


  # warning when using very high background rate?
  expect_warning({
    expect_output({
      add_lowpass_sampled(
        tree = tree,
        phydata = phy_data,
        sample_data = head(lp_data, 1),
        return_details = TRUE,
        max_vaf_bkgr = 0.5
      )
    })
  }, "Do you really want to set a background rate of")


  # failing when not including essential data?
  for (drop in c("alt", "depth", "cn")) {
    lp_data_alt = lp_data

    for (i in seq_along(lp_data_alt))
      lp_data_alt[[i]][, drop] = NULL

    expect_error({
      expect_output({
        add_lowpass_sampled(
          tree = tree,
          phydata = phy_data,
          sample_data = head(lp_data_alt, 2)
        )
      })
    }, "missing")

  }


  # test all plotting methods with all result sets
  generated_result_sets = list(
    lp_adding_results,
    lp_adding_results_plus_loss,
    lp_adding_results_no_optim,
    lp_adding_results_no_details,
    lp_adding_results_alt_names,
    lp_adding_results_bootstrap
  )

  plotting_functions = list(
    plot_lp_loglik = plot_lp_loglik,
    plot_lp_loglik_edge = plot_lp_loglik_edge,
    plot_sample_data = plot_sample_data,
    plot_tree = plot_tree,
    plot_lp_position_edge = plot_lp_position_edge
  )

  labeller_functions =
    list(
      NULL,
      function(x) gsub("EPICC_", "", x)
    )

  subsets =
    list(
      NULL,
      names(lp_data)[2]
    )

  expect_fail =
    matrix(
      FALSE,
      length(generated_result_sets),
      length(plotting_functions)
    )

  expect_fail[4, -4] = TRUE
  expect_fail[3, 3] = TRUE

  for (i in seq_along(generated_result_sets)) {
    for (j in seq_along(plotting_functions)) {

      d = generated_result_sets[[i]]
      f = plotting_functions[[j]]

      if (expect_fail[i, j]) {
        if (!inherits(d, "phylo"))
          expect_null(f(d))
      } else {

        for (k in seq_along(labeller_functions)) {
          for (l in seq_along(subsets)) {
            expect_s3_class(
              suppressWarnings(
                f(d,
                  labeller_function = labeller_functions[[k]],
                  subset = subsets[[l]]
                )),
              "ggplot"
            )
          }
        }

        if (names(plotting_functions)[j] == "plot_sample_data") {
          expect_s3_class(
            plot_sample_data(
              d,
              labeller_function = function(x) gsub("EPICC_", "", x),
              purity_vals[gsub("_+$", "", d$per_sample_results$sample)]
            ), "ggplot")
        }
      }
    }
  }

})

test_that("lp assignment is consistent", {

  load("data.rda")

  tree_with_lp_added_expected =
    readRDS("expected_example_results.rds") %>%
    (function(x) {x$inital_values$n_cores = 1; return(x) })

  set.seed(123)

  tree_with_lp_added_obs =
    with(example_lp_data,  {
      MLLPT::add_lowpass_sampled(
        tree = tree,
        phydata = phydata,
        sample_data = samples,
        return_details = TRUE,
        n_bootstraps = 10,
        purity_estimates = c(1,1)
      )
    })


  # compare the inputs used
  expect_equal(
    tree_with_lp_added_obs$inital_values,
    tree_with_lp_added_expected$inital_values
  )

  # compare the mutation data used
  expect_equal(
    tree_with_lp_added_obs$mutations_per_edge,
    tree_with_lp_added_expected$mutations_per_edge
  )

  expect_equal(
    tree_with_lp_added_obs$mutation_data,
    tree_with_lp_added_expected$mutation_data %>%
      lapply(lapply, dplyr::filter, N > 0)
  )

  # compare the lik data
  expect_equal(
    tree_with_lp_added_obs$ll_per_edge,
    tree_with_lp_added_expected$ll_per_edge
  )

  expect_equal(
    tree_with_lp_added_obs$max_ll_per_edge,
    tree_with_lp_added_expected$max_ll_per_edge
  )

  # compare the lik data
  expect_equal(
    tree_with_lp_added_obs$ll_per_edge,
    tree_with_lp_added_expected$ll_per_edge
  )

  expect_equal(
    tree_with_lp_added_obs$max_ll_per_edge,
    tree_with_lp_added_expected$max_ll_per_edge
  )


  # compare the results
  expect_equal(
    ape::dist.topo(
      ape::unroot(tree_with_lp_added_obs$tree),
      ape::unroot(tree_with_lp_added_expected$tree)
    ) %>% as.numeric(), 0
  )

  expect_equal(
    tree_with_lp_added_obs$per_sample_results %>% dplyr::select(-idx),
    tree_with_lp_added_expected$per_sample_results
  )

})


test_that("mutation count data files can be read", {

  # objects for cn annotation:
  cn_df = data.frame(
    id = c("chr1:1_A/T", "chr1:2_C/G", "chr1:3_G/C", "chr1:4_T/A"),
    cn = 1:4,
    mm = 1:4
  )

  cn_gr = GenomicRanges::GRanges("chr1", IRanges::IRanges(1:4, 1:4))
  S4Vectors::elementMetadata(cn_gr) = data.frame(cn = 1:4, mm = 1:4)

  # try loading with default
  tab_file = system.file("extdata", "test_mutation.table", package = "MLLPT")
  res = expect_warning(load_genotyping_file(tab_file), ".*CN data missing.*")
  expect_known_hash(res, "437c8be539")

  # modify result so it is equal to the expected using the above annotations
  res$copy_number = 1:4
  res$mm = 1:4

  # 1) Test using a Granges object as input
  expect_equal(
    load_genotyping_file(tab_file, cn_data = cn_gr),
    res
  )

  expect_warning(
    load_genotyping_file(
      file = tab_file,
      cn_data = {
        cn_gr$mm = NULL
        cn_gr
      }
    ),
    "mutation multiplicity data missing! Assuming mm = 1."
  )


  # 2) Test using a data frame as input
  expect_equal(
    load_genotyping_file(tab_file, cn_data = cn_df),
    res
  )

  expect_equal(
    load_genotyping_file(tab_file, cn_data = head(cn_df, 3)),
    dplyr::filter(res, pos %in% 1:3)
  )

  expect_warning(
    load_genotyping_file(tab_file, cn_data = dplyr::select(cn_df, -mm)),
    "mutation multiplicity data missing! Assuming mm = 1."
  ) %>% expect_equal(dplyr::mutate(res, mm = 1))


  # 3) Test using a bad object as input
  expect_error(
    load_genotyping_file(tab_file, cn_data = list(a = 1)),
    "'cn_data' has to be a data[.]frame or GenomicRanges object[.]"
  )

  tmp_file = tempfile()
  readr::write_csv(data.frame(a = 1:5), tmp_file)
  expect_warning(expect_error(
    load_genotyping_file(tmp_file, cn_data = cn_df),
    "Wrong file format[?]"
  ))
  on.exit(file.remove(tmp_file))

})

test_that("additonal tests of plotting functions", {

  load("data.rda")
  tree$tip.label[1] = "EPICC_C518_A1_G4_L1 (Added)"

  cn_data_el =
    data.frame(
      chromosome = c("chr1", "chr2"),
      start = 1,
      end = 1000,
      cn = 2
    )

  cn_data = rep(list(cn_data_el), length(tree$tip.label))
  names(cn_data) = gsub(" [(]Added.*", "", tree$tip.label)

  color_by =
    magrittr::set_names(
      substr(tree$tip.label, 12, 12),
      tree$tip.label
    )

  plot_tree(tree, HI = 0.5) %>% expect_s3_class("ggtree")
  plot_tree(tree, color_by = color_by) %>% expect_s3_class("ggtree")
  plot_tree(tree, color_by = factor(color_by)) %>% expect_s3_class("ggtree")
  MLLPT:::remove_root_tip(tree, "GL") %>% plot_tree() %>% expect_s3_class("ggtree")
  MLLPT:::set_lp_tiplength(tree, 0) %>% plot_tree() %>% expect_s3_class("ggtree")
  plot_tree_and_cn_data(tree, cn_data) %>% expect_s3_class("ggplot")

})
