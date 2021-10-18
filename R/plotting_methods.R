#' Plot showing log-likelihood of all possible edges.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#'
plot_lp_loglik =
  function(d, labeller_function=NULL, subset=NULL) {


    if ("min_confidence" %in% names(d$inital_value)) {
      min_confidence = d$inital_values$min_confidence
    } else {
      min_confidence = 1
    }

    d_edge = d$max_ll_per_edge

    if (!is.null(subset)) {
      d_edge = d_edge[names(d_edge) %in% subset]
    }

    d_plot =
      lapply(d_edge, function(x) {p = exp(x-max(x)); p / sum(p)}) %>%
      do.call(what=cbind) %>%
      reshape2::melt() %>%
      magrittr::set_colnames(c("edge","sample","p")) %>%
      dplyr::mutate(is_max = c(p == tapply(p, sample, max)[as.character(sample)])) %>%
      dplyr::mutate(is_accept = p > min_confidence) %>%
      dplyr::mutate(edge=factor(edge, unique(sort(edge)), ordered = TRUE)) %>%
      dplyr::mutate(shape=ifelse(is_accept & is_max, 8, ifelse(is_max, 20, NA)))

    if (!is.null(labeller_function)) {
      d_plot$sample = labeller_function(as.character(d_plot$sample))
    }

    plot_p_edge =
      d_plot %>%
      ggplot(aes(x=sample, y=edge, fill=p+1e-16)) +
      cowplot::theme_cowplot() +
      geom_tile() +
      geom_point(aes(shape=shape), size=1.2) +
      #geom_text(aes(label=ifelse(value < 1e-16, "", sprintf("%1.0e", value)))) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
      xlab("") +
      ylab("Tree edge") +
      labs(fill="p(Edge)") +
      scale_fill_continuous(trans="log10", low="gray98", high="red") +
      scale_shape_identity()

    plot_p_edge
  }


#' Plot showing log-likelihood of fraction mutated from most likely edge.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#'
plot_lp_loglik_edge =
  function(d, labeller_function=NULL, subset=NULL) {

    per_sample_res = d$per_sample_results

    if (!is.null(subset)) {
      per_sample_res = per_sample_res %>% filter(sample %in% subset)
    }

    if (!is.null(labeller_function)) {
      per_sample_res$sample_mod = labeller_function(per_sample_res$sample)
    } else {
      per_sample_res$sample_mod = per_sample_res$sample
    }

    # extract ll data of best edge
    lld = lapply(seq_len(NROW(per_sample_res)), function(i) {
      smp = per_sample_res$sample[i]
      smp_m = per_sample_res$sample_mod[i]
      edge =  per_sample_res$edge[i]
      lapply(names(d$ll_per_edge[[smp]]), function(edge_c) {
        data.frame(
          sample_mod = smp_m,
          edge = edge_c,
          best = (edge_c == edge),
          frac = seq_along(d$ll_per_edge[[smp]][[edge_c]])/length(d$ll_per_edge[[smp]][[edge_c]]),
          llik = d$ll_per_edge[[smp]][[edge_c]]
        )
      }) %>% do.call(what=rbind)
    }) %>% do.call(what=rbind)

    #names(lld) = with(per_sample_res, paste0("Edge ", edge, " - ", sample_mod))

    #
    edge_data_per_smp =
      lld %>%
      split(.$sample) %>%
      lapply(function(x) { dplyr::mutate(x, p = exp(llik - max(llik))) }) %>%
      lapply(function(x) { dplyr::mutate(x, p = p / sum(p)) }) %>%
      do.call(what=rbind)

    per_case_max_pos =
      per_sample_res %>%
      dplyr::mutate(sample_label = paste0("Edge ", edge, " - ", sample_mod)) %>%
      dplyr::select(frac=pi, sample_mod=sample_label)

    anc_list = phangorn::Ancestors(d$inital_values$tree, type = "all")
    for (i in seq_along(anc_list)) anc_list[[i]] = c(rev(anc_list[[i]]), i)
    node_order = unique(unlist(anc_list))
    edge_order = order(match(d$inital_values$tree$edge[,2], node_order), decreasing = FALSE)

    edge_label =
      paste0(
        seq_along(d$inital_values$tree$edge[,1]),
        " (", d$inital_values$tree$edge[,1],
        "-> ", d$inital_values$tree$edge[,2], ")"
      ) %>% magrittr::set_names(seq_along(.))

    edge_label = edge_label[edge_order]

    per_sample_res =
      per_sample_res %>%
      dplyr::mutate(edge=factor(edge, names(edge_label), edge_label, ordered=TRUE))

    d_plot =
      edge_data_per_smp %>%
      dplyr::mutate(edge=factor(edge, names(edge_label), edge_label, ordered=TRUE)) %>%
      ggplot(aes(x=frac, y=p+1e-16, color=best)) +
      theme_cowplot() +
      geom_line() +
      facet_grid(sample_mod~edge, scales = "free_y") +
      xlab("Position on edge") +
      ylab("Likelihood") +
      geom_vline(data=per_sample_res, aes(xintercept=pi), linetype=3) +
      scale_color_manual(values = c("black","red"), breaks = c("FALSE","TRUE")) +
      guides(color="none") +
      scale_x_continuous(n.breaks = 2)

    d_plot
  }


#' Plot showing ML estimate of sample parameters.
#'
#' @param d Object returned from 'add_lowpass_sampled'.
#' @param labeller_function (optional) function returning a modified version of sample labels.
#' @param external_purity_estimate (optional) numeric vector of external purity estimates of each sample.
#' @param label_external_purity (optional) label to be used in the legend of the external purity estimates.
#' @param subset (optional) character vector defining a subset of samples to plot.
#'
#' @return A ggplot object
#' @export
#'
plot_sample_data =
  function(d, labeller_function=function(x) return(x), external_purity_estimate=NULL, label_external_purity="Independent estimate", subset=NULL) {

    if (!d$inital_values$optimize_values) {
      return(NULL)
    }


    per_sample_results = d$per_sample_results

    if (!is.null(external_purity_estimate)) {
      stopifnot(length(external_purity_estimate) == NROW(per_sample_results))
      per_sample_results$external_estimate = external_purity_estimate
    }

    per_sample_results =
      per_sample_results %>%
      cbind(., initial_purity = d$inital_values$purity) %>%
      dplyr::filter(sample %in% subset | is.null(subset)) %>%
      dplyr::mutate(sample=labeller_function(sample)) %>%
      dplyr::mutate(sample=factor(sample, rev(sort(unique(sample))), ordered = TRUE))

    bkgr_vaf_data = data.frame(x=d$inital_values$vaf_bkgr)
    limit_y = min(c(per_sample_results$background_vaf, bkgr_vaf_data$x), na.rm=TRUE) / 10

    plot_vaf =
      per_sample_results %>%
      ggplot(aes(x=sample, y=background_vaf+1e-6)) +
      cowplot::theme_cowplot() +
      geom_point() +
      geom_hline(data=bkgr_vaf_data, aes(yintercept = x, linetype="Initial value")) +
      coord_flip() +
      scale_y_log10(limits=c(limit_y+1e-6, d$inital_values$max_vaf_bkgr+1e-6), n.breaks = 3) +
      scale_linetype_manual(values=4) +
      theme(legend.position = "bottom") +
      labs(linetype="") +
      ylab("Error rate") +
      xlab("") +
      cowplot::background_grid()



    plot_purity =
      per_sample_results %>%
      ggplot(aes(x=sample, y=purity)) +
      theme_cowplot() +
      geom_segment(aes(xend=sample, yend=initial_purity), linetype=3) +
      geom_point(size=2) +
      geom_point(aes(shape="Initial value", color="Initial value", y=initial_purity), size=2) +
      scale_y_continuous(n.breaks = 4, limits = c(0,1)) +
      coord_flip() +
      labs(color="", shape="") +
      ylab("Purity") +
      xlab("") +
      theme(legend.position = "bottom") +
      scale_shape_manual(values=c(4,3), breaks=c("Initial value", label_external_purity)) +
      scale_color_manual(values=c("red","red"), breaks=c("Initial value", label_external_purity)) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
      cowplot::background_grid()

    if (!is.null(external_purity_estimate)) {
      plot_purity = plot_purity +
        geom_segment(aes(xend=sample, yend=external_estimate), linetype=3) +
        geom_point(aes(shape=label_external_purity, color=label_external_purity, y=external_estimate), size=2)
    }


    plot_loss_frac =
      per_sample_results %>%
      ggplot(aes(x=sample, y=loss_frac)) +
      theme_cowplot() +
      geom_point(size=2) +
      coord_flip() +
      labs(color="", shape="") +
      ylab("Loss rate") +
      xlab("") +
      theme(legend.position = "bottom") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
      cowplot::background_grid() +
      scale_y_continuous(n.breaks = 3, limits = c(0,NA))


    plot_dp =
      d$mutation_data %>%
      lapply(do.call, what=rbind) %>%
      lapply(function(x) weighted.mean(x$N, x$weight)) %>%
      reshape2::melt() %>%
      dplyr::mutate(L1=labeller_function(L1)) %>%
      dplyr::mutate(sample=factor(L1, levels(per_sample_results$sample), ordered=TRUE)) %>%
      dplyr::filter(!is.na(sample)) %>%
      ggplot(aes(x=sample, y=value)) +
      theme_cowplot() +
      geom_point() +
      ylim(0, NA) +
      coord_flip() +
      ylab("Coverage") +
      xlab("") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      cowplot::background_grid() +
      scale_y_continuous(n.breaks = 3, limits = c(0,NA))


    if (d$inital_values$max_loss_frac > 0) {
      pg =
        cowplot::plot_grid(
          cowplot::as_grob(plot_vaf),
          cowplot::as_grob(plot_loss_frac),
          cowplot::as_grob(plot_purity),
          cowplot::as_grob(plot_dp),
          align = "h",
          axis = "b",
          rel_widths = c(0.75, 0.45, 0.45, 0.45),
          nrow = 1
        )
    } else {
      pg =
        cowplot::plot_grid(
          cowplot::as_grob(plot_vaf),
          cowplot::as_grob(plot_purity),
          cowplot::as_grob(plot_dp),
          align = "h",
          axis = "b",
          rel_widths = c(0.8, 0.5, 0.5),
          nrow = 1
        )
    }

    return(pg)
  }


#' A internal function used for the plotting of trees
#'
#' @param tree object of class phylo to plot.
#' @param HI (optional) HI of the tree (unused).
#' @param color_by (optional) lookup table with color labels (NULL).
#' @param linewidth (optional) numeric value defining linewidth of tree edges (1).
#' @param pointsize (optional)  numeric value defining point size of tree tips  (1).
#' @param textsize (optional)  numeric value defining linewidth of tree edges (2.5).
#' @param labeller_function (optional) HI of the tree (unused).
#'
#' @return ggplot object
#' @export
#'
#' @examples plot_tree(ape::rtree(5))
plot_tree = function(tree, HI=NULL, color_by=NULL, linewidth=1, pointsize=1, textsize=2.5, labeller_function=NULL) {

  if ("tree" %in% names(tree)) tree = tree$tree
  checkmate::assertClass(tree, "phylo", null.ok = FALSE)
  checkmate::assertNumber(HI, null.ok = TRUE)
  checkmate::assertCharacter(color_by, names = TRUE, null.ok = TRUE)
  checkmate::assertNumber(linewidth, null.ok = FALSE, na.ok = FALSE, lower = 0, finite = TRUE)
  checkmate::assertNumber(pointsize, null.ok = FALSE, na.ok = FALSE, lower = 0, finite = TRUE)
  checkmate::assertNumber(textsize, null.ok = FALSE, na.ok = FALSE, lower = 0, finite = TRUE)
  checkmate::assertFunction(labeller_function, null.ok = TRUE)

  if (is.null(labeller_function)) {
    labeller_function = function(x) x
  }

  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plot tree
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  tree$tip.label = unlist(tree$tip.label)
  added_node = grepl("Added", tree$tip.label)
  #tree$tip.label = gsub(" [(]Added.*[)]", " ", tree$tip.label)
  shape_lookup = if_else(added_node, "A", "U") %>% magrittr::set_names(tree$tip.label)
  new_label = gsub(" [(]Added.*", " *", labeller_function(tree$tip.label))  %>% magrittr::set_names(tree$tip.label)

  if (is.null(color_by)) {
    color_by = rep("gray10", length(tree$tip.label)) %>% magrittr::set_names(tree$tip.label)
    color_scale = scale_color_identity()
  } else {
    color_scale = scale_color_brewer(palette = "Set1")
  }

  plt =
    tree %>%
    ggtree::ggtree(layout="rectangular", color="gray10", size=linewidth) +
    ggtree::geom_tiplab(aes(label=new_label[label], color=color_by[label]), size=textsize, hjust=-0.25) +
    ggtree::geom_tippoint(aes(shape=shape_lookup[label], color=color_by[label]), size=pointsize) +
    ggtree::geom_treescale(y=-1, x=0, fontsize=textsize, linesize=linewidth, offset = 0.9) +
    scale_shape_manual(breaks=c("U","A"), values=c(U=16, A=8)) +
    labs(color="") +
    guides(shape="none") +
    guides(color = guide_legend(nrow=1)) +#, override.aes=aes(label="A"))) +
    theme(title=element_text(size=8, color="gray20")) +
    theme(legend.position="bottom") +
    theme(panel.background = element_rect(fill="transparent", colour = NA)) +
    theme(plot.background = element_rect(fill="transparent", colour = NA)) +
    color_scale

  plt = plt + xlim(0, 1.6 * max(plt$data$x))

  return(plt)
}


#' Internal function that removes the reference tip from a tree
#'
#' @param tree object of class phylo.
#' @param out character vector with the tip name.
#'
#' @return modified tree
#'
remove_root_tip = function(tree, out) {
  data = tree %>% tidytree::as_tibble()

  whO = which(data$label == out)
  whEP = which(data$parent == data$parent[whO] & data$node != data$node[whO] & data$parent != data$node)

  data$branch.length[whEP] = data$branch.length[whEP] + data$branch.length[whO]

  data = data[-whO,]

  data$parent = match(data$parent, data$node)
  data$node = match(data$node, data$node)

  tidytree::as.phylo(data)
}


#' Scale edge length of added LP samples for plotting
#'
#' @param x A phylo object.
#' @param frac_height The edge length relative to the tree height.
#'
#' @return A phylo object.
#' @export
#'
set_lp_tiplength = function(x, frac_height=0.01) {

  checkmate::assertClass(x, "phylo")
  checkmate::assertNumber(frac_height, lower = 0, upper = 1)

  wh_added = grep("Added", x$tip.label)
  idx_added = x$edge[,2] %in% wh_added
  x$edge.length[idx_added] = max(ape::node.depth.edgelength(x)) * frac_height

  return(x)
}

