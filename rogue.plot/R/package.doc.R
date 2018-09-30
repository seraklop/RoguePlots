#' rogue.plot: A package to create visualizations of rogue taxon placement.
#'
#' This package contains functions for plotting rogue graphs as introduced in
#'   Klopfstein & Spasojevic (2018). The \code{\link{create.rogue.plot}} function offers
#'   a convenient wrapper around the functions \code{\link{get.rogue.placement}},
#'   \code{\link{get.edge.formats}}, and \code{\link{phylogramm.plot.twocoloured}}, which can
#'   also be used individually to create rogue plots.
#'   Additionally, the package allows manipulations of morphological matrices: it reads in
#'   nexus-formated standard data that can include polymorphisms with
#'   (\code{\link{read.nexus.morpho}}), restrict the state space of the characters to those
#'   states observed in the data (\code{\link{restrict.state.space}}), or expand the state
#'   space to account for unobserved states (\code{\link{expand.state.space}}). Furthermore,
#'   for instance for use with an asymmetric transition model as the one implemented in
#'   MrBayes, the state labes can be randomized (\code{\link{randomize.state.labels}}).
#'   All these functions take care of polymorphisms properly.
#'
#' @docType package
#' @name rogue.plot
#-------------------------------------------------------------------------------
NULL
