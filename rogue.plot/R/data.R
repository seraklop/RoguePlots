#' Bayesian Trees of Fossil and Extant Parasitoids
#'
#' A set of 100 phylogenetic trees from a Bayesian analysis of a morphological
#'   dataset including extant and fossil parasitoid wasps. See Klopfstein &
#'   Spasojevic 2018 for details.
#'
#' @format Phylogenetic trees in the \code{multiphylo} format defined in the
#'   \code{ape} package.
#' @references
#' Klopfstein & Spasojevic (2018)
#-------------------------------------------------------------------------------
"fossil.trees"
#-------------------------------------------------------------------------------
#' Consensus Tree of Extant Parasitoids
#'
#' A majority-rule consensus tree of a Bayesian analysis of a morphological
#'   dataset including extant and fossil parasitoid wasps, with the fossil
#'   terminals dropped. See Klopfstein & Spasojevic 2018 for details.
#'
#' @format Phylogenetic tree in the \code{phylo} format defined in the
#'   \code{ape} package.
#' @references
#' Klopfstein & Spasojevic (2018)
#-------------------------------------------------------------------------------
"consensus.tree"
#-------------------------------------------------------------------------------
#' A List of Taxon Groups in the Parasitoids Dataset
#'
#' List of two vectors each containing important taxon sets in the parasitoid
#'   data. See Klopfstein & Spasojevic 2018 for details.
#'
#' @format List containing two vectors of mode \code{character}.
#' \describe{
#'   \item{outgroup}{outgroup taxa to root the trees in the parasitoid dataset}
#'   \item{fossils}{fossil taxa in the dataset}
#'   }
#'
#' @references
#' Klopfstein & Spasojevic (2018)
#-------------------------------------------------------------------------------
"taxon.lists"
#-------------------------------------------------------------------------------
#' The Original Dataset from Klopfstein & Spasejovic (2018)
#'
#' Morphological matrix in its original form as it was used in the parasitoid
#'   fossil study. See Klopfstein & Spasojevic 2018 for details.
#'
#' @format List with as elements for each taxon a vector of mode
#'   \code{character} containing the coding of 222 morphological characters.
#'   The names of the list items correspond to the taxa in the analysis.
#'
#' @references
#' Klopfstein & Spasojevic (2018)
#-------------------------------------------------------------------------------
"orig.morpho"
#-------------------------------------------------------------------------------
#' Expanded Number of States
#'
#' Vector containing for each character in the morphological matrix the number
#'   of states observed in a matrix with more taxa. Used to inform the total
#'   number of states in \code{expand.state.space}.
#'
#' @format vector with as many integers as morphological characters.
#'
#' @references
#' Klopfstein & Spasojevic (2018)
#-------------------------------------------------------------------------------
"tot.nb.states"
