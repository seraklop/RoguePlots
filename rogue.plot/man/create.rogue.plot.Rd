% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rogue.plot.R
\name{create.rogue.plot}
\alias{create.rogue.plot}
\title{Create Rogue Plots}
\usage{
create.rogue.plot(tree, reftrees, rogues, outgroup = 1,
  type = "greyscale", col = NULL,
  outfile.table = "Rogue_taxon_placements.txt",
  outfile.plot = "Rogue_taxon_placements.pdf", min.prob = 0.01,
  cex.tips = par("cex"), tip.color = "black")
}
\arguments{
\item{tree}{an object of class \code{phylo} (defined in the \code{ape}
package).}

\item{reftrees}{an object of class \code{multiphylo} containing a set of
trees, e.g., from a bootstrap or Bayesian analysis; they all have to
contain the same tips, and the same tips as \code{tree}, potentially plus
the rogue taxa to analyze.}

\item{rogues}{a vector of mode character containing the names of the tip or
tips whose placement is to be calculated.}

\item{outgroup}{a vector of mode numeric or character specifying the new
outgroup.}

\item{type}{a character string specifying the type of formatting to use for
the rogue plot. Options are "greyscale", "user", or any of the 11-colour
palettes from the \code{RColorBrewer} package.}

\item{col}{in case of \code{type = "user"}, a vector of mode character
giving the colours used to reflect the frequencies of attachment of the
rogue taxa to different edges in the tree. Has to contain at least
eleven colours.}

\item{outfile.table}{string specifying the name of the file to which the
table with the taxon placements is saved. Defaults to
"Rogue_taxon_placements.txt".}

\item{outfile.plot}{name of the PDF file to which the plots are saved.
An empty quote (\code{""}) causes the plots to be directed to the standard
plotting device. Defaults to "Rogue_taxon_placements.pdf".}

\item{min.prob}{minimum probability of attachment, below which an edge is
considered as zero for plotting. Defaults to 1\% of reference trees.}

\item{cex.tips}{a numeric value giving the scaling factor of the tip labels.
The default is to take the current value from the graphical parameters.}

\item{tip.color}{the colours used for the tip labels, recycled if shorter
then the number of tips (defaults to "black").}
}
\value{
\code{create.rogue.plot} returns the number of rogue taxa for which
  the plots were generated.
}
\description{
This function can be used to illustrate the placement of rogue taxa in a
  partially resolved (e.g., consensus) tree. The rogue plot has been
  introduced in Klopfstein & Spasojevic (2018) to illustrate the placement of
  fossil taxa on a consensus tree resulting from a phylogenetic analysis of
  a morphological dataset, but can be used to visualize taxon placements in
  various contexts.
}
\details{
\code{create.rogue.plot} is a wrapper function for the other functions in the
  package. It creates rogue plots given a tree to plot, a set of reference
  trees, such as bootstrap or Bayesian posterior trees, and a (list of)
  taxa of interest, such as rogue taxa or fossils.

A rectangular phylogram is used to illustrate the placement of rogue taxa on
  a fully or partially resolved (e.g., consensus) tree, where horizontal
  edges represent the probability of placement on the actual edge, while
  vertical edges are coloured according to the sum of the probabilities that
  the rogue attaches to a branch not present in the tree, but of which the
  current node is the most recent common ancestor.

To arrange the tree that is to be plotted before passing it to this
  function, \code{\link[ape]{root}} and \code{\link[ape]{ladderize}} can be
  used. The reference trees all have to contain the same taxa, and the tree
  to plot needs to have the same taxa or the same taxa minus those passed to
  the function as \code{rogues}. If \code{rogues} is not specified, then
  they are inferred from the difference in the list of tips of the tree
  versus the reference trees.

The value of \code{type} determines the formatting of the rogue plot.
  "greyscale" represents edges with zero to \code{min.prob}
  probability of rogue attachment by black, dotted lines, while solid lines
  in grey scale are used to paint edges with above \code{min.prob}
  probabilities. If a colour version is preferred, one can either choose the
  \code{type = "user"} option and pass a colour vector of at least length 11
  via the \code{col} argument, or use one of the pallet names in the package
  \code{RColorBrewer} with the required length. Examples are "Spectral" or
  "RdYlBu" - see the function \code{display.brewer.all} in the
  \code{RColorBrewer} package.
}
\examples{
## example using random trees
reftrees <- ape::rmtree(5, 21, tip.label = c("rogue.taxon", paste("t", 1:20,
  sep = "")))
reftrees <- c(rep(list(reftrees[[1]]), 60), rep(list(reftrees[[2]]),
  25), rep(list(reftrees[[3]]), 10), reftrees)
tree <- ape::root(reftrees[[1]], outgroup = "t1")
create.rogue.plot(tree, reftrees, rogues = "rogue.taxon", outgroup =
  "t1", type = "Spectral", outfile.table = "Rogue.txt", outfile.plot = "")
create.rogue.plot(tree, reftrees, rogues = "rogue.taxon", outgroup =
  "t1", type = "greyscale", outfile.table = "Rogue.txt", outfile.plot = "")
tree <- ape::root(reftrees[[61]], outgroup = "t1")
create.rogue.plot(tree, reftrees, rogues = "rogue.taxon", outgroup =
  "t1", type = "Spectral", outfile.table = "Rogue.txt", outfile.plot = "")

## example from Klopfstein & Spasojevic 2018
data(consensus.tree)
data(fossil.trees)
data(taxon.lists)
create.rogue.plot(consensus.tree, fossil.trees, rogues = taxon.lists$fossils,
  outgroup = taxon.lists$outgroup, type = "Spectral", col = NULL,
  min.prob = 0.01, cex.tips = 1.2, outfile.table = "Rogue_placements.txt",
  outfile.plot = "Rogue_plots.pdf")
}
\references{
Klopfstein, S., Spasojevic, T. (2018). Illustrating phylogenetic placement
of fossils using RoguePlots: An example from ichneumonid parasitoid wasps
(Hymenoptera, Ichneumonidae) and an extensive morphological matrix. bioRxiv
https://doi.org/10.1101/425090
}
