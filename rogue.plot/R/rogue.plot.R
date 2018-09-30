#
#-------------------------------------------------------------------------------
#' Read NEXUS-formated morphological data with polymorphisms
#'
#' This function allows reading in 'standard' data in NEXUS format that can
#'   include polymorphisms.
#'
#' @param file a file name specified either by a variable of mode character or
#'   by a double-quoted string.
#' @param return.as.list a logical specifying whether the alignment should be
#'   returned as a list (if \code{TRUE}) or as a matrix (if \code{FALSE}).
#' @details
#' Morphological data in NEXUS format is read in from the first data block in
#'   the specified file (all following data blocks are ignored). The data has
#'   to be on a single line, after the taxon identifier followed by spaces or
#'   tabs - the 'interleaved' data format is not supported.
#' Polymorphic character states have to be specified in brackets and separated
#'   by commas, for example \code{(0,1,4)}.
#'
#' @return
#' \code{read.nexus.morpho} returns either a list with each entry representing
#'   a taxon (if \code{return.as.list} equals \code{TRUE}), or a matrix with
#'   each row representing a taxon.
#'
#' @references
#' Klopfstein, S., Spasojevic, T. (2018).
#'
#' @examples
#' ## first writing the dataset from Klopfstein & Spasojevic (2018) to a file
#' ## using the \code{write.nexus.data} function from the package \code{ape},
#' ## then reading it in using the current function which handles polymorphies.
#' data(orig.morpho)
#' ape::write.nexus.data(orig.morpho, file = "morpho.temp", format = "standard",
#'    interleaved = FALSE)
#' morpho <- read.nexus.morpho(file = "morpho.temp", return.as.list = TRUE)
#' str(morpho)
#' unlink("morpho.temp")
#-------------------------------------------------------------------------------
#' @export
read.nexus.morpho <- function(file, return.as.list = TRUE)
{
    lines <- readLines(file)
    if (lines[1] != "#NEXUS"){
        stop("File is not marked as NEXUS file in first line. Probable file format error. Exiting...")
    }
    lines <- lines[which(lines != "")]      #remove empty lines
    lines <- lines[-grep("^\\[", lines)]    #remove lines which start with a comment
    data.start <- grep("begin data", tolower(lines))[1]
    data.end <- grep("end;", tolower(lines))
    data.end <- data.end[which(data.end > data.start)][1]
    lines <- lines[(data.start +1):(data.end -1)]
    dims <- gsub("\t", "", gsub(";", "", tolower(strsplit(lines[grep("dimensions", tolower(lines))[1]], split = " ")[[1]])))
    ntax <- as.numeric(gsub("ntax=", "", dims[grep("ntax=", dims)]))
    nchar <- as.numeric(gsub("nchar=", "", dims[grep("nchar=", dims)]))
    formats <- gsub("\t", "", gsub(";", "", tolower(strsplit(lines[grep("format", tolower(lines))[1]], split = " ")[[1]])))
    datatype <- gsub("datatype=", "", formats[grep("datatype=", formats)])
    if(datatype != "standard"){
        stop("Data format not 'standard'. Exiting...")
    }
    missChar <- gsub("missing=", "", formats[grep("missing=", formats)])
    gapChar <- gsub("gap=", "", formats[grep("gap=", formats)])
    if(length(gapChar) == 0){
        gapChar <- "-"
    }

    #now identify and go through matrix part
    mat.start <- grep("matrix", tolower(lines)) +1
    mat.end <- grep(";", lines)[which(grep(";", lines) > mat.start)[1]]
    lines <- lines[mat.start:mat.end]
    lines[length(lines)] <- gsub(";", "", lines[length(lines)])
    lastline <- gsub("[[:blank:]]", "", lines[length(lines)])
    if (lastline == ""){
        lines <- lines[-length(lines)]
    }
    char.matrix <- matrix(nrow = ntax, ncol = nchar, dimnames = list(1:ntax, NULL))
    firstRound <- T
    for (i in 1:length(lines)){
        this.dat <- strsplit(lines[i], "[[:blank:]]")[[1]]
        this.dat <- this.dat[which(this.dat != "")]
        if (firstRound){
            rownames(char.matrix)[i] <- this.dat[1]
            d <- paste(this.dat[2:length(this.dat)], collapse = "")
            d <- gsub(gapChar, missChar, d, fixed = T)
            d <- strsplit(strsplit(d, "(", fixed = T)[[1]], ")", fixed = T)
            d <- d[which(lapply(d, length) > 0)]
            counter <- 1
            for (j in 1:length(d)){
                if(length(d[[j]]) == 1){
                    if (length(grep(",", d[[j]], fixed = T)) == 1){
                        char.matrix[i, counter] <- paste("(", d[[j]][1], ")", sep = "")
                        counter <- counter +1
                    } else {
                        new.chars <- strsplit(d[[j]], split = "")[[1]]
                        char.matrix[i, counter:(counter + length(new.chars) -1)] <- new.chars
                        counter <- counter + length(new.chars)
                    }
                } else {
                    char.matrix[i, counter] <- paste("(", d[[j]][1], ")", sep = "")
                    counter <- counter +1
                    new.chars <- strsplit(d[[j]][2], split = "")[[1]]
                    if ((counter + length(new.chars) -1) > nchar){
                        write.table(char.matrix, file = "error-matrix.tmp", quote = F, sep = "\t", row.names = T, col.names = F)
                        stop(paste("Error: too many characters encountered for taxon '", this.dat[1], "'. Next chars: '", paste(new.chars, collapse = ""), "'. Printing erroneous matrix to file 'error-matrix.tmp'..."))
                    }
                    char.matrix[i, counter:(counter + length(new.chars) -1)] <- new.chars
                    counter <- counter + length(new.chars)
                }
            }
            if (i == ntax){
                firstRound <- F
            }
        } else {
            curr.taxon <- which(rownames(char.matrix) == this.dat[1])
            if(length(curr.taxon) == 0){
                stop(paste("Error: taxon '", this.dat[1], "' not found in taxon list. Format proper interleaved? Error occurred on the following line: '", lines[i], "'. Exiting..."))
            }
            d <- paste(this.dat[2:length(this.dat)], collapse = "")
            d <- strsplit(strsplit(d, "(", fixed = T)[[1]], ")", fixed = T)
            counter <- which(is.na(char.matrix[curr.taxon, ]))[1]
            for (j in 1:length(d)){
                if(length(d[[j]]) == 1){
                    if (grep(",", d[[j]], fixed = T) == 1){
                        char.matrix[i, counter] <- d[[j]][1]
                        counter <- counter +1
                    } else {
                        new.chars <- strsplit(d[[j]], split = "")[[1]]
                        char.matrix[i, counter:(counter + length(new.chars) -1)] <- new.chars
                        counter <- counter + length(new.chars)
                    }
                } else {
                    char.matrix[i, counter] <- d[[j]][1]
                    counter <- counter +1
                    new.chars <- strsplit(d[[j]][2], split = "")[[1]]
                    char.matrix[i, counter:(counter + length(new.chars) -1)] <- new.chars
                    counter <- counter + length(new.chars)
                }
            }
        }
    }
    if(any(is.na(char.matrix))){
        write.table(char.matrix, file = "error-matrix.tmp", quote = F, sep = "\t", row.names = T, col.names = F)
        stop(paste("Error: Matrix of ", ntax, " taxa and ", nchar, "characters not entirely filled. Printint erroneous matrix to file 'error-matrix.tmp'..."))
    }
    if (return.as.list){
        new.mat <- split(char.matrix, 1:nrow(char.matrix))
        names(new.mat) <- rownames(char.matrix)
        char.matrix <- new.mat
    }
    return(char.matrix)
}
#-------------------------------------------------------------------------------
#' Restrict state space
#'
#' Removing unobserved states from an alignment of standard characters,
#'   preserving polymorphisms.
#'
#' @param alignment An alignment of characters of type "standard", in the
#'   form of either a list of the same length as there are taxa in the
#'   alignment, or as a matrix with as many rows as taxa and as many columns
#'   as characters.
#' @param return.as.list a logical specifying whether the new alignment should
#'   be returned as a list (if \code{TRUE}) or as a matrix (if \code{FALSE}).
#' @param log.name Name of the log file to which a summary of the changes
#'   is written. Defaults to"StateChanges_restrict.log".
#'
#' @details
#' The "standard" data type as it is defined in the NEXUS format is used to
#'   represent discrete characters such as they result from morphological
#'   examination. They are often represented by numbers as state labels ("0",
#'   "1", "2", etc.) and their evolution can be modelled using Lewis' (2001)
#'   likelihood approach.
#' Phylogenetic inference programs such as MrBayes that can handle such
#'   discrete characters, assume that the true number of states that are
#'   possible in each character is given by the largest number observed in
#'   the corresponding column of the data matrix. This function removes empty
#'   states and recodes the remaining ones to ascertain that the smallest state
#'   label is "0" and the largest equals the number of states minus 1. Changes
#'   to state labels are recorded in a log file.
#'
#' @return
#' \code{restrict.state.space} returns either a list with each element
#'   representing a taxon (if \code{return.as.list} equals \code{TRUE}), or a
#'   matrix with each row representing a taxon.
#'
#' @references
#' Lewis (2001) A likelihood approach to estimating phylogeny from discrete
#'   morphological character data. Systematic Biology, 50, 913-925.
#'
#' @examples
#' data(orig.morpho)
#' restr.morpho <- restrict.state.space(orig.morpho, return.as.list = TRUE)
#' str(unlist(lapply(orig.morpho, "[", 222)))
#' str(unlist(lapply(restr.morpho, "[", 222)))
#-------------------------------------------------------------------------------
#' @export
restrict.state.space <- function(alignment, return.as.list = F, log.name = "StateChanges_restrict.log"){

    align.mat <- check.alignment.format(alignment, function.name = "restrict.state.space")
    write(paste("character_#", "removed_state(s)", sep = "\t"), file = log.name, append = F)

    for (i in 1:dim(align.mat)[2]) {
        chars <- paste(unlist(strsplit(align.mat[ ,i], split = ",")), collapse = "")
        chars <- gsub("?", "", gsub(")", "", gsub("(", "", chars, fixed = T), fixed = T), fixed = T)
        chars <- unlist(strsplit(chars, split = ""))
        nb.states <- length(unique(chars))

        if (nb.states < max(as.numeric(unique(chars))) + 1) {
            max.label <- max(as.numeric(unique(chars)))
            s <- 0
            repl <- F
            nb.minus <- 1
            repl.states <- NULL
            while (s < max.label){
                s.label <- as.character(s)
                pos <- grep(s.label, align.mat[, i])

                if (repl){
                    if (length(pos) > 0){
                        for (pp in pos){
                            s.new.label <- as.character( s - nb.minus )
                            align.mat[pp, i] <- gsub(pattern = s.label, replacement = s.new.label, align.mat[pp, i])
                        }
                    } else {
                        nb.minus <- nb.minus + 1
                    }
                }

                if (length(pos) == 0){
                    repl <- T
                    repl.states <- c(repl.states, s.label)
                }
                s <- s +1
            }
            write(paste(i, repl.states, sep = "\t"), file = log.name, append = T)
        }
    }

    if (return.as.list) {
        new.mat <- split(align.mat, 1:nrow(align.mat))
        names(new.mat) <- rownames(align.mat)
        align.mat <- new.mat
    }

    return(align.mat)

}
#-------------------------------------------------------------------------------
#' Expand state space
#'
#' Adding unobserved states to an alignment of standard characters.
#'
#' @param alignment An alignment of characters of type "standard", in the
#'   form of either a list of the same length as there are taxa in the
#'   alignment, or as a matrix with as many rows as taxa and as many columns
#'   as characters.
#' @param nb.states a single value or vector of type \code{numeric} giving the
#'   number of states each characters should contain. These number have to be
#'   equal or higher than the current number of states in each character;
#'   otherwise the matrix is returned unchanged.
#' @param return.as.list a logical specifying whether the new alignment should
#'   be returned as a list (if \code{TRUE}) or as a matrix (if \code{FALSE}).
#' @param log.name name of the log file to which details about the conducted
#'   changes are written.
#'
#' @details
#' The "standard" data type as it is defined in the NEXUS format is used to
#'   represent discrete characters such as they result from morphological
#'   examination. They are often represented by numbers as state labels ("0",
#'   "1", "2", etc.) and their evolution can be modelled using Lewis' (2001)
#'   likelihood approach.
#' Phylogenetic inference programs such as MrBayes that can handle such
#'   discrete characters, assume that the true number of states that are
#'   possible in each character is given by the largest number observed in
#'   the corresponding column of the data matrix. This function adds empty
#'   states to ascertain that the largest state label equals the target
#'   number of states as specified in \code{nb.states}. Changes to state labels
#'   are recorded in a log file.
#' The target number of states can for instance stem from theoretical
#'   considerations or from character states in a more inclusive matrix.
#'   \code{nb.states}) can either contain a single value, which is then assumed
#'   to be the same for all the characters, or a vector of the same length as
#'   the number of characters in the matrix. Only single-digit numbers ("0"-"9")
#'   are currently allowed, a restriction also imposed by software such as
#'   MrBayes.
#'
#' @return
#' \code{expand.state.space} returns either a list with each element
#'   representing a taxon (if \code{return.as.list} equals \code{TRUE}), or a
#'   matrix with each row representing a taxon.
#'
#' @references
#' Lewis (2001) A likelihood approach to estimating phylogeny from discrete
#'   morphological character data. Systematic Biology, 50, 913-925.
#'
#' @examples
#' ## compare the state labels
#' data(orig.morpho)
#' data(tot.nb.states)
#' expanded.morpho <- expand.state.space(orig.morpho, nb.states = tot.nb.states,
#'   return.as.list = TRUE)
#' str(unlist(lapply(orig.morpho, "[", 36)))
#' str(unlist(lapply(expanded.morpho, "[", 36)))
#-------------------------------------------------------------------------------
#' @export
expand.state.space <- function(alignment, nb.states, return.as.list = F,
                               log.name = "StateChanges_expanded.log"){

    align.mat <- check.alignment.format(alignment, function.name =
                                            "expand.state.space")

    if(length(nb.states) < 1){
        stop("Error in 'expand.state.space': no vector with state dimensions
             passed. Exiting...")
    }
    if(length(nb.states) == 1){
        nb.states <- rep(nb.states, dim(align.mat)[2])
    } else if(length(nb.states) != dim(align.mat)[2]){
        stop("Error in 'expand.state.space': vector with state dimensions is not
             of the same length as the number of characters (or 1). Exiting...")
    }
    write(paste("character_#", "new_max", sep = "\t"), file = log.name, append = F)

    for (i in 1:dim(align.mat)[2]) {
        chars <- paste(unlist(strsplit(align.mat[ ,i], split = ",")), collapse = "")
        chars <- gsub("?", "", gsub(")", "", gsub("(", "", chars, fixed = T), fixed = T), fixed = T)
        chars <- unlist(strsplit(chars, split = ""))
        max.label <- max(as.numeric(unique(chars)))

        if (max.label < nb.states[i] - 1) {
            pos <- grep(max.label, align.mat[, i])
            max.new.label <- as.character(nb.states[i] - 1)
            for (pp in pos){
                align.mat[pp, i] <- gsub(pattern = max.label, replacement = max.new.label, align.mat[pp, i])

            }
            write(paste(i, max.new.label, sep = "\t"), file = log.name, append = T)
        }
    }

    if (return.as.list) {
        new.mat <- split(align.mat, 1:nrow(align.mat))
        names(new.mat) <- rownames(align.mat)
        align.mat <- new.mat
    }
    return(align.mat)
    }
#-------------------------------------------------------------------------------
#' Randomize state labels
#'
#' State labels are randomized, perserving the size of the state space,
#'    polymorphisms, and ordered characters.
#'
#' @param alignment An alignment of characters, either standard or DNA, in the
#'   form of either a list of the same length as there are taxa in the
#'   alignment, or as a matrix with as many rows as taxa and as many columns
#'   as characters.
#' @param ordered.chars an optional vector giving the indices of those
#'   characters that should be treated as ordered. Defaults to \code{NULL}, in
#'   which case all characters are assumed as being unordered.
#' @param return.as.list a logical specifying whether the new alignment should
#'   be returned as a list (if \code{TRUE}) or as a matrix (if \code{FALSE}).
#'
#' @details
#' State labels typically have different meanings among different characters in
#'   a morphological matrix. Nevertheless, the state label "0" is more often
#'   used for "absence" of a particular state; state labels are thus arbitrary,
#'   but not entirely random. The asymmetrical model (Lewis 2001) implemented
#'   in software such as MrBayes allows for asymmetry in the transitions
#'   between states, assuming an underlying binomial (for binary characters) or
#'   Dirichlet distribution (for multi-state characters), which implies that
#'   state labels should be random, and not only arbitrary.
#' This function randomizes the state labels of all the characters so that they
#'   are truly random. It can cope with ordered characters by only randomizing
#'   the direction of the ordering, while preserving the sequence. It is also
#'   treating polymorphisms properly.
#'
#' @return
#' \code{randomize.state.labels} returns either a list with each element
#'   representing a taxon (if \code{return.as.list} equals \code{TRUE}), or a
#'   matrix with each row representing a taxon.
#'
#' @references
#' Lewis (2001) A likelihood approach to estimating phylogeny from discrete
#'   morphological character data. Systematic Biology, 50, 913-925.
#'
#' @examples
#' data(orig.morpho)
#' rand.morpho <- randomize.state.labels(orig.morpho, return.as.list = TRUE)
#' str(unlist(lapply(orig.morpho, "[", 14)))
#' str(unlist(lapply(rand.morpho, "[", 14)))
#-------------------------------------------------------------------------------
#' @export
randomize.state.labels <- function(alignment, ordered.chars = NULL, return.as.list = F)
{
    align.mat <- check.alignment.format(alignment, function.name = "randomize.state.labels")

    ntax <- dim(align.mat)[1]
    nchar <- dim(align.mat)[2]

    for ( i in 1:nchar ){
        curr.states <- gsub("(", "", gsub(")", "", align.mat[,i], fixed = T), fixed = T)
        states <- unique(unlist(strsplit(curr.states, ",")))
        states <- states[which(states != "?")]
        if (i %in% ordered.chars) {
            if (rbinom(1, size = 1, prob = 0.5)){
                states <- sort(states)
                newstates <- sort(states, decreasing = T)
            } else {
                newstates <- states
            }
        } else {
            newstates <- sample(states, size = length(states), replace = F)
        }
        for ( j in 1:ntax ){
            prevChar <- strsplit(curr.states[j], split = ",", fixed = T)[[1]]
            if ( prevChar[1] != "?" ){
                if (length(prevChar) == 1){
                    align.mat[j, i] <- newstates[which(states == prevChar)]
                } else {
                    new.char <- newstates[which(states == prevChar[1])]
                    for (c in 2:length(prevChar)){
                        new.char <- c(new.char, newstates[which(states == prevChar[c])])
                    }
                    new.char <- paste("(", paste(new.char, collapse = ","), ")", sep = "")
                    align.mat[j, i] <- new.char
                }
            }
        }
    }
    if (return.as.list){
        new.mat <- split(align.mat, 1:nrow(align.mat))
        names(new.mat) <- rownames(align.mat)
        align.mat <- new.mat
    }
    return(align.mat)
}
#-------------------------------------------------------------------------------
#' Create Rogue Plots
#'
#' This function can be used to illustrate the placement of rogue taxa in a
#'   partially resolved (e.g., consensus) tree. The rogue plot has been
#'   introduced in Klopfstein & Spasejovic 2018 to illustrate the placement of
#'   fossil taxa on a consensus tree resulting from a phylogenetic analysis of
#'   a morphological dataset, but can be used to visualize taxon placements in
#'   various contexts.
#' @param tree an object of class \code{phylo} (defined in the \code{ape}
#'   package).
#' @param reftrees an object of class \code{multiphylo} containing a set of
#'   trees, e.g., from a bootstrap or Bayesian analysis; they all have to
#'   contain the same tips, and the same tips as \code{tree}, potentially plus
#'   the rogue taxa to analyze.
#' @param rogues a vector of mode character containing the names of the tip or
#'   tips whose placement is to be calculated.
#' @param outgroup a vector of mode numeric or character specifying the new
#'   outgroup.
#' @param type a character string specifying the type of formatting to use for
#'   the rogue plot. Options are "greyscale", "user", or any of the 11-colour
#'   palettes from the \code{RColorBrewer} package.
#' @param col in case of \code{type = "user"}, a vector of mode character
#'   giving the colours used to reflect the frequencies of attachment of the
#'   rogue taxa to different edges in the tree. Has to contain at least
#'   eleven colours.
#' @param min.prob minimum probability of attachment, below which an edge is
#'   considered as zero for plotting. Defaults to 1\% of reference trees.
#' @param outfile.table string specifying the name of the file to which the
#'   table with the taxon placements is saved. Defaults to
#'   "Rogue_taxon_placements.txt".
#' @param outfile.plot name of the PDF file to which the plots are saved.
#'   An empty quote (\code{""}) causes the plots to be directed to the standard
#'   plotting device. Defaults to "Rogue_taxon_placements.pdf".
#' @param cex.tips a numeric value giving the scaling factor of the tip labels.
#'   The default is to take the current value from the graphical parameters.
#' @param tip.color the colours used for the tip labels, recycled if shorter
#'   then the number of tips (defaults to "black").
#'
#' @details
#' \code{create.rogue.plot} is a wrapper function for the other functions in the
#'   package. It creates rogue plots given a tree to plot, a set of reference
#'   trees, such as bootstrap or Bayesian posterior trees, and a (list of)
#'   taxa of interest, such as rogue taxa or fossils.
#'
#' A rectangular phylogram is used to illustrate the placement of rogue taxa on
#'   a fully or partially resolved (e.g., consensus) tree, where horizontal
#'   edges represent the probability of placement on the actual edge, while
#'   vertical edges are coloured according to the sum of the probabilities that
#'   the rogue attaches to a branch not present in the tree, but of which the
#'   current node is the most recent common ancestor.
#'
#' To arrange the tree that is to be plotted before passing it to this
#'   function, \code{\link[ape]{root}} and \code{\link[ape]{ladderize}} can be
#'   used. The reference trees all have to contain the same taxa, and the tree
#'   to plot needs to have the same taxa or the same taxa minus those passed to
#'   the function as \code{rogues}. If \code{rogues} is not specified, then
#'   they are inferred from the difference in the list of tips of the tree
#'   versus the reference trees.
#'
#' The value of \code{type} determines the formatting of the rogue plot.
#'   "greyscale" represents edges with zero to \code{min.prob}
#'   probability of rogue attachment by black, dotted lines, while solid lines
#'   in grey scale are used to paint edges with above \code{min.prob}
#'   probabilities. If a colour version is preferred, one can either choose the
#'   \code{type = "user"} option and pass a colour vector of at least length 11
#'   via the \code{col} argument, or use one of the pallet names in the package
#'   \code{RColorBrewer} with the required length. Examples are "Spectral" or
#'   "RdYlBu" - see the function \code{display.brewer.all} in the
#'   \code{RColorBrewer} package.
#'
#' @return \code{create.rogue.plot} returns the number of rogue taxa for which
#'   the plots were generated.
#'
#' @references Klopfstein, S., Spasojevic, T. (2018)
#'
#' @examples
#' ## example using random trees
#' reftrees <- ape::rmtree(5, 21, tip.label = c("rogue.taxon", paste("t", 1:20,
#'   sep = "")))
#' reftrees <- c(rep(list(reftrees[[1]]), 60), rep(list(reftrees[[2]]),
#'   25), rep(list(reftrees[[3]]), 10), reftrees)
#' tree <- ape::root(reftrees[[1]], outgroup = "t1")
#' create.rogue.plot(tree, reftrees, rogues = "rogue.taxon", outgroup =
#'   "t1", type = "Spectral", outfile.table = "Rogue.txt", outfile.plot = "")
#' create.rogue.plot(tree, reftrees, rogues = "rogue.taxon", outgroup =
#'   "t1", type = "greyscale", outfile.table = "Rogue.txt", outfile.plot = "")
#' tree <- ape::root(reftrees[[61]], outgroup = "t1")
#' create.rogue.plot(tree, reftrees, rogues = "rogue.taxon", outgroup =
#'   "t1", type = "Spectral", outfile.table = "Rogue.txt", outfile.plot = "")
#'
#' ## example from Klopfstein & Spasojevic 2018
#' data(consensus.tree)
#' data(fossil.trees)
#' data(taxon.lists)
#' create.rogue.plot(consensus.tree, fossil.trees, rogues = taxon.lists$fossils,
#'   outgroup = taxon.lists$outgroup, type = "Spectral", col = NULL,
#'   min.prob = 0.01, cex.tips = 1.2, outfile.table = "Rogue_placements.txt",
#'   outfile.plot = "Rogue_plots.pdf")
#-------------------------------------------------------------------------------
#' @export
create.rogue.plot <- function (tree, reftrees, rogues, outgroup = 1, type = "greyscale", col = NULL, outfile.table = "Rogue_taxon_placements.txt", outfile.plot = "Rogue_taxon_placements.pdf", min.prob = 0.01, cex.tips = par("cex"), tip.color = "black"){

    if (missing(rogues) ){
        rogues <- setdiff(reftrees[[1]]$tip.label, tree$tip.label)
    }
    if (setequal(tree$tip.label, reftrees[[1]]$tip.label)){
        tree <- ape::drop.tip(tree, rogues)
    }

    taxa.total <- union(tree$tip.label, rogues)
    if (!setequal(taxa.total, reftrees[[1]]$tip.label)){
        stop("Error: taxa in the first reference tree do not match those in the tree to plot (excluding the rogue taxon/taxa). Check tip labels for consistency. Exiting...")
    }

    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        avail.palettes <- rownames(RColorBrewer::brewer.pal.info)[which(RColorBrewer::brewer.pal.info$maxcolors == 11)]
    } else {
        avail.palettes <- NULL
    }


    if (tolower(type) == "greyscale"){
        zero.LWD <- 4
        zero.COL <- "grey85"
        LWD <- 4
        COL <- paste("grey", seq(from = 70, to = 5, by = -6), sep = "")
    } else if (is.element(type, avail.palettes) ){
        zero.LWD <- 3
        zero.COL <- "grey70"
        LWD <- 3
        COL <- rev(RColorBrewer::brewer.pal(11, type)[-6])
    } else if (tolower(type) == "user"){
        LWD <- 3
        zero.LWD <- 3
        if (length(col) >= 11){
            zero.COL <- col[1]
            COL <- col[2:11]
        } else {
            warning("Warning: the colour scheme 'type' was chosen as 'user', but less than 11 colours were passed in argument 'col'. Remember to also pass a colour for zero attachmen probability as the first colour in the vector. Using a custom palette instead.")
            zero.COL <- "grey70"
            COL <- terrain.colors(10)
        }
    }

    placements <- get.rogue.placement(tree, reftrees, rogues, outgroup)
    plcmnt.table <- NULL
    column.rogues <- NULL
    for (l in 1:length(placements)){
        plcmnt.table <- rbind(plcmnt.table, placements[[l]])
        column.rogues <- c(column.rogues, rep(names(placements)[l], dim(placements[[l]])[1]))
    }
    colnames(plcmnt.table) <- colnames(placements[[1]])
    plcmnt.table <- cbind(column.rogues, plcmnt.table)
    colnames(plcmnt.table)[1] <- "rogue"
    write.table(plcmnt.table, file = outfile.table, quote = F, row.names = F)

    edge.formats <- get.edge.formats(tree, placements, zero.col = zero.COL, col = COL, zero.lwd = zero.LWD, lwd = LWD, min.prob = min.prob)

    if (outfile.plot != "") {
        pdf(file=outfile.plot, width=15, height=21, onefile = T)
    }

    for (f in 1:length(placements)){
        rogue <- names(placements)[f]
        ff <- edge.formats[[f]]
        layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(5,1))
        par(plt = c(0.05, 0.9, 0, 0.9))
        phylogramm.plot.twocoloured (tree, edge.color = ff$cols.horizontal, vertical.edge.color = ff$cols.vertical, edge.lwd = ff$lwd.horizontal, vertical.edge.lwd = ff$lwd.vertical, tip.color = tip.color, edge.lty = 1, no.margin = FALSE, cex = cex.tips)
        title(main = rogue, cex = 1.5 * cex.tips)
        par(plt = c(0, 0.3, 0.4, 0.7))
        prob.scale(col = c(zero.COL, COL), breaks = c(0, min.prob, seq(0.1, 1, by = 0.1)), cex = cex.tips, title = "")
    }
    if (outfile.plot != "") {
        dev.off()
    }
    return(length(rogues))
}
#-------------------------------------------------------------------------------
#' Get Rogue Taxon Placements
#'
#' Function calculating the frequency of different placements of one or several
#' rogue taxa in a set of reference trees.
#'
#' @param tree an object of class \code{phylo}.
#' @param reftrees an object of class \code{multiphylo} containing a set of
#'   treesfrom a bootstrap or Bayesian analysis; they all have to contain the
#'   same tips, and the same tips as \code{tree} PLUS any rogue taxa to analyze.
#' @param rogues a vector of mode character containing the names of the tip or
#'   tips whose placement is to be calculated.
#' @param outgroup a vector of mode numeric or character specifying the new
#'   outgroup.
#'
#' @details
#' This function lists all groups appearing as sister groups to the rogue taxon
#'   (or taxa) in any of the reference trees, their respective frequencies,
#'   and whether or not they are present in the tree that is to be plotted.
#'
#' Make sure that the reference trees contain exactly the taxa in the tree,
#'   plus the rogue taxon or taxa.
#'
#' The reference trees are rooted before analysis, using either the specified
#'   outgroup or the first taxon in the tiplabels. A vector of outgroup taxa,
#'   either in numeric or character format, can be passed, in which case the
#'   function checks whether the specified outgroup is monophyletic. If not,
#'   a warning is issued and the first taxon among the outgroups is used to
#'   root the tree.
#'
#' @return
#' Returns a list whose elements correspond to a dataframe for each rogue,
#'   which contain comma-separated lists of sister taxa, their frequency among
#'   the reference trees, whether or not the sister group is monophyletic in
#'   the consensus tree, and finally its node index in the tree to plot. If a
#'   sister group does not appear in the tree, then the node corresponding to
#'   the most recent common ancestor is recorded.
#'
#' @references
#' Klopfstein, S., Spasojevic, T.(2018) ...
#'
#' @examples
#' ## example using random trees
#' reftrees <- ape::rmtree(5, 21, tip.label = c("rogue.taxon", paste("t",
#'   1:20, sep = "")))
#' reftrees <- c(rep(list(reftrees[[1]]), 60), rep(list(reftrees[[2]]), 25),
#'   rep(list(reftrees[[3]]), 10), reftrees)
#' tree <- ape::root(reftrees[[1]], outgroup = "t1")
#' placements <- get.rogue.placement(tree, reftrees, rogues = "rogue.taxon",
#'   outgroup = "t1")
#' str(placements)
#'
#' ## example from Klopfstein & Spasojevic 2018
#' data(consensus.tree)
#' data(fossil.trees)
#' data(taxon.lists)
#' placements <- get.rogue.placement(consensus.tree, fossil.trees,
#'   rogues = c(taxon.lists$fossils[1], taxon.lists$fossils[7]),
#'   outgroup = taxon.lists$outgroup[1])
#' str(placements)
#'
#-------------------------------------------------------------------------------
#' @export
get.rogue.placement <- function (tree, reftrees, rogues, outgroup = 1){

    if (is.numeric(outgroup)){
        if (any(outgroup) > length(tree$tip.label)){
            warning(paste("Warning: in 'get.rogue.placement', the outgroup
                          was specified using (a) numerical value(s) higher
                          than the number of taxa in the tree. Using the first
                          taxon as outgroup instead!"))
            outgroup <- tree$tip.label[1]
        } else {
            outgroup <- tree$tip.label[outgroup]
        }
    }
    rogue.placement <- NULL
    message("Calculations of rogue taxon placement in progress...")
    total <- length(rogues) * 2
    progress <- txtProgressBar(min = 0, max = total, style = 3)

    taxaToDrop <- setdiff(setdiff(reftrees[[1]]$tip.label, tree$tip.label), rogues)

    for (k in 1:length(rogues)){
        rogue <- rogues[k]
        rogue.sisters <- vector()

        for (i in 1:length(reftrees)){
            treeN <- reftrees[[i]]
            if (length(taxaToDrop) > 0){
                treeN <- ape::drop.tip(treeN, taxaToDrop)
            }
            if (is.monophyletic(treeN, outgroup)){
                treeN <- ape::root(treeN, outgroup = outgroup,
                                   resolve.root = TRUE)
            } else {
                warning("Warning: outgroup taxa not monophyletic in reference
                        tree number ", i, ". Using only first outgroup taxon
                        for rooting.")
                treeN <- ape::root(treeN, outgroup = outgroup[1],
                                   resolve.root = TRUE)
            }
            if (length(rogues) > 1){
                treeN <- ape::drop.tip(treeN, rogues[-k])
            }
            rogue.index <- match(rogue, treeN$tip.label)
            rogue.ancesNode <- treeN$edge[which(treeN$edge[,2] == rogue.index)]
            rogue.sis <- setdiff(list.descendents(treeN, rogue.ancesNode,
                                                  return.numeric = F), rogue)
            rogue.sis.contrary <- setdiff(setdiff(treeN$tip.label, rogue.sis), rogue)
            if(length(rogue.sis.contrary) == 1){
                rogue.sis <- rogue.sis.contrary
            }
            rogue.sisters[i] <- paste(sort(rogue.sis), collapse = ",")
        }
        Sys.sleep(0.1)
        setTxtProgressBar(progress, 2*k - 1)

        uniq.sisters <- unique(rogue.sisters)
        tempMat <- matrix(nrow = length(uniq.sisters), ncol = 4,
                          dimnames = list(NULL, c("sister.group", "frequency",
                          "present_in_tree", "containing.group.index")) )

        for (j in 1:length(uniq.sisters)){
            freq.sis.group <- sum(rogue.sisters == uniq.sisters[j]) / length(reftrees)
            groupN <- as.vector(unlist(strsplit(uniq.sisters[j], ",")))
            if (length(groupN) == 1){
                sis.index <- match(groupN, tree$tip.label)
                present.in.tree <- T
            } else {
                sis.index <- getMrca(tree, groupN)
                MRCA.desc <- list.descendents(tree, sis.index, return.numeric = F)
                if (setequal(groupN, MRCA.desc) ){
                    present.in.tree <- T
                } else {
                    present.in.tree <- F
                }
            }
            tempMat[j, ] <- c(uniq.sisters[j], freq.sis.group, present.in.tree, sis.index)
        }

        rogue.placement[[k]] <- as.data.frame(tempMat)
        names(rogue.placement)[k] <- rogue
        Sys.sleep(0.1)
        setTxtProgressBar(progress, 2*k)
    }
    close(progress)
    return(rogue.placement)
}
#------------------------------------------------------------------------------
#' Get Edge Formatting for Rogue Plots
#'
#' Given a tree and a list of taxon placements (as generated by the
#'   \code{\link{get.rogue.placement}}) function), this function generates a
#'   list with colours and line widths for the horizontal and vertical edges of
#'   a phylogram.
#'
#' @param tree an object of class \code{phylo}.
#' @param placements a list containing, for each rogue taxon, the sister
#'   groups of the rogue and their frequencies in the reference trees. For
#'   details about the format see the help of the
#'   \code{get.rogue.placement}) function.
#' @param zero.col colour used to reflect a placement probability below
#'   \code{min.prob}.
#' @param col a vector with (at least) 10 colours which reflect placement
#'   probabilities between \code{min.prob} and 1.
#' @param zero.lwd line width used to reflect a placement probability below
#'   \code{min.prob}. Defaults to 1.
#' @param lwd line width used to reflect a placement probability above
#'   \code{min.prob}. Defaults to 3.
#' @param min.prob minimum probability of attachment, below which an edge is
#'   considered as zero for plotting. Defaults to 1\% of reference trees.
#'
#' @details
#' \code{get.edge.formats} compares a list containing rogue taxon
#'   placement information (in the format generated by
#'   \code{\link{get.rogue.placement}}) with the tree to plot and generates a
#'   list of formating vectors for the horizontal and vertical branches of a
#'   phylogram. Both colour and line type of an edge can be set according to
#'   the rogue taxon placement probabilites.
#'
#' @return
#' Returns a list containing, for each rogue taxon, the horizontal and vertical
#'   colour vectors and horizontal and vertical line widths. Beware that
#'   the lengths of the horizontal vectors corresponds to the number of edges
#'   in the tree, while the vertical ones contain as many elements as the tree
#'   has internal nodes.
#'
#' @references
#' Klopfstein, S., Spasojevic, T.(2018) ...
#'
#' @examples
#' ## example using random trees
#' reftrees <- ape::rmtree(5, 21, tip.label = c("rogue.taxon", paste("t",
#'   1:20, sep = "")))
#' reftrees <- c(rep(list(reftrees[[1]]), 60), rep(list(reftrees[[2]]),
#'   25), rep(list(reftrees[[3]]), 10), reftrees)
#' tree <- reftrees[[1]]
#' placements <- get.rogue.placement(tree, reftrees, rogues =
#'   "rogue.taxon", outgroup = "t1")
#' edge.formats <- get.edge.formats(tree, placements, zero.col = "black",
#'   col = paste("grey", seq(from = 95, to = 20, by = -8), sep = ""),
#'   zero.lwd = 1, lwd = 3, min.prob = 0.01)
#' str(edge.formats)
#'
#' ## example from Klopfstein & Spasojevic 2018
#' data(consensus.tree)
#' data(fossil.trees)
#' data(taxon.lists)
#' placements <- get.rogue.placement(consensus.tree, fossil.trees,
#'   rogues = taxon.lists$fossils[7], outgroup = taxon.lists$outgroup[1])
#' edge.formats <- get.edge.formats(tree, placements, zero.col = "black",
#'   col = paste("grey", seq(from = 95, to = 20, by = -8), sep = ""),
#'   zero.lwd = 1, lwd = 3, min.prob = 0.01)
#' str(edge.formats)
#-------------------------------------------------------------------------------
#' @export
get.edge.formats <- function (tree, placements, zero.col = "black",
                              col = paste("grey", seq(from = 95, to = 15,
                              by = -8), sep = ""), zero.lwd = 1, lwd = 3,
                              min.prob = 0.01){

    edge.formats <- NULL

    for (i in 1:length(placements)){
        rogue <- names(placements)[i]
        test_table <- placements[[i]]
        test_table <- test_table[which(as.numeric(as.character(
                              test_table$frequency)) >= min.prob), ]

        cols.horizontal <- rep.int(zero.col, length(tree$edge[,1]))
        lwd.horizontal <- rep.int(zero.lwd, length(tree$edge[,1]))

        in.cons <- test_table[which(test_table$present_in_tree == T), ]
        if (length(in.cons) > 0){
            for (n in 1:dim(in.cons)[1]) {
                curr.nodeN <- in.cons$containing.group.index[n]
                pos <- match(curr.nodeN, tree$edge[,2])
                prob <- ceiling(as.numeric(as.character(in.cons$frequency[n])) * 10)
                cols.horizontal[pos] <- col[prob]
                lwd.horizontal[pos] <- lwd
            }
        }

        not.in.cons <- test_table[which(test_table$present_in_tree == F), ]
        cols.vertical <- rep.int(zero.col, tree$Nnode)
        lwd.vertical <- rep.int(zero.lwd, tree$Nnode)

        nodes.in.order <- sort(unique(tree$edge[,1]))

        if (length(not.in.cons) > 0) {
            tt <- not.in.cons[, c(2,4)]
            tt$frequency <- as.numeric(as.character(tt$frequency) )
            if (length(unique(tt$containing.group.index)) < dim(tt)[1]){
                tt <- aggregate(tt$frequency, by = list(corresp.node =
                                    tt$containing.group.index), FUN = sum)
            }
            for (cn in 1:dim(tt)[1]) {
                curr.nodeN <- tt$corresp.node[cn]
                pos <- match(curr.nodeN, nodes.in.order)
                prob <- ceiling(as.numeric(as.character(tt$x[cn])) * 10)
                cols.vertical[pos] <- col[prob]
                lwd.vertical[pos] <- lwd
            }
        }
        edge.formats[[i]] <- list(cols.horizontal, lwd.horizontal,
                                  cols.vertical, lwd.vertical)
        names(edge.formats[[i]] ) <- c("cols.horizontal", "lwd.horizontal",
                                       "cols.vertical", "lwd.vertical")
        names(edge.formats)[i] <- rogue

    }
    return(edge.formats)
}
#-------------------------------------------------------------------------------
#' Draw a Two-Coloured Phylogram
#'
#' Draw a phylogram with horizontal and vertical branches formated separately,
#'   for instance to create a rogue plot as introduced in Klopfstein &
#'   Spasojevic (2018).
#'
#' @param tree an object of class \code{phylo} to be plotted.
#' @param edge.color a vector of mode character giving the colours used to draw
#'   the HORIZONTAL edges of the phylogram. These are taken to be in the same
#'   order as the \code{edge} component of the \code{phylo} object. If fewer
#'   colours are given than the length of \code{edge}, then they are recycled.
#'   Defaults to 'black'.
#' @param vertical.edge.color same for the VERTICAL edges in the phylogram. If
#'   fewer colours are given than the number of internal nodes, then they are
#'   recycled. Defaults to 'grey50'.
#' @param edge.lwd a numeric vector giving the line widths of the edges of the
#'   plotted phylogeny. These are taken to be in the same order as the
#'   component \code{edge} of the \code{tree}. If fewer values are given than
#'   the length of \code{edge}, then these are recycled.
#' @param vertical.edge.lwd same as previous argument, but for the VERTICAL
#'   edges. If fewer values are given than the number of internal nodes, then
#'   these are recycled.
#' @param use.edge.length a logical indicating whether to use the edge lengths
#'   of the phylogeny to draw the branches (the default) or not (if FALSE).
#'   This option has no effect if the object of class 'phylo' has no
#'   \code{edge.length} element.
#' @param show.tip.label a logical indicating whether to show the tip labels
#'   on the phylogeny (defaults to TRUE, i.e. the labels are shown).
#' @param show.node.label a logical indicating whether to show the node labels
#'   on the phylogeny (defaults to FALSE, i.e. the labels are not shown).
#' @param tip.color the colours used for the tip labels, recycled if needed.
#' @param edge.lty a numeric vector giving the line type of the edges of the
#'   plotted phylogeny. Options are 1: plain, 2: dashed, 3: dotted, 4: dotdash,
#'   5: longdash, 6: twodash.
#' @param font an integer specifying the type of font for the labels: 1: plain
#'   text, 2: bold, 3: italic (the default), 4: bold italic.
#' @param cex a numeric value giving the scaling factor of the tip and node
#'   labels (Character EXpansion). The default is to take the current value
#'   from the graphical parameters.
#' @param adj a numeric specifying the justification of the text strings of the
#'   labels: 0: left-justification, 0.5: centering, or 1: right-justification.
#'   If NULL (the default), the value is set with respect to the orientation of
#'   the phylogram (see details).
#' @param srt a numeric giving how much the labels are rotated in degrees
#'   (negative values are allowed resulting in clock-like rotation); the value
#'   has an effect respectively to the value of direction.
#' @param no.margin a logical. If TRUE, the margins are set to zero and the
#'   plot uses all the space of the device.
#' @param root.edge a logical indicating whether to draw the root edge
#'   (defaults to FALSE); this has no effect if \code{use.edge.length = FALSE}.
#' @param label.offset a numeric giving the space between the nodes and the
#'   tips of the phylogeny and their corresponding labels.
#' @param underscore a logical specifying whether the underscores in tip labels
#'   should be written as spaces (the default) or left as they are (if TRUE).
#' @param direction a character string specifying the direction of the tree.
#'   Four values are possible: "rightwards" (the default),
#'   "leftwards", "upwards", and "downwards".
#' @param plot a logical controlling whether to draw the tree. If FALSE, the
#'   graphical device is set as if the tree was plotted, and the coordinates
#'   are saved as well.
#' @param node.depth an integer value (1 or 2) used if edge lengths are not
#'   used to plot the tree; 1: the node depths are proportional to the number
#'   of tips descending from each node (the default), 2: they are evenly
#'   spaced.
#' @param align.tip.label a logical value or an integer. If TRUE, the tips are
#'   aligned and dotted lines are drawn between the tips of the tree and the
#'   labels. If an integer, the tips are aligned and this gives the type of the
#'   lines (lty).
#' @param x.lim a numeric vector of length one or two giving the limit(s) of
#'   the x-axis. If NULL, this is computed with respect to various parameters
#'   such as the string lengths of the labels and the branch lengths. If a
#'   single value is given, this is taken as the upper limit. Defaults to
#'   \code{NULL}.
#' @param y.lim same as above, but for the y-axis.
#' @param ... further arguments passed to \code{\link[ape]{plot.phylo}}.
#'
#' @details
#' This function plots a tree in 'rectangular phylogram' format, allowing to
#'   pass colour and line typevectors for the horizontal and vertical edges
#'   separately. The colour vector for horizontal edges is passed in the same
#'   order as the edges ('edge.color'), and the vector for the vertical edges
#'   in the same order as the internal node indices.
#'
#' Typically, one would first run \code{\link{get.rogue.placement}} on a tree
#'   and corresponding reference trees, then obtain edge formatting vectors
#'   from \code{\link{get.edge.formats}}, and finally pass those to the current
#'   function to generate the plot. Instead the wrapper function
#'   \code{\link{create.rogue.plot}} can be used.
#'
#' This function is based on the \code{\link[ape]{plot.phylo}} function in
#'   the \code{ape} package version 5.1 by E. Paradis.
#'
#' @return
#' \code{phylogramm.plot.twocoloured} returns invisibly a list with
#'    components whose values are those used for the current plot (see the
#'    help of \code{plot.phylo} for details).
#'
#' @references
#' Klopfstein, S., Spasojevic, T. (2018).
#'
#' @examples
#' ## general examples
#' tree <- ape::rtree(12, rooted = TRUE)
#' cols.horiz <- rep("black", length(tree$edge.length))
#' cols.verti <- rainbow(tree$Nnode)
#' phylogramm.plot.twocoloured(tree, edge.color = cols.horiz,
#'   vertical.edge.color = cols.verti)
#' cols.horiz <- rainbow(length(tree$edge.length))
#' cols.verti <- paste("grey", seq(from = 95, to = 20,
#'   by = -floor(75 / tree$Nnode)), sep = "")
#' phylogramm.plot.twocoloured(tree, edge.color = cols.horiz,
#'   vertical.edge.color = cols.verti)
#'
#' ## example of a Rogue Plot
#' reftrees <- ape::rmtree(5, 21, tip.label = c("rogue.taxon", paste("t",
#'   1:20, sep = "")))
#' reftrees <- c(rep(list(reftrees[[1]]), 60), rep(list(reftrees[[2]]),
#'   25), rep(list(reftrees[[3]]), 10), reftrees)
#' tree <- reftrees[[1]]
#' placements <- get.rogue.placement(tree, reftrees, rogues =
#'   "rogue.taxon", outgroup = "t1")
#' tree <- ape::drop.tip(tree, tip = "rogue.taxon")
#' edge.formats <- get.edge.formats(tree, placements, zero.col = "black",
#'   col = paste("grey", seq(from = 95, to = 20, by = -8), sep = ""),
#'   zero.lwd = 1, lwd = 3, min.prob = 0.01)
#' ff <- edge.formats[[1]]
#' phylogramm.plot.twocoloured(tree, edge.color = ff$cols.horizontal,
#'   vertical.edge.color = ff$cols.vertical, edge.lwd = ff$lwd.horizontal,
#'   vertical.edge.lwd = ff$lwd.vertical, tip.color = "black", edge.width = 2,
#'   no.margin = TRUE, cex = 0.8)
#'
#-------------------------------------------------------------------------------
#' @export
phylogramm.plot.twocoloured <- function (tree, edge.color = "black",
                                         vertical.edge.color = "grey50", edge.lwd = 2, vertical.edge.lwd = 2,
                                         use.edge.length = TRUE, show.tip.label = TRUE,
                                         show.node.label = FALSE, tip.color = "black", edge.lty = 1,
                                         font = 3, cex = par("cex"), adj = NULL, srt = 0, no.margin = FALSE,
                                         root.edge = FALSE, label.offset = 0, underscore = FALSE,
                                         direction = "rightwards",  plot = TRUE, node.depth = 1,
                                         align.tip.label = FALSE, x.lim = NULL, y.lim = NULL, ...)
{
    Ntip <- length(tree$tip.label)
    if (Ntip < 2) {
        stop("found less than 2 tips in the tree")
    }
    if (any(tabulate(tree$edge[, 1]) == 1))
        stop("there are single (non-splitting) nodes in your tree; you may need
             to use collapse.singles()")

    .nodeHeight <- function (edge, Nedge, yy) .C(ape::node_height,
                                                 as.integer(edge[,1]), as.integer(edge[, 2]),
                                                 as.integer(Nedge), as.double(yy))[[4]]
    .nodeDepth <- function (Ntip, Nnode, edge, Nedge, node.depth)
                                            .C(ape::node_depth, as.integer(Ntip),
                                            as.integer(edge[,1]), as.integer(edge[, 2]),
                                            as.integer(Nedge), double(Ntip + Nnode),
                                            as.integer(node.depth))[[5]]
    .nodeDepthEdgelength <- function (Ntip, Nnode, edge, Nedge, edge.length)
        .C(ape::node_depth_edgelength, as.integer(edge[, 1]),
           as.integer(edge[, 2]), as.integer(Nedge),
           as.double(edge.length), double(Ntip + Nnode))[[5]]
    Nedge <- dim(tree$edge)[1]
    Nnode <- tree$Nnode
    if (any(tree$edge < 1) || any(tree$edge > Ntip + Nnode))
        stop("tree badly conformed; cannot plot. Check the edge matrix.")

    ROOT <- Ntip + 1
    direction <- match.arg(direction, c("rightwards", "leftwards", "upwards",
                                        "downwards"))
    if (is.null(tree$edge.length))
        use.edge.length <- FALSE
    if (is.numeric(align.tip.label)) {
        align.tip.label.lty <- align.tip.label
        align.tip.label <- TRUE
    } else {
        if (align.tip.label)
            align.tip.label.lty <- 3
    }
    if (align.tip.label) {
        if (!use.edge.length || ape::is.ultrametric(tree))
            align.tip.label <- FALSE
    }
    if (!use.edge.length || is.null(tree$root.edge) || !tree$root.edge)
        root.edge <- FALSE
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- tree$edge
    phyOrder <- attr(tree, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
        tree <- ape::reorder.phylo(tree)
        if (!identical(tree$edge, xe)) {
            ereorder <- match(tree$edge[, 2], xe[, 2])
            if (length(edge.color) > 1) {
                edge.color <- rep(edge.color, length.out = Nedge)
                edge.color <- edge.color[ereorder]
            }
            if (length(edge.lty) > 1) {
                edge.lty <- rep(edge.lty, length.out = Nedge)
                edge.lty <- edge.lty[ereorder]
            }
            if (length(edge.lwd) > 1) {
                edge.lwd <- rep(edge.lwd, length.out = Nedge)
                edge.lwd <- edge.lwd[ereorder]
            }
        }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- tree$edge[tree$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip

    z <- ape::reorder.phylo(tree, order = "postorder")
    yy <- .nodeHeight(z$edge, Nedge, yy)

    if (!use.edge.length) {
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth) - 1
        xx <- max(xx) - xx
    } else {
        xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, z$edge.length)
    }

    if (!horizontal) {
        tmp <- yy
        yy <- xx
        xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
        if (direction == "rightwards")
            xx <- xx + tree$root.edge
        if (direction == "upwards")
            yy <- yy + tree$root.edge
    }

    if (no.margin)
        par(mai = rep(0, 4))
    if (show.tip.label)
        nchar.tip.label <- nchar(tree$tip.label)
    max.yy <- max(yy)
    if (is.null(x.lim)) {
        if (horizontal) {
            x.lim <- c(0, NA)
            pin1 <- par("pin")[1]
            strWi <- strwidth(tree$tip.label, "inches", cex = cex)
            xx.tips <- xx[1:Ntip] * 1.04
            alp <- try(uniroot(function (a) max(a * xx.tips + strWi) - pin1,
                               c(0, 1e+06))$root, silent = TRUE)
            if (is.character(alp)) {
                tmp <- max(xx.tips)
                if (show.tip.label)
                    tmp <- tmp * 1.5
            } else {
                tmp <- if (show.tip.label)
                    max(xx.tips + strWi/alp)
                else max(xx.tips)
            }
            if (show.tip.label)
                tmp <- tmp + label.offset
            x.lim[2] <- tmp
        } else x.lim <- c(1, Ntip)

    } else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (!horizontal)
            x.lim[1] <- 1
        else -1
    }
    if (direction == "leftwards")
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (horizontal) {
            y.lim <- c(1, Ntip)
        } else {
            y.lim <- c(0, NA)
            pin2 <- par("pin")[2]
            strWi <- strwidth(tree$tip.label, "inches", cex = cex)
            yy.tips <- yy[1:Ntip] * 1.04
            alp <- try(uniroot(function (a) max(a * yy.tips + strWi) - pin2,
                               c(0, 1e+06))$root, silent = TRUE)
            if (is.character(alp)) {
                tmp <- max(yy.tips)
                if (show.tip.label)
                    tmp <- tmp * 1.5
            }
            else {
                tmp <- if (show.tip.label)
                    max(yy.tips + strWi/alp)
                else max(yy.tips)
            }
            if (show.tip.label)
                tmp <- tmp + label.offset
            y.lim[2] <- tmp
        }

    } else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (horizontal)
            y.lim[1] <- 1
        else -1
    }
    if (direction == "downwards")
        yy <- y.lim[2] - yy
    if (root.edge) {
        if (direction == "leftwards")
            x.lim[2] <- x.lim[2] + tree$root.edge
        if (direction == "downwards")
            y.lim[2] <- y.lim[2] + tree$root.edge
    }
    asp <- NA
    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "",
                 ylab = "", axes = FALSE, asp = asp, ...)
    if (plot) {
        if (is.null(adj)) {
            adj <- if (direction == "leftwards") {
                1
            } else 0
        }
        if (show.tip.label) {
            MAXSTRING <- max(strwidth(tree$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
            }
            if (!horizontal) {
                psr <- par("usr")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                    loy <- -loy
                    srt <- 180 + srt
                }
            }
        }

        nodes <- (Ntip + 1):(Ntip + Nnode)
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp
        }
        x0v <- xx[nodes]
        y0v <- y1v <- numeric(Nnode)
        for (i in 1:Nnode) {
            descs <- tree$edge[,2][which(tree$edge[, 1] == i + Ntip)]
            tmp <- range(yy[descs])
            y0v[i] <- tmp[1]
            y1v[i] <- tmp[2]
        }
        x0h <- xx[tree$edge[, 1]]
        x1h <- xx[tree$edge[, 2]]
        y0h <- yy[tree$edge[, 2]]

        #format of horizontal branches
        Nedge <- dim(tree$edge)[1]
        edge.color <- rep(edge.color, length.out = Nedge)
        edge.lty <- rep(edge.lty, length.out = Nedge)
        edge.lwd <- rep(edge.lwd, length.out = Nedge)

        #format of vertical branches
        if (length(vertical.edge.lwd) == Nnode ){
            lwd.v <- vertical.edge.lwd
        } else {
            lwd.v <- rep(vertical.edge.lwd, length.out = Nnode)
        }
        if (length(vertical.edge.color) == Nnode ){
            color.v <- vertical.edge.color
        } else {
            color.v <- rep(vertical.edge.color, length.out = Nnode)
        }

        if (horizontal) {
            segments(x0h, y0h, x1h, y0h, col = edge.color, lwd = edge.lwd,
                     lty = edge.lty)
            segments(x0v, y0v, x0v, y1v, col = color.v, lwd = lwd.v, lty = edge.lty)
        } else {
            segments(y0h, x0h, y0h, x1h, col = edge.color, lwd = edge.lwd,
                     lty = edge.lty)
            segments(y0v, x0v, y1v, x0v, col = color.v, lwd = lwd.v, lty = edge.lty)
        }


        if (root.edge) {
            rootcol <- if (length(edge.color) == 1) {
                edge.color
            } else "black"

            switch(direction, rightwards = segments(0, yy[ROOT], tree$root.edge,
                                                    yy[ROOT], col = rootcol),
                   leftwards = segments(xx[ROOT], yy[ROOT],
                                        xx[ROOT] + tree$root.edge, yy[ROOT],
                                        col = rootcol),
                   upwards = segments(xx[ROOT], 0, xx[ROOT], tree$root.edge,
                                      col = rootcol),
                   downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT],
                                        yy[ROOT] + tree$root.edge, col = rootcol))
        }
        if (show.tip.label) {
            if (is.expression(tree$tip.label))
                underscore <- TRUE
            if (!underscore)
                tree$tip.label <- gsub("_", " ", tree$tip.label)

            if (align.tip.label) {
                xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]),
                                 leftwards = min(xx[1:Ntip]),
                                 upwards = xx[1:Ntip], downwards = xx[1:Ntip])
                yy.tmp <- switch(direction, rightwards = yy[1:Ntip],
                                 leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]),
                                 downwards = min(yy[1:Ntip]))
                segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp,
                         lty = align.tip.label.lty)
            }
            else {
                xx.tmp <- xx[1:Ntip]
                yy.tmp <- yy[1:Ntip]
            }
            text(xx.tmp + lox, yy.tmp + loy, tree$tip.label, adj = adj, font = font, srt = srt, cex = cex, col = tip.color)

        }
        if (show.node.label)
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)],
                 tree$node.label, adj = adj, font = font, srt = srt,
                 cex = cex)
    }
    L <- list(type = "phylogram", use.edge.length = use.edge.length,
              node.pos = 1, node.depth = node.depth,
              show.tip.label = show.tip.label,
              show.node.label = show.node.label, font = font, cex = cex,
              adj = adj, srt = srt, no.margin = no.margin,
              label.offset = label.offset, x.lim = x.lim, y.lim = y.lim,
              direction = direction, tip.color = tip.color, Ntip = Ntip,
              Nnode = Nnode, root.time = tree$root.time,
              align.tip.label = align.tip.label)
    #assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)),
     #      envir = .PlotPhyloEnv)
    invisible(L)
}
#-------------------------------------------------------------------------------
#' Create Scale for Probabilities
#'
#' This function plots a scale from zero to one next to a vertical colour code.
#'
#' @param col a vector of mode character giving the colours used to reflect
#'   probabilities. At least two colours have to be passed.
#' @param breaks (optional): the limits of the categories of the scale. If
#'   \code{breaks} are not provided, they are spaced equally from zero to one,
#'   in order to match the number of colours passed in \code{col}.
#' @param title a title for the scale. Pass \code{'""'} to suppress the title.
#' @param cex a numeric value giving the scaling factor of the axis labels
#'   (Character EXpansion). The default is to take the current value from the
#'   graphical parameters. If printed, the title of the scale is chosen to be
#'   1.2x as large.
#'
#' @details
#' \code{prob.scale} prints a probability scale using the colour vector
#'   passed to it. The scale is not added to the current plot but printed
#'   separately. Use \code{\link{layout}} or similar to create subplot areas if
#'   intending to plot the scale next to a graph.
#'
#' @examples
#' ## example of an evenly spaced grey scale:
#' prob.scale(col = paste("grey", 10*(9:2), sep = ""), cex = 0.9)
#'
#' ## example of an evenly spaced rainbow scale:
#' prob.scale(col = rainbow(n = 8), cex = 0.9)
#'
#' ## example of an unevenly scale with user-defined breaks:
#' prob.scale(col = c("grey80", rainbow(n = 5, start = 0.35, end = 0.7)),
#'   breaks = c(0, 0.05, seq(from = 0.2, to = 1.0, by = 0.2)), cex = 0.9)
#-------------------------------------------------------------------------------
#' @export
prob.scale <- function (col, breaks, title = "probability", cex = par("cex")) {

    if (missing(col)){
        stop("Scale cannot be drawn unless a vector of colours of at least
             length 2 is passed to the function. Exiting...")
    }
    if (length(col) < 2){
        stop("Scale cannot be drawn unless a vector of colours of at least
             length 2 is passed to the function. Exiting...")
    }
    if (!missing(breaks)){
        if (length(breaks) != (length(col) +1)){
            warning("There have to be exactly one more breaks than colours for
                    the scale. Using equally spaced breaks instead...")
            breaks <- seq(0, 1, by = 1 / length(col))
        }
    } else {
        breaks <- seq(0, 1, by = 1 / length(col))
    }

    col.rect <- vector(mode="list", length(col))
    for (i in 1:length(col.rect)){
        col.rect[[i]] <- c(i -1, i, i, i -1) / length(col)
    }
    ylim <- c(0, 1.1)
    xlim <- c(0,1)

    plot(1.1, 1.1, type = "n", ylim = ylim, xlim = xlim, xaxt = "n", yaxt = "n",
         xlab = "", ylab = "", main = title, cex.main = cex * 1.2, bty = "n")
    axis(4, at = seq(from = 0, to = 1, by = 1 / length(col)),
         labels = sprintf("%.2f", breaks), tick = F, las = 1, cex.axis = cex)

    for (i in 1:length(col.rect)){
        polygon(c(0,0,1,1), col.rect[[i]], col=col[i], border=NA)
    }
}
#-------------------------------------------------------------------------------
# internal
is.monophyletic <- function (tree, tips) {
    if (length(tips) == 1){
        return(T)
    }
    MRCA <- getMrca(tree, tips)
    descs <- list.descendents(tree, MRCA)
    if (length(tips) == length(descs)) {
        return(T)
    } else {
        return(F)
    }
}
#-------------------------------------------------------------------------------
# internal
list.descendents <- function (tree, nodeIndex, return.numeric = T){
    descs = NULL
    if (nodeIndex <= length(tree$tip.label) ){
        descs <- c(descs, nodeIndex)
    } else {
        children <- tree$edge[which(tree$edge[, 1] == nodeIndex), 2]
        for (i in 1:length(children)){
            descs <- c(descs, list.descendents(tree, children[i]))
        }
    }
    if (return.numeric) {
        return(descs)
    } else {
        return(tree$tip.label[descs])
    }
}
#-------------------------------------------------------------------------------
# internal
getMrca <- function (tree, tips){
    if (length(tips) < 2){
        stop("Error: at least two tips needed to calculate MRCA for. Exiting...")
    }
    ancestors <- rep(list(NULL), length(tips))
    root.ind <- as.integer(tree$edge[, 1][!match(tree$edge[, 1], tree$edge[, 2], 0)][1])
    for(i in 1:length(tips)){
        index <- which(tree$tip.label == tips[i])
        if (length(index) != 1){
            stop(paste("Error: no tip corresponding to taxon '", tips[i], " found in tree. Exiting...", sep = ""))
        }
        while (index != root.ind){
            index <- tree$edge[which(tree$edge[, 2] == index), 1]
            ancestors[[i]] <- c(ancestors[[i]], index)
        }
    }
    inters <- ancestors[[1]]
    for( j in 2:length(ancestors)){
        inters <- intersect(inters, ancestors[[j]])
    }
    return(inters[1])
}
#-------------------------------------------------------------------------------
#internal
check.alignment.format <- function(alignment, function.name = "function"){
    if (is.list(alignment)){
        align.mat <- matrix(nrow = length(alignment), ncol = length(alignment[[1]]), dimnames = list(names(alignment), NULL))
        for (a in 1:length(alignment)){
            align.mat[a, ] <- alignment[[a]]
        }
    } else if (is.matrix(alignment)){
        align.mat <- alignment
    } else {
        stop(paste("Error in '", function.name, "': alignment not passed in right format - has to be a list or matrix. Exiting...", sep = ""))
    }
    return(align.mat)
}
#-------------------------------------------------------------------------------
