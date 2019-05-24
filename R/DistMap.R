DistMap <- setClass(Class = "DistMap",
                    slots = c(raw.data = "matrix",
                              data = "matrix",
                              binarized.data = "matrix",
                              insitu.matrix = "matrix",
                              geometry = "matrix",
                              mcc.scores = "matrix")
                    )

setMethod("initialize",
          "DistMap",
          function(.Object, raw.data=matrix(), data=matrix(), binarized.data=matrix(),
                   insitu.matrix=matrix(), geometry=matrix(), mcc.scores=matrix()) {
            .Object@raw.data = raw.data
            .Object@data = data
            .Object@binarized.data = binarized.data
            .Object@insitu.matrix = insitu.matrix
            .Object@geometry = geometry
            .Object@mcc.scores = mcc.scores
            .Object
          })

#' Binarize single cell data
#'
#' Binarizes the single cell data contained in the \code{DistMap} object. This is done in three steps:
#' (i) gene thresholds are chosen for a given quantile, (ii) the gene correlations of the binarized
#' data are compared against those of the \code{insitu.matrix} and a score is computed, (iii) the
#' binarized data based on the minimal score is computed and added to the \code{DistMap} object.
#'
#' @param object A \code{DistMap} object.
#' @param quantiles A range of quantiles.
#' @return A \code{DistMap} object with binarized data.
setGeneric(name = "binarizeSingleCellData",
           def = function(object, quantiles=seq(0.15, 0.16, 0.01)) {
             standardGeneric("binarizeSingleCellData")
           })
setMethod(f = "binarizeSingleCellData",
          signature = "DistMap",
          function(object, quantiles) {
            rsme.cor <- function(x, y) {
              return (sqrt(mean((x[lower.tri(x)] - y[lower.tri(y)])^2)))
            }

            insitu.genes <- colnames(object@insitu.matrix)
            insitu.cor <- cor(object@insitu.matrix)
            rsme.scores = data.frame("quantile" = quantiles,
                                     "score" = 0)

            for (quantile in quantiles) {
              gene.thresholds <- sapply(insitu.genes, function(gene) quantile(as.numeric(object@data[gene, object@data[gene, ] > 0]), quantile))
              names(gene.thresholds) <- gsub("\\..*", "", names(gene.thresholds))
              binarized.data <- object@data[insitu.genes, ]
              binarized.data[binarized.data >= 0] <- 0
              binarized.data[apply(object@data[insitu.genes, ], 2, function(cell) cell > gene.thresholds)] <- 1
              rsme.scores$score[which(rsme.scores$quantile == quantile)] <- rsme.cor(insitu.cor, cor(t(binarized.data)))
            }
            best.quantile <- rsme.scores$quantile[which.min(rsme.scores$score)]
            best.gene.thresholds <- sapply(insitu.genes, function(gene) quantile(as.numeric(object@data[gene, object@data[gene, ] > 0]),
                                                                                 best.quantile))
            names(best.gene.thresholds) <- gsub("\\..*", "", names(best.gene.thresholds))

            binarized.data <- object@data[insitu.genes, ]
            binarized.data[binarized.data >= 0] <- 0
            binarized.data[apply(object@data[insitu.genes, ], 2, function(cell) cell > best.gene.thresholds)] <- 1
            object@binarized.data <- binarized.data

            return (object)
          })


#' Compute Matthews Correlation Coefficients
#'
#' Computes MCC scores and adds them to the \code{DistMap} object. Requires the \code{binarized.data}
#' to be computed.
#'
#' @param object A \code{DistMap} object.
#' @param cell A number denoting the cell.
#' @return The MCC scores between the given cell and all bins.
setGeneric(name = "computeMCC",
           def = function(object, cell) {
             standardGeneric("computeMCC")})
setMethod(f = "computeMCC",
          signature = "DistMap",
          function(object, cell) {
            mtrx <- -sweep(object@insitu.matrix, 2, 2*object@binarized.data[, cell], '-')
            FP = rowSums(mtrx == 2)
            TP = rowSums(mtrx == 1)
            TN = rowSums(mtrx == 0)
            FN = rowSums(mtrx == -1)
            MCC <- (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            MCC[is.na(MCC)] <- 0

            return (MCC)
          })


#' Map cells
#'
#' Maps single cells onto the reference database.
#'
#' @param object A \code{DistMap} object.
#' @return A \code{DistMap} object with computed MCC scores.
setGeneric(name = "mapCells",
           def = function(object) {
             standardGeneric("mapCells")})
setMethod(f = "mapCells",
          signature = "DistMap",
          function(object) {
            mcc.scores <- sapply(1:(dim(object@data)[2]), function(cell) computeMCC(object, cell))
            mcc.scores <- exp(mcc.scores)
            object@mcc.scores <- mcc.scores
            return (object)
          })


#' Compute vISH
#'
#' Computes a virtual in situ hybridization for a given gene and a given threshold.
#'
#' @param object A \code{DistMap} object.
#' @param gene A gene symbol.
#' @param threshold A threshold between 0 and 1.
#' @return A vector of gene expression in every bin.
setGeneric(name = "computeVISH",
           def = function(object, gene, threshold=0.75) {
             standardGeneric("computeVISH")})
setMethod(f = "computeVISH",
          signature = "DistMap",
          function(object, gene, threshold) {
            gene.expr <- as.numeric(object@data[gene, ])
            b1 <- sweep(object@mcc.scores, 2, gene.expr, '*')
            gene.expr[gene.expr > 0] <- 1
            b2 <- sweep(object@mcc.scores, 2, gene.expr, '*')
            q <- rowSums(b1)/rowSums(b2)
            q[is.na(q)] <- 0
            q <- q/(1+q)
            q[q < quantile(q, threshold)] <- 0
            return (as.numeric(q))
          })

#' Compute gene gradient
#'
#' Computes a gradient of a gene.
#'
#' @param object A \code{DistMap} object.
#' @param gene The gene symbol.
#' @param threshold A threshold between 0 and 1.
setGeneric(name = "computeGeneGradient",
           def = function(object, gene, threshold=0.75, type='DV',
                          slide.x=0, slide.y=0) {
             standardGeneric("computeGeneGradient")})
setMethod(f = "computeGeneGradient",
          signature = "DistMap",
          function(object, gene, threshold, type, slide.x, slide.y) {
            vISH = computeVISH(object, gene, threshold)
            mtrx <- sweep(object@mcc.scores - min(object@mcc.scores),
                          MARGIN = 2, as.numeric(object@raw.data[gene, ]), '*')
            umis = rowMeans(mtrx) * vISH
            df <- data.frame("x"=100*(object@geometry[, 'x'] - min(object@geometry[, 'x']))/max(object@geometry[, 'x'] - min(object@geometry[, 'x'])),
                             "y"=100*(object@geometry[, 'z'] - min(object@geometry[, 'z']))/max(object@geometry[, 'z'] - min(object@geometry[, 'z'])),
                             "u"=exp(umis))

            df <- df[df$x >= slide.x & df$x <= 100 - slide.x, ]
            df <- df[df$y >= slide.y & df$y <= 100 - slide.y, ]

            if (type == 'DV') {
              ggplot(df, aes(y=u, x=y)) + geom_smooth() + coord_flip() + xlab('DV') + ylab('UMIs') + ggtitle(gene) + scale_x_continuous(limits = c(0, 100))
            }
            else {
              ggplot(df, aes(y=u, x=x)) + geom_smooth() + ylab('UMIs') + xlab('AP') + ggtitle(gene) + scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) + scale_y_continuous(expand=c(0, 0))
            }
          })

#' In Silico Dissect
#'
#' Performs in silico dissection of the embryo according to lisst of genes.
#'
#' @param object A \code{DistMap} object.
#' @param gene.sets A list of gene sets according to which the in silico dissection
#' will be performed. Each of the gene sets is given as a vector of gene symbols of
#' minimal length 2.
#' @return A \code{matrix} with scores for every cell against each gene set.
setGeneric(name = "inSilicoDissect",
           def = function(object, gene.sets) {
             standardGeneric("inSilicoDissect")
           })
setMethod(f = "inSilicoDissect",
          signature = "DistMap",
          function(object, gene.sets) {
            scores <- lapply(gene.sets, function(genes) apply(object@data[genes, ], 2, function(x) mean(x, na.rm = T)))
            scores <- matrix(unlist(scores), nrow=1297)
            scores <- t(apply(apply(scores, 2, scale), 1, scale))
            row.names(scores) <- names(object@raw.data)
            return(scores)
          })



