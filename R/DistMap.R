DistMap <- setClass(Class = "DistMap",
                    slots = c(data = "matrix",
                              insitu.matrix = "matrix",
                              mcc.scores = "matrix")
                    )

setMethod("initialize",
          "DistMap",
          function (.Object, data=matrix(), insitu.matrix=matrix(), mcc.scores=matrix()) {
            .Object@data = data
            .Object@insitu.matrix = insitu.matrix
            .Object@mcc.scores = mcc.scores
            .Object
          })

