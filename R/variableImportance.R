#' Coefficent of correlation between observed
#' and predicted values
#'
#' The correlation coefficent is computed
#' between observed values (if present) or
#' against prediction obtained from fitted
#' values.
#'
#' @param object An object of \code{\link{Modelling-class}}.
#' @param obs an optional vector of observed values
#' @param newdata an optional data.frame containing all variables
#'  used in the model. If omitted, the fitted values are used
#' @param method the method used for correlation. By default, is used "spearman" method.
#' For more detail see {\code{\link[stats]{cor}}} function.
#'
#' @return the R value
#' @export
#'
#' @examples

setGeneric("getR", function(object, obs = NULL, newdata = NULL,
                            method = "spearman")
  standardGeneric("getR"))

setMethod("getR", "Modelling", function(object, obs = NULL, newdata = NULL,
                                        method = "spearman"){
  if(is.null(obs))
    obs <- object@yobs

  pred <- predict(object, newdata)

  cor(obs, pred, method = method)
})


#' Variable importance through shuffle method
#'
#' The importance of each variable is computed
#' by randomizing its values and predicting.
#' A reference RMSE is computed using the
#' full-model whereas the RMSE of the shuffle
#' model is used for estimating variable importance.
#'
#' @param object An object of \code{\link{Modelling-class}}.
#' @param obs optional vector of real values
#' @param newdata an optional data.frame containing all variables
#'  used in the object. If omitted, the fitted values are used
#' @param shuffles number of shuffles per variable
#' @param res.fun function for returning a single value per
#'  variable after shuffling
#'
#' @return a vector containing the importance of each variable
#'  normalized by the RMSE of the full-model
#'
#' @export
#'
#' @examples

setGeneric("variableImp", function(object, obs = NULL, newdata = NULL,
                                        shuffles = 500,
                                        res.fun = mean)
  standardGeneric("variableImp"))

setMethod("variableImp", "Modelling", function(object, obs = NULL, newdata = NULL,
                        shuffles = 500,
                        res.fun = mean){
  if(is.null(newdata))
    newdata <- object@varobs
  # Get reference error
  refRMSE <- error(object, obs = obs, newdata = newdata, type = "rmse")
  vars <- colnames(object@varobs)

  # Looping through all variables
  imp <- sapply(vars, function(v){

    # Repeted shuffling
    res <- sapply(1:shuffles, function(s){
      val.mod <- newdata
      val.mod[,v] <- sample(val.mod[,v], nrow(val.mod))

      # Calculating RMSE with shuffled variable
      shuffledRMSE <- error(object, obs = obs, newdata = val.mod, type = "rmse")

      # Returning variable importance
      (shuffledRMSE - refRMSE)/refRMSE
    })
    res.fun(res)
  })
  imp
}
)

#' Fitness dual model function
#'
#' @description The fitness dual model function has as input moldels variable
#' representing a potential solution and returns a numerical value
#' describing its ''fitness''. This function could be used in a genetic algorithm (
#' for more details see \code{\link[GA]{ga}})
#'
#' @param vars models variables
#'
#' @param x A data.frame where samples are in rows
#'  and features are in columns.
#'
#' @param y A data.frame containing validation value
#'
#' @param p A numeric value corresponding to randomly sample
#' a percentage of rows within a data frame
#'
#' @param cntr Control the computational nuances of the train function
#'  (for more detail see \code{\link{trainmod}}).
#' By default, is used the cross-validation method \code{"LOOCV"}.
#'
#' @param method A string specifying which classification and/or
#'  regression model to use (for details see \code{\link{trainmod}}.
#'  By default, is used the method \code{"rf"}.
#'  Possible values are found using names (\code{\link[caret]{getModelInfo}}).
#'
#' @return numeric string describing putative variables ''fitness''.
#' @export
#'
#' @examples

dualModelFitness <- function(vars, x, y, p = .8,
                             cntr = trainControl(method = "LOOCV"),
                             method = "rf"){
   if(all(vars == 0))
       return(0)

     t <- sample(1:nrow(x), size = floor(p*nrow(x)))

     vars <- which(vars == 1)
     n <- colnames(x)[vars]

     t.x <- setNames(data.frame(x[t, vars]), n)
     t.y <- as.vector(y[t])

     v.x <- setNames(data.frame(x[-t, vars]), n)
     v.y <- as.vector(y[-t])

     model <- trainmod(x = t.x, y = t.y, cntr = cntr,
                       method = method)

     1/(error(model, obs = v.y, newdata = v.x, type = "rmse"))
}

#' Genetic Algorithm
#'
#' @description for more details see \code{\link[GA]{ga}}
#'
#' @param data_x a data.frame
#' @param data_y numeric string
#' @param train a df.
#' @param type The type of genetic algorithm depending
#' on the nature of decision variables. They cpul be "binary" (default), "real-valued"
#' and "permutation". For more details see \code{\link[GA]{ga}}
#' @param fitness
#' @param nBits
#' @param popSize The population size.
#' @param elitism The number of best fitness individuals to survive at each generation
#' @param run
#' @param maxiter
#' @param monitor
#' @param parallel
#' @param names
#' @param seed
#' @param ... Additional arguments to be passed to the fitness function
#'
#' @return
#'
#' @export
#'
#' @import GA
#'
#' @examples
gadualModel <- function(data_x, data_y, train, type = "binary",
                        fitness = dualModelFitness,
                        nBits = ncol(data_x), popSize = 20,
                        elitism = 2, run = 10, maxiter = 500,
                        monitor = plot, parallel = T,
                        names = colnames(data_x),
                        seed = 11, ...) {

  set.seed(seed)
  ga1 <- ga(type = type, fitness = fitness, x = data_x, y = data_y,
            nBits = nBits, popSize = popSize, elitism = elitism,
            run = run, maxiter = maxiter, monitor = monitor,
            parallel = parallel, names = names, seed = seed, ...)

  sol <- which(ga1@solution[1,] == 1)
  n <- colnames(train)[sol]
  data_x <- setNames(data.frame(train[,sol]), n)

  model <- trainmod(x = data_x, y = data_y)

  v <- setNames(data.frame(validation[,sol]), n)

  imp <- variableImp(model, validation[[i]], newdata = v)

  r2 <- getR(model, validation[[i]], newdata = v)

  pred <- predict(model, v)

  return(list(imp=imp, r2 =r2, pred = pred))

}

#' Computing variable importance for each
#  variable in the dataset
#'
#' @description Extracting variable importance
#' and building a matrix. The importance of
#' predictors is reported in each column.
#' The name of the column is the name of
#' the predicted variables whereas the values
#' in each row correspond to the importance of
#' the variable. If I want to know the importance
#' of the 5th variable in predicting the
#' 1st i have to extract the value in the 1st column
#' at the 5th row.
#'
#' @param train
#' @param validation
#'
#' @return importance matrix
#' @export
#'
#' @examples
ImpMatrix <- function(train, validation) {

  library(foreach)
  library(doMC)

  old <- getDoParWorkers()
  registerDoMC(2)
  toPredict <- colnames(train)
  res <- foreach(v=toPredict) %dopar% {
    # Computing variable importance for each
    # variable in the dataset
    y <- train[,v]
    x <- train[,colnames(train) != v]

    model <- trainmod(x = x, y = y)

    # Get test data
    obs <- validation[,v]
    newdata <- validation[,colnames(validation) != v]

    # Get variable importance based on test data
    imp <- variableImp(model, obs, newdata, shuffles = 100)
    r2 <- getR(model, obs, newdata)
    list(imp=imp, r2 = r2)
  }
  registerDoMC(old)

  # Extracting variable importance
  # and building a matrix. The importance of
  # predictors is reported in each column.
  # The name of the column is the name of
  # the predicted variables whereas the values
  # in each row correspond to the importance of
  # the variable. If I want to know the importance
  # of the 5th variable in predicting the
  # 1st i have to extract the value in the 1st column
  # at the 5th row.


  cost <- sapply(seq_along(res), function(i) {
    r <- res[[i]]$imp[colnames(train)]
    names(r) <- colnames(train)
    r
  })
  colnames(cost) <- rownames(cost)
  return(cost)

}

#
# # Overall importance as predictors
# suggestion <- order(rowMeans(cost, na.rm = T), decreasing = T)
#
# sort(rowMeans(cost, na.rm = T))
#
# # Getting R-squared values
# r2 <- sapply(res, "[[", "r2")
# names(r2) <- colnames(cost)
#

#' Custom fitness function
#'
#' @param i
#' @param cost
#'
#' @return
#' @export
#'
#' @examples

simple_fitness <- function(i, cost){
  c <- cost[i,i]
  gain <- sum(c[upper.tri(c, diag = F)])
  loss <- sum(c[lower.tri(c, diag = F)])
  gain - loss
}

# fitness_custom <- function(i, cost, r2){
#   ordered <- data.frame(t(combn(colnames(cost)[i], 2)), stringsAsFactors = F)
#   keep <- apply(ordered, 1, function(x) cost[x[1], x[2]] > 0)
#   ordered <- ordered[keep,]
#
#   from <- unique(ordered[[1]])
#   to <- unique(ordered[[2]])
#
#   roots <- from[!from %in% to]
#
#   penalty <- sum(log(r2[!names(r2) %in% roots]))
#   c <- cost[i,i]
#
#   gain <- sum(c[upper.tri(c, diag = F)])
#   loss <- sum(c[lower.tri(c, diag = F)])
#   (gain - loss) + penalty
# }
#
# ga1 <- ga(fitness = simple_fitness, # custom fitness function
#           type = "permutation", # optimization data type
#           lower = 1,
#           upper = ncol(train),
#           cost = cost,
#           # r2 = r2,
#           monitor = plot,
#           suggestions = matrix(suggestion, nrow=1),
#           popSize = 10000,
#           maxiter = 100,
#           elitism = 1000,
#           keepBest = T,
#           parallel = T,
#           pmutation = 0.2,
#           run = 20)
#
# best <- data.frame(taxa=colnames(train)[ga1@solution[1,]],
#                    position=1:ncol(train))
#
# lapply(ga1@bestSol, function(x){
#   data.frame(taxa = colnames(train)[x],
#              position = 1:ncol(train))
# }) %>%
#   bind_rows() %>%
#   group_by(taxa, position) %>%
#   tally() %>%
#   ungroup() %>%
#   group_by(taxa) %>%
#   mutate(perc = n/sum(n),
#          ord = median(rep(position, n))) %>%
#   ungroup() %>%
#   mutate(taxa = reorder(taxa, ord, unique)) %>%
#   ggplot(aes(x = position, y = taxa, fill = perc)) +
#   geom_tile() +
#   geom_text(data = best, aes(x = position, y = taxa, label = "+"),
#             inherit.aes = F) +
#   scale_fill_distiller(palette = "Spectral") +
#   coord_cartesian(expand = F) +
#   theme(panel.background = element_rect(fill = "black"),
#         panel.grid = element_blank()) +
#   scale_x_continuous(breaks = 1:ncol(train))
#
#
# do.call(rbind, positions) %>%
#   group_by(taxa, position) %>%
#   tally() %>%
#   ungroup() %>%
#   group_by(taxa) %>%
#   mutate(perc = n/sum(n),
#          ord = mean(rep(position, n))) %>%
#   ungroup() %>%
#   mutate(taxa = reorder(taxa, ord, unique)) %>%
# ggplot(aes(x = position, y = taxa, fill = perc)) +
#   geom_tile() +
#   scale_fill_distiller(palette = "Spectral") +
#   coord_cartesian(expand = F) +
#   theme(panel.background = element_rect(fill = "black"),
#         panel.grid = element_blank()) +
#   scale_x_continuous(breaks = 1:ncol(train))
#
#
# do.call(rbind, positions) %>%
#   group_by(taxa, position) %>%
#   tally() %>%
#   ungroup() %>%
#   group_by(taxa) %>%
#   mutate(perc = n/sum(n),
#          ord = median(rep(position, n))) %>%
#   ungroup() %>%
#   mutate(taxa = reorder(taxa, ord, unique)) %>%
# ggplot(aes(x = position, y = perc)) +
#   geom_col() +
#   facet_grid(taxa ~ .)
#
# ga1 <- ga(fitness = function(i, costs) {
#   c <- costs[i,i]
#   gain <- sum(c[upper.tri(c, diag = F)])
#   loss <- sum(c[lower.tri(c, diag = F)])
#   gain - loss
# }, # custom fitness function
# type = "permutation", # optimization data type
# lower = 1,
# upper = ncol(train),
# costs = cost.mult,
# monitor = plot,
# popSize = 1000,
# maxiter = 100,
# elitism = 100,
# keepBest = T,
# parallel = T,
# pmutation = 0.2,
# run = 10)
#
#
#
# data.frame(taxa=colnames(train)[ga1@solution[1,]],
#            position=1:ncol(train))
#
# ordered <- colnames(train)[ga1@solution[1,]]
# ordered <- setNames(data.frame(t(combn(ordered, 2)), stringsAsFactors = F),
#                     c("from", "to"))
#
# colnames(cost) <- rownames(cost)
# names(r2) <- colnames(cost)
# apply(ordered, 1, function(x) {
#   penalty <- log(r2[x[2]])
#   c <- cost[x[1], x[2]]
#   c
#   # c + penalty
#   })
# ordered <- ordered[apply(ordered, 1, function(x) cost[x[1] ,x[2]]) > 0,]
#
# graph <- empty.graph(colnames(train))
# arcs(graph) <- ordered
#
# graphviz.plot(graph)


