setOldClass("train")
setClassUnion("trainOrNULL", c("train", "NULL"))

#' The S4 class for Classification and Regression information.
#'
#' @description This function generates a valid references to a class.
#'
#' @slot varobs data.frame. A data.frame containing observated values.
#' @slot yobs numeric. A numeric vector of observated values used for training
#' data.frame.
#' @slot classification train. Classification model.
#' @slot regression train. Regression models.
#'
#' @return
#' @export
#'
#' @examples

Modelling <- setClass("Modelling",
                      slots =c(varobs = "data.frame",
                               yobs = "numeric",
                               classification = "trainOrNULL",
                               regression = "trainOrNULL"))
setValidity("Modelling",
            function(object){
              if (nrow(object@varobs) == length(object@yobs))
                TRUE
              else "Error: y must be a numeric vector of length 'nrow(x)'"
            })


#' Fit and Tune Predective Models and Select the "best" one
#'
#' @description This function sets up a grid of tuning parametres for
#' classification and regression model. This function return an object of
#' \code{\link{Modelling-class}}.
#'
#' @param x  A data.frame where samples are in rows
#'  and features are in columns.
#'
#' @param y A numeric vector containing the observated outcome for each sample.
#'
#' @param cntr Control the computational nuances of the train function
#'  (for more detail see {\code{\link[caret]{trainControl}}}).
#' By default, is used the cross-validation method \code{"LOOCV"}.
#'
#' @param method A string specifying which classification and/or
#'  regression model to use (for details see \code{\link{caret}}.
#'  By default, is used the method \code{"rf"}.
#'  Possible values are found using names (\code{\link[caret]{getModelInfo}}).
#'
#' @param ... Arguments passed to the classification and regression
#'  algorithm (such as \code{\link{randomForest}})
#'
#' @return a Modelling S4 Class is returned;
#' for details see \code{\link{Modelling-class}}
#'
#' @import caret
#' @export
#'
#' @examples
#'

trainmod <- function(x, y,
                     cntr = caret::trainControl(method = "LOOCV"),
                     method = "rf", ...){
  varobs <- x
  yobs <- y

  y.fct <- factor(yobs > 0)

  classification <-  NULL

  if(nlevels(y.fct) > 1){
    classification <- caret::train(x = varobs, y = y.fct,
                                   method = method,
                                   trControl = cntr, ...)
    var <- varobs[yobs > 0,,drop = F]
    y <- yobs[yobs > 0]
  }

  regression <- caret::train(x = varobs, y = yobs,
                             method = method,
                             trControl = cntr, ...)

  res <- Modelling(regression = regression,
                   classification = classification,
                   varobs = varobs, yobs = yobs)
  return(res)
}

#' Classifciation and/or Regression Models Prediction
#'
#' @description These functions can be used for a single train object or
#' to loop through a number of \code{\link{Modelling-class}} objects to calculate
#' the training and test data predictions and class probabilities.
#' The newdata (an optional set of data to predict on - see
#' {\code{\link[stats]{predict}}}) is NULL.
#' If not NULL, newdata training data are used.
#' For more detail see {\code{\link[stats]{predict}}}
#'
#' @param Modelling For predict, an object of class \code{\link{Modelling-class}}.
#'
#' @return A numeric vector of predicted outcomes.
#'
#' @export
#'
#' @examples

setMethod("predict", "Modelling", function(object, newdata = NULL, ...) {

  if(!is.null(newdata)){

    newdata = as.data.frame(newdata)

    newdata <- newdata[, colnames(newdata) %in% object@regression$finalModel$xNames,
                       drop = F]

    pred <- predict(object@regression, newdata = newdata, ...)

    if(!is.null(object@classification)){

      newdata <- newdata[, colnames(newdata) %in% object@classification$finalModel$xNames,
                         drop = F]

      zeroes <- predict(object@classification, newdata = newdata, ...)

      pred[zeroes == F] <- 0
    }
    return(pred)
  }

  else {

    pred <- predict(object@regression, newdata = NULL, ...)

    if(!is.null(object@classification)){

      zeroes <- predict(object@classification, newdata = NULL, ...)

      pred[zeroes == F] <- 0
    }
    return(pred) }

  return(pred)
})

#' Measurement error models
#'
#' @description These functions are appropriate for cases where the model outcome is a number.
#' Given an object of \code{\link{Modelling-class}}, standard error (se),
#' mean absolute error (mae), Root Mean Squared Error (rmse), or R-squared (rsq)
#' is calculated.
#'
#' @param object An object of \code{\link{Modelling-class}}.
#' This method use a numeric vector of observated outcomes.
#'
#' @param obs  a numeric vector of observated outcomes.
#'
#' @param newdata An optional data-frame in which to look for variables
#' with which to predict. If omitted, the fitted values are used.
#'
#' @param type A string of character that defines a . Current possibilities
#' are "se", "mae, "rmse" and "rsq". Type must be specified
#'
#' @param ... Additional parameters to be passed the the S4 method
#'
#' @return
#' @export
#'
#' @examples

setGeneric("error", function(object, obs = NULL, newdata = NULL,
                             type = c("se", "mae", "rmse", "rsq"), ...)
  standardGeneric("error"))

setMethod("error", "Modelling",

          function(object, obs = NULL, newdata = NULL,
                   type = c("se", "mae","rmse", "rsq"), ...){

            match.arg(type)

            if(length(type) == 1){
              type %in% type}
            else stop("Error: 'type' must be of length 1")

            pred <- predict(object, newdata = newdata, ...)

            if(length(obs) == length(pred))
              TRUE
            else
              stop("Error: observation must ba a numeric vector of lenght prediction")

            if(is.null(obs) && is.null(newdata)){
              obs <- object@yobs }

            switch(type,
                   se = sum((obs - pred)^2),
                   mae = mean(abs(obs - pred)),
                   rmse = sqrt(sum((obs - pred)^2)/length(pred)),
                   rsq = (1 - (sum((obs - pred)^2))/sum((obs - mean(obs))^2))
            )
          }
)
