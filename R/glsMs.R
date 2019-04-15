#'GLS model selection for time series
#'
#'@param formula An object of class \code{formula} describing the model to be fitted. Should be
#'in the form \code{series ~ time}.
#'@param data Input data.frame to be analyzed.
#'@param diagnostic \code{logical}. Returns best fitting model from selection process if \code{TRUE}.
#'@param list_models \code{logical}. Returns coefficients and AICc scores of all fitted models.
#'
#'@param ... Other arguments may be passed to the gls.
#'
#'
#'@export
#'
#'@return gls model object
#'
#'
#'@examples
#'#Generate series
#'
#'m <- 0.1
#'x <- 1:30
#'
#'data <- data.frame(x = x,
#'                   y = m*x + rnorm(30, sd = 0.35))
#'glsMs(formula = y ~ x, data = data, diagnostic = T)

glsMs <- function(formula, data,
                  diagnostic = F,
                  list_models = F,...) {

  data <- model.frame(formula, data = data)
  names(data)[1] <- "y"
  names(data)[2] <- "x"

  if(nrow(data) < 30){
    message("Tests for trend are biased at low sample sizes.")
  }

  #Add quadratic term for quadratic models
  data$x2 <- data$x^2


  #Model fitting -------------------------------------------------------
  constant_norm <-
    nlme::gls(y ~ 1,
              data = data,
              na.action = na.omit,
              method = "ML", ...,
              ...)
  print("constant_norm")
  constant_ar1 <-
    try(nlme::gls(y ~ 1,
                  data = data,
                  correlation = nlme::corAR1(form = ~x),
                  na.action = na.omit,
                  method = "ML",
                  ...))
  print("constant_ar1")
  if (class(constant_ar1) == "try-error"){
    return(best_lm <- data.frame(model = NA,
                                 aicc  = NA,
                                 coefs..Intercept = NA,
                                 coefs.x = NA,
                                 coefs.x2 = NA,
                                 pval = NA))
  }

  #Null model with AR(2) error structure
  constant_ar2 <-
    try(nlme::gls(y ~ 1,
                  data = data,
                  correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                  method = "ML",
                  ...))
  print("constant_ar2")
  if (class(constant_ar2) == "try-error"){
    message("BFGS optimizer has failed, defaulting to Nelder-Mead routine (NULL AR2)")

    constant_ar2 <- try(nlme::gls(y ~ 1,
                                  data = data,
                                  control = nlme::glsControl(opt = "optim",
                                                             optimMethod = "Nelder-Mead"),
                                  correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                                  method = "ML", ...))
    if (class(constant_ar2) == "try-error"){

      message("Nelder-Mead optimizer has failed, defaulting to SANN routine (NULL AR2)")

      constant_ar2 <- try(nlme::gls(y ~ 1,
                                    data = data,
                                    control = nlme::glsControl(opt = "optim",
                                                               optimMethod = "SANN"),
                                    correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                                    method = "ML", ...))

      if (class(constant_ar2) == "try-error"){


        return(message("NUll AR2 model has failed!!!!!"))
      }

    }
  }


  # Linear model with normal error
  linear_norm <- nlme::gls(y ~ x,
                           data = data,
                           na.action = na.omit,
                           method = "ML", ...)

  # Linear model with AR1 error
  linear_ar1 <-
    try(nlme::gls(y ~ x,
                  data = data,
                  correlation = nlme::corAR1(form = ~x),
                  na.action = na.omit,
                  method = "ML", ...))

  if (class(linear_ar1) == "try-error"){
    return(best_lm <- data.frame(model = NA,
                                 aicc  = NA,
                                 coefs..Intercept = NA,
                                 coefs.x = NA,
                                 coefs.x2 = NA,
                                 pval = NA))

  }

  #Linear model with AR2 error
  linear_ar2 <- try(nlme::gls(y ~ x,
                              data = data,
                              correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                              method = "ML", ...))
  if (class(linear_ar2) == "try-error"){

    message("BFGS optimizer has failed, defaulting to Nelder-Mead routine (AR2)")

    linear_ar2 <- try(nlme::gls(y ~ x,
                                data = data,
                                control = nlme::glsControl(opt = "optim",
                                                           optimMethod = "Nelder-Mead"),
                                correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                                method = "ML", ...))
    if (class(linear_ar2) == "try-error"){

      message("Nelder-Mead optimizer has failed, defaulting to SANN routine (AR2)")

      linear_ar2 <- try(nlme::gls(y ~ x,
                                  data = data,
                                  control = nlme::glsControl(opt = "optim",
                                                             optimMethod = "SANN"),
                                  correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                                  method = "ML", ...))

      if (class(linear_ar2) == "try-error"){


        return(message("Linear AR2 model has failed!!!!!"))
      }
    }
  }

  # Polynomial model with normal error
  poly_norm <- nlme::gls(y ~ x + x2,
                         data = data,
                         na.action = na.omit,
                         method = "ML", ...)

  # Polynomial model with AR1 error
  poly_ar1 <-
    try(nlme::gls(y ~ x + x2,
                  data = data,
                  correlation = nlme::corAR1(form = ~x),
                  na.action = na.omit,
                  method = "ML", ...))
  if (class(poly_ar1) == "try-error"){
    return(best_lm <- data.frame(model = NA,
                                 aicc  = NA,
                                 coefs..Intercept = NA,
                                 coefs.x = NA,
                                 coefs.x2 = NA,
                                 pval = NA))
  }

  #Polynomial model with AR2 error
  poly_ar2 <- try(nlme::gls(y ~ x + x2,
                            data = data,
                            correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                            method = "ML", ...))
  if (class(linear_ar2) == "try-error"){

    message("BFGS optimizer has failed, defaulting to Nelder-Mead routine (AR2)")

    linear_ar2 <- try(nlme::gls(y ~ x + x2,
                                data = data,
                                control = nlme::glsControl(opt = "optim",
                                                           optimMethod = "Nelder-Mead"),
                                correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                                method = "ML", ...))
    if (class(linear_ar2) == "try-error"){

      message("Nelder-Mead optimizer has failed, defaulting to SANN routine (AR2)")

      linear_ar2 <- try(nlme::gls(y ~ x + x2,
                                  data = data,
                                  control = nlme::glsControl(opt = "optim",
                                                             optimMethod = "SANN"),
                                  correlation = nlme::corARMA(form = ~x, p = 2, q = 0),
                                  method = "ML", ...))

      if (class(linear_ar2) == "try-error"){


        return(message("Quadratic AR2 model has failed!!!!!"))
      }
    }
  }


  print("here")
  # Calculate AICs for all models
  df_aicc <-
    data.frame(model = c("poly_norm",
                         "poly_ar1",
                         "poly_ar2",

                         "linear_norm",
                         "linear_ar1",
                         "linear_ar2"),
               aicc  = c(AICcmodavg::AICc(poly_norm),
                         AICcmodavg::AICc(poly_ar1),
                         AICcmodavg::AICc(poly_ar2),

                         AICcmodavg::AICc(linear_norm),
                         AICcmodavg::AICc(linear_ar1),
                         AICcmodavg::AICc(linear_ar2)),

               coefs = rbind(coef(poly_norm),
                             coef(poly_ar1),
                             coef(poly_ar2),
                             c(coef(linear_norm), NA),
                             c(coef(linear_ar1),  NA),
                             c(coef(linear_ar2), NA)),

               pval = c(anova(update(constant_norm),
                              update(poly_norm))$`p-value`[2],

                        anova(update(constant_ar1),
                              update(poly_ar1))$`p-value`[2],

                        anova(update(constant_ar1),
                              update(linear_ar1))$`p-value`[2],

                        anova(update(constant_norm),
                              update(linear_norm))$`p-value`[2],

                        anova(update(constant_ar1),
                              update(linear_ar1))$`p-value`[2],

                        anova(update(constant_ar2),
                              update(linear_ar2))$`p-value`[2]),

               stringsAsFactors = FALSE)

  df_aicc <- dplyr::arrange(df_aicc, aicc)


  #Find best model
  best_lm <- df_aicc[df_aicc$aicc == min(df_aicc$aicc),]

  if (best_lm$model == "poly_norm") {
    model <- poly_norm
  } else if (best_lm$model == "poly_ar1") {
    model <- poly_ar1
  } else if (best_lm$model == "poly_ar2") {
    model <- poly_ar2
  } else if (best_lm$model == "linear_norm") {
    model <- linear_norm
  } else if (best_lm$model == "linear_ar1") {
    model <- linear_ar1
  } else if (best_lm$model == "linear_ar2") {
    model <- linear_ar2
  }

  model <- c(model, list("glsMs" = best_lm$model))
  class(model) <- "gls"

  if (!diagnostic){
    return(model)
  } else {
    return(list(model, df_aicc))
  }


}
