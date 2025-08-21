#' @importFrom rlang .data

utils::globalVariables(c("numHyp", "alphaTotal",
                         "pdigits", "idigits", "plotInPercent", "Tmax_atPlot",
                         "mtParam", "enrollmentAll",
                         "inputD", "G"))

# operations with IA timing and IF ----

# derive analysis time since start using input dataset 'D' and 'enrollment'
deriveCalendarTime <- function(D, enrollment) {
  Ts <- vector("list", nrow(D))
  if (!"enrollment" %in% names(D))
    D <-
      D %>% tibble::add_column(enrollment = rep(list(enrollment), nrow(D)))
  for (i in 1:nrow(D)) {
    # use 'iaTime' if specified explicitly
    if (!is.null(D$iaTime[[i]])) {
      if (!any(is.na(D$iaTime[[i]]))) {
        Ts[[i]] <- D$iaTime[[i]]
        next
      }
    }
    # else derive
    Ts[[i]] <- sapply(D$iaSpec[[i]], function(e) {
      n2Time(
        D$endpointParam[[e$H]],
        n = e$atIF * D$hypN[e$H],
        D$enrollment[[e$H]],
        ratio = D$allocRatio[e$H]
      )
    })
  }
  return(Ts)
}

# derivation sample size from calendar time at IA
derive_Ns <- function(D, enrollment, doRounding = TRUE) {
  Ns <- vector("list", nrow(D))
  if (!"enrollment" %in% names(D))
    D <-
      D %>% tibble::add_column(enrollment = rep(list(enrollment), nrow(D)))
  for (i in 1:nrow(D)) {
    # if D$infoFr and D$hypN are given explicitly use them
    # to set Ns otherwise derive using D$iaTime and D$endpointParam
    if (!is.null(D$infoFr[[i]]) && !is.null(D$hypN[[i]]))
      if (!any(is.na(D$infoFr[[i]])) && !is.na(D$hypN[[i]])) {
        Ns[[i]] <- D$infoFr[[i]] * D$hypN[[i]]
        next
      }
    # else derive Ns
    Ns[[i]] <- Time2n(
      x = D$endpointParam[[i]],
      T = D$iaTime[[i]],
      enrollment = D$enrollment[[i]],
      ratio = D$allocRatio[[i]]
    )
  }
  if (doRounding)
    Ns <- lapply(Ns, function(x)
      round(x))
  return(Ns)
}

# derivation of information fractions    # currently not used
deriveIF <- function(D, enrollment, digits = 2) {
  IFs <- list(nrow(D))
  for (i in 1:nrow(D)) {
    IFs[[i]] <- round(
      Time2n(
        x = D$endpointParam[[i]],
        T = D$iaTime[[i]],
        enrollment = enrollment,
        ratio = D$allocRatio[[i]]
      ) /
        D$hypN[[i]],
      digits = digits
    )
  }
  return(IFs)
}

# Obtain the parameter under H1 on the natural scale ----

#' appendMCP internal S3 
#' @description Obtain the parameter under H1 on the natural scale
#' @param x A endpoint object
#' @keywords internal
#' @name getDelta
#' @rdname internalS3
getDelta <- function(x){
  UseMethod("getDelta", x)
}

#' @export
#' @method getDelta default
#' @keywords internal
#' @name getDelta.default
#' @rdname internalS3
getDelta.default <- function(x) {
  x$p1 - x$p2
}

#' @export
#' @method getDelta tte_exp
#' @keywords internal
#' @name getDelta.tte_exp
#' @rdname internalS3
getDelta.tte_exp <- function(x) {
  x$p1 / x$p2     # return hazard ratio
}

#' @export
#' @method getDelta tte_pwe
#' @keywords internal
#' @name getDelta.tte_pwe
#' @rdname internalS3
getDelta.tte_pwe <- function(x) {
  utils::head(x$p1/x$p2, 1)     # return hazard ratio
}

# Format details of the effect size (parameter on the natural scale) ----

#' @description Format a string reporting the effect
#' @param x An endpoint object
#' @keywords internal
#' @name getEffectSizeDetails
#' @rdname internalS3
getEffectSizeDetails <- function (x) {
  UseMethod("getEffectSizeDetails", x)
}

#' @method getEffectSizeDetails normal
#' @keywords internal
#' @name getEffectSizeDetails.normal
#' @rdname internalS3
#' @export
getEffectSizeDetails.normal <- function(x) {
  sprintf("%2.f", x$p1 - x$p2)
}

#' @method getEffectSizeDetails binomial
#' @keywords internal
#' @name getEffectSizeDetails.binomial
#' @rdname internalS3
#' @export
getEffectSizeDetails.binomial <- function(x) {
  sprintf("%.2f (%d%% vs %d%%)",
          x$p1 - x$p2,
          round(x$p1 * 100),
          round(x$p2 * 100))
}

#' @method getEffectSizeDetails binomial_pooled
#' @keywords internal
#' @name getEffectSizeDetails.binomial_pooled
#' @rdname internalS3
#' @export
getEffectSizeDetails.binomial_pooled <- function(x) {
  sprintf("%.2f (%d%% vs %d%%)",
          x$p1 - x$p2,
          round(x$p1 * 100),
          round(x$p2 * 100))
}

#' @method getEffectSizeDetails binomial_unpooled
#' @keywords internal
#' @name getEffectSizeDetails.binomial_unpooled
#' @rdname internalS3
#' @export
getEffectSizeDetails.binomial_unpooled <- function(x) {
  sprintf("%.2f (%d%% vs %d%%)",
          x$p1 - x$p2,
          round(x$p1 * 100),
          round(x$p2 * 100))
}

#' @method getEffectSizeDetails tte_exp
#' @keywords internal
#' @name getEffectSizeDetails.tte_exp
#' @rdname internalS3
#' @export
getEffectSizeDetails.tte_exp <- function(x) {
  sprintf("HR = %.2f (mCntl = %.1f mo)", x$p1 / x$p2,-log(1 / 2) / x$p2)
}

# Median of a PWE distribution
med_pwe <- function(duration, rate) {
  stats::uniroot(
    f        = \(x) gsDesign2::ppwe(x, duration, rate) - 0.5,
    interval = c(1e-6, 1e6)
  )$root
}

#' @method getEffectSizeDetails tte_pwe
#' @keywords internal
#' @name getEffectSizeDetails.tte_pwe
#' @rdname internalS3
#' @export
getEffectSizeDetails.tte_pwe <- function(x) {
  sprintf("HR = %.2f (mCntl = %.1f mo)", utils::head(x$p1 / x$p2, 1),
          med_pwe(x$durations, x$p2))
}

# compute standardization factor ----

# compute standardization factor, psi,  that link natural parameter to the effect size, theta
# that non-centrality parameter is in the form theta*sqrt(n)= psi*delta*sgrt(n)
# e.g, for log-rank test, psi = sqrt(c*(1-c)), where c=1/(1+ratio), ratio is the allocation ratio

#' @description
#'  Compute standardization factor, psi, that links
#' natural parameter to the effect size, theta
#'
#' @param x An endpoint object
#' @param ratio Allocation ratio
#' @keywords internal
#' @name getStandardizingCoef
#' @rdname internalS3
#' 
#' @examplesIf FALSE
#' # Internally dispatched examples (for development or internal tests)
#' getStandardizingCoef(list(list(), class = "tte_exp"), ratio = 1)
#' # should return 0.5
getStandardizingCoef <- function (x, ratio=1) {
  UseMethod("getStandardizingCoef", x )
}

#' @method getStandardizingCoef default
#' @keywords internal
#' @name getStandardizingCoef.default
#' @rdname internalS3
#' @export
getStandardizingCoef.default <- function(x, ratio = 1){
  sqrt(ratio)/(1+ratio) * 1
}

#' @method getStandardizingCoef binomial
#' @keywords internal
#' @name getStandardizingCoef.binomial
#' @rdname internalS3
#' @export
getStandardizingCoef.binomial <- function(x, ratio = 1){
  # x$alpha and x$beta use used to calibrate standardization coefficient
  # when mapping binomial rates to normal effect size using n.I
  if (is.null(x$alpha) | is.null(x$beta)){
    x$alpha <- 0.025
    x$beta  <- 0.85
  }
  nAux <- gsDesign::nBinomial( p1=x$p1, p2=x$p2, ratio=ratio,alpha=x$alpha, beta=x$beta)
  thetaAux <- (stats::qnorm(1-x$alpha) + stats::qnorm(1-x$beta))/sqrt(nAux)
  return(thetaAux/(x$p1-x$p2))
}

# binomial_pooled gives somewhat conservative power
#' @method getStandardizingCoef binomial_pooled
#' @keywords internal
#' @name getStandardizingCoef.binomial_pooled
#' @rdname internalS3
#' @export
getStandardizingCoef.binomial_pooled <- function(x, ratio = 1){
  if(!is.array(x$pPooled)) pPooled <- x$pPooled
  if (!is.null(x$p1) && !is.null(x$p2)) pPooled <- 1/(1+ratio)*(x$p1+ratio*x$p2)
  sqrt(ratio)/(1+ratio) / sqrt( pPooled *(1-pPooled))
}

#' @method getStandardizingCoef binomial_unpooled
#' @keywords internal
#' @name getStandardizingCoef.binomial_unpooled
#' @rdname internalS3
#' @export
getStandardizingCoef.binomial_unpooled <- function(x, ratio = 1){
  1/sqrt(x$p1*(1-x$p1)*(1+ratio) + x$p2*(1-x$p2)*(1+ratio)/ratio)
}


# operation with number and timing of events ----

# combine events from control and treatment arms
eEvents_total <-
  function(hr = 1,
           ratio = 1,
           lambdaC = 1,
           eta = 0,
           gamma = 1,
           R = 1,
           S = NULL,
           T,
           Tfinal = NULL,
           minfup = 0,
           digits = 4,
           target = 0)
# assuming a common dropout
  {
    if (T == 0)
      return(0)
    Qe <- ratio / (1 + ratio)
    eDC <-
      gsDesign::eEvents(
        lambda = lambdaC,
        eta = eta,
        gamma = gamma * (1 - Qe),
        R = R,
        S = S,
        T = T,
        Tfinal = Tfinal,
        minfup = minfup
      )
    eDE <-
      gsDesign::eEvents(
        lambda = lambdaC * hr,
        eta = eta,
        gamma = gamma * Qe,
        R = R,
        S = S,
        T = T,
        Tfinal = Tfinal,
        minfup = minfup
      )
    return(sum(eDC$d + eDE$d) - target)
  }
eEvents_totalVec <- Vectorize(eEvents_total, c("hr", "T"))

# get time when a given number of events will be reached
tEvents <-
  function (n,
            hr = 1,
            ratio = 1,
            lambdaC = 1,
            eta = 0,
            gamma = 1,
            R = 1,
            S = NULL,
            T,
            Tfinal = NULL,
            minfup = 0,
            tol = .Machine$double.eps ^ 0.25)
# n is the target number of events
  {
    # solve eEvents_total for parameter T
    z <-
      stats::uniroot(
        f = eEvents_total,
        interval = c(1e-04, 1e5),
        hr = hr,
        ratio = ratio,
        lambdaC = lambdaC,
        eta = eta,
        gamma = gamma,
        R = R,
        S = S,
        Tfinal = Tfinal,
        minfup = minfup,
        target = n,
        tol = tol
      )
    z$root
  }
tEventsVec <- Vectorize(tEvents, c("hr", "n"))

#' @description
#' Obtain a calendar time given the sample size (or events)
#' @param x An endpoint object
#' @param n A sample size or number of events
#' @param enrollment A enrollement object
#' @param ratio An allocation ration
#' @keywords internal
#' @rdname internalS3
#' @name n2Time
n2Time <- function (x, n, enrollment, ratio) {
  UseMethod("n2Time", x)
}

#' @method n2Time default
#' @keywords internal
#' @name n2Time.default
#' @rdname internalS3
#' @export
n2Time.default <- function(x, n, enrollment, ratio = 1) {
  fun <- function(t) {
    Time2n(x, t, enrollment = enrollment, ratio = ratio) - n
  }
  if (is.null(x$maturityTime))
    x$maturityTime <- 0
  # accrualDuration <- group_by(enrollment, stratum) %>%
  #     dplyr::mutate( totDuration = sum(duration))  %>% ungroup() %>%
  #     dplyr::select(totDuration) %>% max()
  accrualDuration <- sum(enrollment$duration)
  stats::uniroot(fun, c(0, accrualDuration + x$maturityTime))$root
}

#' @method n2Time tte_exp
#' @keywords internal
#' @name n2Time.tte_exp
#' @rdname internalS3
#' @export
n2Time.tte_exp <- function(x, n, enrollment, ratio = 1) {
  tEvents(
    n = n,
    hr = x$p1 / x$p2,
    ratio = ratio,
    lambdaC = x$p2,
    eta = x$dropoutHazard,
    gamma = enrollment$rate,
    R = enrollment$duration
  )
}

#' @method n2Time tte_pwe
#' @keywords internal
#' @name n2Time.tte_pwe
#' @rdname internalS3
#' @export
n2Time.tte_pwe <- function(x, n, enrollment, ratio = 1) {
  tEvents(
    n = n,
    hr = utils::head(x$p1 / x$p2, 1),
    ratio = ratio,
    lambdaC = x$p2,
    eta = x$dropoutHazard,
    gamma = enrollment$rate,
    R = enrollment$duration,
    S = x$durations,
  )
}

#' @description
#'  Calculate the sample size (or the number of events) available
#'  at the given calendar time
#'
#' @param x An endpoint object
#' @param T Calendar time points
#' @param enrollment Enrollment object
#' @param ratio Allocation ratio
#' @keywords internal
#' @rdname internalS3
#' @name Time2n
Time2n <- function (x, T, enrollment, ratio) {
  UseMethod("Time2n", x)
}

#' @method Time2n default
#' @keywords internal
#' @name Time2n.default
#' @rdname internalS3
#' @export
Time2n.default <- function(x, T, enrollment, ratio = 1) {
  if (is.null(x$dropoutHazard))
    eta <- 0
  else
    eta <- x$dropoutHazard
  if (is.null(x$maturityTime))
    maturatyTime <- 0
  else
    maturatyTime <- x$maturityTime

  timeAux <- c(0, cumsum(enrollment$duration))
  N_rand <- cumsum(c(0, enrollment$rate * c(enrollment$duration)))

  stats::approx(
    x = timeAux + maturatyTime,
    y = N_rand * exp(as.numeric(-eta) * (maturatyTime)),
    rule = 2,
    xout = T
  )$y
}

#' @method Time2n tte_exp
#' @keywords internal
#' @name Time2n.tte_exp
#' @rdname internalS3
#' @export
Time2n.tte_exp <- function(x, T, enrollment, ratio = 1) {
  eEvents_totalVec(
    T = T,
    hr = x$p1 / x$p2,
    ratio = ratio,
    lambdaC = x$p2,
    eta = x$dropoutHazard,
    gamma = enrollment$rate,
    R = enrollment$duration
  )
}

#' @method Time2n tte_pwe
#' @keywords internal
#' @name Time2n.tte_pwe
#' @rdname internalS3
#' @export
Time2n.tte_pwe <- function(x, T, enrollment, ratio = 1) {
  eEvents_totalVec(
    T = T,
    hr = utils::head(x$p1 / x$p2, 1),
    ratio = ratio,
    lambdaC = x$p2,
    eta = x$dropoutHazard,
    gamma = enrollment$rate,
    R = enrollment$duration,
    S = x$durations
  )
}

# utils functions ----

#' Enrollment and Data Availability Plot
#'
#' @param D A dataset containing interim analysis information
#' @param enrollment An enrollment object
#' @param Tmax A scalar defining the maximum time on the x-axis
#' @param plotInPercent logical; if TRUE the y-axis reports in percentages
#'
#' @returns A ggplot object representing data availability for each hypothesis
#' @export
#'
#' @examples
#' \dontrun{
#' plot_iaTiming(D,enrollment)
#' }
#' @keywords internal
plot_iaTiming <-
  function(D,
           enrollment,
           Tmax = 80,
           plotInPercent = TRUE) {
    if (!"enrollment" %in% names(D))
      D <-
        D %>% tibble::add_column(enrollment = rep(list(enrollment), nrow(D)))


    timeGrid <- seq(0, Tmax, 1)
    Nmax <- sum(enrollment$rate * enrollment$duration)
    N_Rand <- stats::approx(
      x = c(0, cumsum(enrollment$duration)),
      y = c(0, cumsum(enrollment$rate * enrollment$duration)),
      xout = timeGrid,
      rule = 2
    )$y
    Randomized <- if (plotInPercent)
      N_Rand / Nmax * 100
    else
      N_Rand

    Y <- list(nrow(D))
    for (i in 1:nrow(D)) {
      N <-
        Time2n(
          x = D$endpointParam[[i]],
          T = timeGrid,
          enrollment = D$enrollment[[i]],
          ratio = D$allocRatio[[i]]
        )
      Y[[i]] <- if (plotInPercent)
        N / D$hypN[i] * 100
      else
        N
    }


    maturity <- sapply(D$endpointParam,
                       function(x)
                         ifelse(
                           !is.null(x$maturityTime),
                           paste0(", ", x$maturityTime, " months data"),
                           ""
                         ))
    # names(Y) <- paste0( D$ep, " (",D$id,")", maturity) # to be a plot legend
    names(Y) <-
      paste0(D$id, ": ", D$tag, " ", maturity) # to be a plot legend
    dat <- tibble::tibble(
      Time = timeGrid,
      Randomized,
      as.data.frame(Y, check.names = FALSE),
      .name_repair = "minimal"
    ) %>%
      tidyr::pivot_longer(!.data$Time, names_to = "Type", values_to = "Y") %>%
      dplyr::mutate(Type = factor(.data$Type, levels = c(
        "Randomized",
        setdiff(unique(.data$Type), "Randomized")
      )))

    if (plotInPercent)
      dat <- dat %>% dplyr::filter(Y <= 105)

    ylabStr <-
      ifelse(plotInPercent,
             "Statistical Information, in %",
             "Statistical Information")
    iaTimes <- unique(unlist(D$iaTime))
    ans <-
      ggplot2::ggplot(data = dat, ggplot2::aes(x = .data$Time, y = .data$Y,
                                               col = .data$Type)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(xintercept = iaTimes, linetype = "dashed") +
      ggplot2::labs(x = "Time, in months",
           y = ylabStr) +
      ggplot2::scale_x_continuous(breaks = seq(0, Tmax, 6), limits = c(0, Tmax)) +
      ggplot2::annotate(
        "text",
        x = iaTimes,
        y = 0,
        label = round(iaTimes, 0),
        angle = 45,
        size = 3.5,
        hjust = 0
      )
    if (plotInPercent)
      ans <- ans + ggplot2::scale_y_continuous(breaks = seq(0, 100, 10))

    ans
  }


# aux function that put spending function info from list to char string
print_sfInfo <- function(x, digits = 4)
{
  if (!is.null(x)) {
    if (is.character(x$sfu)) {
      sfName <- dplyr::case_match(
        x$sfu,
        "OF" ~ "O'Brien - Fleming",
        "WT" ~ " Wang - Tsiatis",
        "Pocock" ~ "Pocock",
        .default  = "mis-specified"
      )
      if (sfName == "mis-specified")
        stop("Mis-specified boundary")
    } else
      sfName <- get_sf_name(x$sfu)
    #do.call(x$sfu, list(0.025, 1, x$sfupar))$name # dummy call

    if (length(x$nominal)>0 )
      sfName <- paste0("Nominal spend at IA ",
                       paste(seq_along(x$nominal),collapse = ", "),
                       ", then ", sfName
      )
    if (!is.null(x$sfupar) && length(x$sfupar) <= 4) {
      sfParam <- paste0(", parameter = ",
                        paste(lapply(x$sfupar, function(y) {
                          if (is.numeric(y))
                            round(y, digits)
                          else
                            y
                        }), collapse = " "))
    } else
      sfParam <- ""
    paste0(sfName, sfParam)
  } else {
    "No group sequential testing"
  }
}

# given graph G, for each hypothesis generate possible weights and list what
# rejections lead to that weight
getPossibleWeightsInfo <- function(G, # G is gMCP graph object
                                   numDigitsToRound = 5,
                                   and_sympbol = ", ") {
  # aux function to form a character string of reject Hj
  formRejStr <- function (x) {
    aux <- sapply(x, function(y) {
      if (all(y == 1))
        "Initial allocation"
      else
        paste0("Successful ",
               paste("H",
                     which(y == 0),
                     collapse = and_sympbol, #", ",
                     sep = ""))
    })
    paste(aux, collapse = " or ", sep = "_")
  }  # formRejStr

  K <- length(gMCPLite::getNodes(G)) # number of hypotheses in G
  W <- gMCPLite::generateWeights(G)
  # For each Hj generate possible weights and info of rejection scenarios that lead to that weight
  res <-
    tibble::tibble(Hint = apply(W[, 1:K], 1, c, simplify = FALSE),
                   as.data.frame(W[,-c(1:K)])) %>%
    dplyr::mutate(# number of elementary hyp in 'Hint'
      mj = purrr::map_dbl(.data$Hint, sum)) %>%
    tidyr::pivot_longer(-c(.data$Hint, .data$mj), names_to = "Hj", values_to = "possibleWeight") %>%
    dplyr::mutate(possibleWeight = round(.data$possibleWeight, numDigitsToRound)) %>%
    # dplyr::filter(possibleWeight>0) %>%  ## ?? YT Sep. 2023
    dplyr::group_by(.data$Hj, .data$possibleWeight) %>%
    dplyr::reframe(max_mj = max(.data$mj),
                       mj = .data$mj,
                     Hint = .data$Hint) %>%
    dplyr::group_by(.data$Hj, .data$possibleWeight) %>%
    dplyr::filter(.data$mj == .data$max_mj) %>%
    dplyr::reframe(listHj = list(.data$Hint)) %>%
    dplyr::group_by(.data$Hj) %>%
    dplyr::mutate(rejectedHypInfo = purrr::map_chr(.data$listHj, formRejStr)) %>%
    dplyr::select("Hj", "possibleWeight", "rejectedHypInfo") %>%
    dplyr::mutate(Hj = as.numeric(gsub(".*?([0-9]+).*", "\\1", .data$Hj)))
  return(res)
}

# core function that calculate bounds given graph and GSD info
report_MT_grSeq <- function(
    G,        # graph for H1, ..., Hm
    D,
    # infoFr,   # list of length m setting info fraction for the Hi
    # spendFun,  # list of spending functions for Hi
    sigLevel  = 0.025, # overall significance level usually (1-sided)
    pdigits   = 5,     # p-value digits
    ddigits   = 3,     # delta digits
    idigits   = 3,     # info fraction digits
    powdigits = 2      # power digits
)
{
  nodes <- gMCPLite::getNodes(G)
  m <- length(nodes)

  # possibleWeight <- apply(generateWeights(G)[,-c(1:m)],2,unique,simplify = FALSE)
  wInfo <- getPossibleWeightsInfo(G)
  possibleWeight_aux <- split(wInfo$possibleWeight, wInfo$Hj)
  scenarioInfo   <-
    wInfo %>%
    dplyr::group_by(.data$Hj, .data$possibleWeight) %>%
    #    summarise(wInfo = rejectedHypInfo) %>% ?? why I had this line
    dplyr::mutate(hypNames = nodes[.data$Hj])

  get_sfu <- function(sfu, nominal) {
    if ( length(nominal) > 0) {
      gsDesign::sfPoints
    } else {
      sfu
    }
  }

  get_sfupar <- function(timing, sfu, sfupar, nominal, alpha) {
    if (length(nominal)>0) {
      c(cumsum(nominal), sfu(alpha, timing[-seq_along(nominal)], sfupar)$spend)/alpha
    } else {
      sfupar
    }
  }

  paramData <- tibble::tibble(
    hypNames        = nodes,
    test.type       = 1,
    k               = sapply(D$infoFr, length),
    timing          = D$infoFr,
    sfu             = sapply(D$grSeqTesting, function(x) x$sfu),
    sfupar          = sapply(D$grSeqTesting, function(x) x$sfupar),
    nominal         = sapply(D$grSeqTesting, \(x) x$nominal),
    sfInfo          = sapply(D$grSeqTesting, print_sfInfo),
    possibleWeight  = possibleWeight_aux
  ) %>% tidyr::unnest(.data$possibleWeight) %>%
    dplyr::mutate(
      alpha = .data$possibleWeight * sigLevel,
      weight_alpha = paste0(.data$possibleWeight, " (", .data$alpha, ")")
    ) %>%
    dplyr::filter(.data$alpha > 0) |>
    dplyr::mutate(sfupar = purrr::pmap(.l = list(.data$timing, .data$sfu, .data$sfupar, .data$nominal, .data$alpha),
                                       .f = get_sfupar),
                  sfu    = purrr::map2(.data$sfu, .data$nominal, get_sfu)) |>
    dplyr::mutate(spend = purrr::pmap(.l = list(sfu = .data$sfu, sfupar = .data$sfupar,
                                                timing = .data$timing, alpha = .data$alpha),
                                      .f = function(sfu, sfupar, timing, alpha) {
                                        if (is.null(sfu)) {
                                          NULL
                                        } else {
                                          sfu(alpha, timing, sfupar)$spend
                                        }
                                      }))

  gsDesign_m <- function(...) {
    args <- list(...)
    if (args$k == 1) {
      # if a single analysis spend all alpha at once
      return(list(upper = list(
        bound = stats::qnorm(args$alpha, lower.tail = FALSE)
      )))
    } else{
      sf_name <- try(get_sf_name(args$sfu), silent = TRUE)
      if (sf_name == "Bonferroni adjustment")
        return(list(upper = list(
          bound = stats::qnorm(args$alpha * diff(c(0,args$sfupar)),lower.tail = FALSE)
        )))
      return(do.call(gsDesign::gsDesign, args))
    }
  }
  gsDesignArgs <-
    paramData %>% dplyr::select("alpha", "test.type", "k", "timing", "sfu", "sfupar")
  gsDesignList <- purrr::pmap(gsDesignArgs, gsDesign_m)
  # extract nominal p-values
  nominalPval <-
    lapply(gsDesignList, function(x)
      stats::pnorm(x$upper$bound, lower.tail = FALSE))

  res <-
    paramData %>% tibble::add_column(nominalPval) %>%
    dplyr::left_join(scenarioInfo)
  res$Analysis <-
    sapply(res$k, function(x)
      paste(1:x, collapse = "<br>"))
  res$timing_vec <-
    sapply(res$timing, function(x)
      paste(round(x, digits = idigits), collapse = "<br>"))
  res$nominalPvalx2 <- sapply(res$nominalPval, function(x)
    paste(round( 2*x, digits = pdigits), collapse = "<br>") )
  res$nominalPval   <-
    sapply(res$nominalPval, function(x)
      paste(round(x, digits = pdigits), collapse = "<br>"))
  res$scenarioInfo <- dplyr::filter(scenarioInfo, .data$possibleWeight > 0)
  # run power evaluation if given effect size and sample size
  if (!is.null(D$theta) && !is.null(D$hypN)) {
    # extract nominal p-values
    aux <- tibble::tibble(
      hypNames = paramData$hypNames,
      k = sapply(paramData$timing, length),
      timing =    paramData$timing,
      a = lapply(gsDesignList, function(x) {
        if (!is.null(x$lower$bound))
          x$lower$bound
        else
          rep(-20, length(x$upper$bound))
      }),
      b = lapply(gsDesignList, function(x)
        x$upper$bound)
    ) %>%
      dplyr::left_join(dplyr::select(D, "hypNames", "theta", "hypN", "standFactor", "logDelta"))
    aux$n.I <-
      mapply(FUN = "*", aux$timing, aux$hypN, SIMPLIFY = FALSE)

    aux$thetaHat <-
      mapply(
        FUN = function(z, n)
          z / sqrt(n),
        aux$b,
        aux$n.I,
        SIMPLIFY = FALSE
      )
    deltaHat <- mapply(function(th, c, e) {
      ans <- th / c
      if (e)
        exp(-ans)
      else
        ans
    },
    aux$thetaHat,
    aux$standFactor,
    aux$logDelta,
    SIMPLIFY = FALSE)

    gsProbabilityArgs <-
      aux %>%  dplyr::select("k", "theta", "n.I", "a", "b")
    gsProbabilityList <- purrr::pmap(gsProbabilityArgs, gsDesign::gsProbability)
    pow <-
      lapply(gsProbabilityList, function(x)
        cumsum(x$upper$prob))

    res <- res %>% tibble::add_column(deltaHat)
    res$deltaHat <-
      sapply(deltaHat, function(x)
        paste(round(x, digits = ddigits), collapse = "<br>"))

    res <- res %>% tibble::add_column(pow)
    res$pow <-
      sapply(res$pow, function(x)
        paste(round(x, digits = powdigits), collapse = "<br>"))

  }
  res %>% dplyr::arrange(as.numeric(gsub(".*?([0-9]+).*", "\\1", .data$hypNames)), .data$alpha)
}

# knit functions for tables ----

#' Local alpha scenarios
#'
#' @param hyp_testing_dataset A dataset defining testing hypotheses
#' @param digits A scalar defining the number of digits to report in a table
#'
#' @returns A table that reports all possible scenarios for the local significance level
#' @export
#'
#' @examples
#' \dontrun{
#' knit_MT_table(hyp_testing_dataset)
#' }
#' @keywords internal
knit_MT_table <- function(hyp_testing_dataset, digits = 5) {
  df <- hyp_testing_dataset %>%
    dplyr::select("hypNames", "alpha", "possibleWeight", "rejectedHypInfo") %>%
    dplyr::rename(
      'Local alpha level' = .data$alpha,
      'Weight'            = .data$possibleWeight,
      'Testing Scenario'  = .data$rejectedHypInfo,
    )

  if (knitr::is_html_output()) {
    df[, -1] %>%
      knitr::kable("html", escape  = F, digits  = digits,
            caption = "List of possible local alpha levels following the graphical testing procedure") %>%
      kableExtra::kable_styling() %>%
      kableExtra::pack_rows( index         = table(forcats::fct_inorder(df$hypNames)),
                 label_row_css = "text-align: left;"
      )
    # %>%
    # column_spec(1, latex_valign = "m")
    # collapse_rows(columns = 1, valign = "top")
  } else if (knitr::is_latex_output()) {
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    df[, -1] %>%
      dplyr::mutate_all(kableExtra::linebreak) %>%
      knitr::kable(
        "latex",
        booktabs = T,
        escape = F,
        longtable = TRUE,
        caption = "Efficacy p-value Boundaries"
      ) %>%
      kableExtra::kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
      kableExtra::pack_rows(index = table(forcats::fct_inorder(df$hypNames)))
  } else if (knitr::pandoc_to("docx")) {
    #require(flextable)
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    flextable::flextable(df)
  }
}

# prepare table the cut-off nominal p-vals for all Hi at all analyses

#' Nominal p-value boundary scenarios
#'
#' @param hyp_testing_dataset A dataset defining testing hypotheses
#' @param digits A scalar defining the number of digits to report in a table
#' @param include_nominalPvalx2 logical; if TRUE a column with 2-sided p-values is reported
#'
#' @returns A table reporting boundaries to compare with observed
#' p-values calculated for the test statistics at the corresponding analyses
#' @export
#'
#' @examples
#' \dontrun{
#' knit_MT_grSeq_table(hyp_testing_dataset)
#' }
#' @keywords internal
knit_MT_grSeq_table <- function(hyp_testing_dataset, digits = 5, include_nominalPvalx2 = TRUE) {
  df <- hyp_testing_dataset %>%
    dplyr::select("hypNames", "alpha",
                  "Analysis", "timing_vec", "nominalPval", "nominalPvalx2", "deltaHat", "pow") %>%
    dplyr::rename(
      'Local alpha level'       = .data$alpha,
      'Info fraction'           = .data$timing_vec,
      'Nominal p-val (1-sided)' = .data$nominalPval,
      '2 x Nominal p-val'       = .data$nominalPvalx2,
      'Hurdle delta'            = .data$deltaHat,
      'Power'                   = .data$pow
    )
  if (include_nominalPvalx2==FALSE)
    df <- df %>% dplyr::select( - '2 x Nominal p-val')

  if (knitr::is_html_output()) {
    df[,-1] %>%
      knitr::kable("html", escape = F,
            align   = "c",
            digits  = digits,
            caption = "Efficacy p-value Boundaries") %>%
      kableExtra::kable_styling() %>%
      kableExtra::pack_rows(
        index         = table(forcats::fct_inorder(df$hypNames)),
        label_row_css = "text-align: left;"
      )  %>%
      kableExtra::column_spec(column = 1, underline = TRUE, width = '3cm')
    # collapse_rows(columns = 1, valign = "middle")
  } else if (knitr::is_latex_output()) {
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    df[,-1] %>%
      dplyr::mutate_all(kableExtra::linebreak) %>%
      knitr::kable(
        "latex",
        booktabs  = T,
        escape    = F,
        longtable = TRUE,
        caption   = "Efficacy p-value Boundaries"
      ) %>%
      kableExtra::kable_styling(latex_options = c("hold_position", "repeat_header")) %>%
      kableExtra::pack_rows(index = table(forcats::fct_inorder(df$hypNames)))
  } else if (knitr::pandoc_to("docx")) {
    #require(flextable)
    df <-
      data.frame(lapply(df, function(x) {
        gsub("<br>", "\n", x)
      }), stringsAsFactors = F)
    flextable::flextable(df)
  }
}

# spending functions for delayed alpha recycling ----

# Need to set a special spending function to account for the delayed recycling
# code Xi and Tamhane (2013) delayed alpha recycling for MT in grSeq
sfRecycle <- function (alpha, t, param)
{
  x <- list(
    name = "User-defined 2-stage",
    param = param,
    parname =  c(
      "initSigLevel",
      "tStar",
      "sfInit",
      "sfInitParam",
      # spend fun before tStar
      "sfStar",
      "sfStarParam"
    ),
    # spend fun after  tStar
    sf = sfRecycle,
    spend = NULL,
    bound = NULL,
    prob = NULL
  )
  class(x) <- "spendfn"

  gsDesign::checkScalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  gsDesign::checkVector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  if (length(param) != 6)
    stop("sfRecycle : parameter must be of length 6")
  gsDesign::checkScalar(param$initSigLevel, "numeric", c(0, Inf), c(FALSE, FALSE))
  gsDesign::checkVector(param$tStar, "numeric", c(0, Inf), c(TRUE, FALSE))

  t[t > 1] <- 1
  newSigLevel  <- alpha
  initSigLevel <- param$initSigLevel
  tStar        <- param$tStar
  sfInit       <- param$sfInit
  sfInitParam  <- param$sfInitParam
  sfStar       <- param$sfStar
  sfStarParam  <- param$sfStarParam


  if (sfStar(alpha, t = tStar, sfStarParam)$spend == sfStar(alpha, 1, sfStarParam)$spend) {
    stop("No increase of sfStar from tStar to 1\n")
    # x$spend <- ifelse(t<=tStar, sfInit(alpha,t=t,param=sfInitParam)$spend,alpha)
    # return(x)
  }

  # to implement (8), (9) equations in Xi (2015)
  getGammaStar <- function(y) {
    # Note that step type spending function might not be solved for gammaStar
    y - do.call(sfStar, list(
      alpha = y,
      t = tStar,
      param = sfStarParam
    ))$spend - newSigLevel + alpha1
  }

  if (alpha <= param$initSigLevel) {
    # no delayed recycling if no increase in available alpha
    x$spend <-
      do.call(sfInit, list(
        t     = t,
        alpha = alpha,
        param = sfInitParam
      ))$spend
    return(x)
  }

  alpha1 <-
    do.call(sfInit,
            list(alpha = initSigLevel, t = tStar, param = sfInitParam))$spend
  gammaStar <-
    stats::uniroot(getGammaStar, c(.Machine$double.eps, 1e16))$root
  # gammaStar <- (newSigLevel - alpha1) / (1-log(1+(exp(1)-1)*tStar))
  x$spend <- ifelse(
    t >= tStar,
    do.call(
      sfInit,
      list(t = tStar, alpha = initSigLevel, param = sfInitParam)
    )$spend +
      do.call(sfStar, list(
        t = t,     alpha = gammaStar, param = sfStarParam
      ))$spend -
      do.call(sfStar, list(
        t = tStar, alpha = gammaStar, param = sfStarParam
      ))$spend,
    do.call(sfInit, list(
      t = t, alpha = initSigLevel, param = sfInitParam
    ))$spend
  )
  x
}

sfSuperPower <- function (alpha, t, param)
{
  gsDesign::checkScalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  par <-  param[1] - param[2] * alpha
  gsDesign::checkScalar(par, "numeric", c(0, 15), c(FALSE, TRUE))
  gsDesign::checkVector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  x <-
    list(
      name    = "Augmented Power",
      param   = param,
      parname = c("A", "B"),
      sf      = gsDesign::sfPower,
      spend   = alpha * t ^ par,
      bound   = NULL,
      prob    = NULL
    )
  class(x) <- "spendfn"
  x
}

sfBonferroni <- function (alpha, t, param)
{
  x <- list(name  = "Bonferroni adjustment",
            param = param, parname = "Weights",
            sf    = sfBonferroni,
            spend = NULL, bound = NULL, prob = NULL)
  class(x) <- "spendfn"
  gsDesign::checkScalar(alpha, "numeric", c(0, Inf), c(FALSE, FALSE))
  gsDesign::checkVector(t, "numeric", c(0, Inf), c(TRUE, FALSE))
  t[t > 1] <- 1
  k <- length(t)
  j <- length(param)
  if (j == k - 1) {
    x$param <- c(param, 1)
    j <- k
  }
  if (j != k) {
    stop("Cumulative user-specified proportion of spending must be specified for each interim analysis")
  }
  if (!is.numeric(param)) {
    stop("Numeric user-specified spending levels not given")
  }
  incspend <- x$param - c(0, x$param[1:k - 1])
  if (min(incspend) < 0) {
    stop("Cumulative user-specified spending levels must be non-decreasing with each analysis")
  }
  if (max(x$param) > 1) {
    stop("Cumulative user-specified spending must be >= 0 and <= 1")
  }
  x$spend <- alpha * x$param
  x
}

# Define the function to extract a name of spending function picking on code
get_sf_name <- function(func) {
  if (is.character(func)) {
    func <- get(func)  # Get the function object by its name
  }
  if (!is.function(func)) {
    stop("Input must be a function name or a function.")
  }
  func_body <- body(func) # Extract the function body
  # Convert the function body to a single character string
  func_code_string <- paste(deparse(func_body), collapse = " ")

  # Look for the pattern "name = " and extract the value
  pattern <- "name\\s*=\\s*\"([^\"]+)\""
  match <- regmatches(func_code_string, regexec(pattern, func_code_string))
  res <- match[[1]][2]
  if (is.na(res)){
    res <- do.call(func, list(0.1, 1, NULL))$name # dummy call
  }
  return(res)
}

# Declare the operator explicitly to avoid notes
# (add this in your R script or package setup)
`%m+%` <- lubridate::`%m+%`

# generate gtable of timelines of IA by hypotheses
timeline_gtable <- function(D, startDate = "2022-10-12", lpi = NULL) {
  D <-
    dplyr::filter(D, .data$regiment == unique(D$regiment)[1]) # limit plotting to just a single regiment
  # Set hypothesis names in single vector
  hypothesis_names       <- D$hypNames
  J                      <- length(hypothesis_names)
  times <- sort(unique(unlist(D$iaTime)))
  # Set number of analyses for each hypothesis
  K_j                    <- sapply(D$infoFr, length)
  # Set maximum number of analyses
  K                      <- max(K_j)
  # Ensure times has at least K entries to prevent index out of bounds
  if (length(times) < K) {
    # Pad with the last time value to ensure we have K entries
    times <- c(times, rep(times[length(times)], K - length(times)))
  }
  # Set colours for each stage to use in final plot
  if (K == 2) {
    colours              <- c("#9ECAE1FF", "#3182BDFF")
  } else {
    colours              <-
      grDevices::adjustcolor(RColorBrewer::brewer.pal(n = K, name = "Blues"))
  }
  times_with_zero        <- c(0, times)
  # Set vector containing stage names
  event                  <- paste0("Stage ", 1:K)
  # Start at year 0 as will only plot time information in months from FPI (could
  # modify to use actual dates if desired)
  if (is.null(startDate))
    startDate <- "0000-01-01"

  start <- end <- as.POSIXct(rep(startDate, K))

  # Set start and end times of each stage
  for (k in 1:K) {
    start[k]             <- start[k] %m+% months(times_with_zero[k])
    end[k]               <-
      end[k] %m+% months(times_with_zero[k + 1])
  }
  # Group for the above information is the Stage and colours are those specified
  # earlier
  group                  <- rep("Stages", K)
  color                  <- colours
  # Add information on LPI
  if (!is.null(lpi)) {
    event                  <- c(event, paste0(lpi, " mo"))
    lpi_time               <-
      as.POSIXct(startDate) %m+% months(lpi)
    start                  <- c(start, lpi_time)
    end                    <- c(end, lpi_time)
    group                  <- c(group, "LPI")
    color                  <- c(color, "black")
  }
  # Loop over the hypotheses and add the information for each of them
  for (j in 1:J) {
    # Extract sample size / event information up to data maturity
    ss_j                 <-
      D$hypN[[j]] * D$infoFr[[j]]              # hypotheses[[j]]$ss[1:K_j[j]]
    # Convert to IF (rounding so it plots better)
    if_j                 <-
      round(D$infoFr[[j]] * 100) # round(100*ss_j/ss_j[K_j[j]])
    # The group is just the hypothesis name
    group                <-
      c(group, rep(D$hypNames[j], K_j[j])) # c(group, rep(hypothesis_names[j], K_j[j]))
    color                <- c(color, colours[1:K_j[j]])
    # Event information is based on combining various strings together
    if (K_j[j] > 1) {
      event_j          <- c(paste0("IA", 1:(K_j[j] - 1), ": "), "FA: ")
    } else {
      event_j          <- "FA: "
    }
    if (D$endpointType[j] == "Binomial") {
      text_j           <- "pts"
    } else if (D$endpointType[j] == "TTE") {
      text_j           <- "ev"
    }
    if (K_j[j] > 1) {
      event_j          <-
        paste0(
          event_j,
          times_with_zero[2:(K_j[j] + 1)],
          " mo,\n",
          ss_j,
          " ",
          text_j,
          c(rep(", ", K_j[j] - 1), ""),
          c(if_j[-K_j[j]], ""),
          c(rep("%IF", K_j[j] - 1), "")
        )
    } else {
      event_j          <-
        paste0(
          event_j,
          times_with_zero[2:(K_j[j] + 1)],
          " mo,\n",
          ss_j,
          " ",
          text_j,
          c(rep(", ", K_j[j] - 1), ""),
          c(if_j[-K_j[j]], ""),
          c(rep("%IF", K_j[j] - 1), "")
        )
    }
    event                <- c(event, event_j)
    start                <- c(start, end[1:K_j[j]])
    end                  <- c(end, end[1:K_j[j]])
  }
  # Build data frame for plotting
  timeline_data  <- data.frame(
    event = event,
    start = start,
    group = group,
    end   = end,
    color = color
  )

  # Set x limit for the plot based on FA time
  x_limit <- lubridate::interval(startDate,
                                 max(as.POSIXct(timeline_data$end)) %m+% months(1)) %/%  months(12)
  # Build initial timeline plot
  p_timeline             <-
    vistime::gg_vistime(timeline_data, linewidth = 5) +
    ggplot2::scale_x_datetime(
      breaks =
        seq(as.POSIXct(min(startDate)),
            max(as.POSIXct(
              timeline_data$end
            )) %m+%
              months(1), "years"),
      labels = 12 * seq(0, x_limit, 1)
    ) +
    ggplot2::xlab("Months") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  # Not always needed, but ggplot_build() can help avoid overlapping labels
  p_build                <- ggplot2::ggplot_build(p_timeline)
  p_build$data[[4]]$size <- 3
  p_build$data[[5]]$size <- 3
  p_timeline             <- ggplot2::ggplot_gtable(p_build)
  return(p_timeline)
}

# main functions ----

#' Check input specifications for consistency (Work in Progress!)
#'
#' @param inputD An input dataset
#' @param G A graph object
#' @param enrollment An enrollment object
#'
#' @returns A list of objects for output tables and plots
#' @export
#'
#' @examples
#' \dontrun{
#' checkInput(inputD,G)
#' }
#' @keywords internal
checkInput <- function(inputD, G, enrollment = NULL) {
  N_rand <- cumsum(c(0, enrollment$rate * c(enrollment$duration)))
  numHyp <- nrow(inputD)
  Hs  <-
    inputD$iaSpec |> unlist(recursive = FALSE) |>
      purrr::map_dbl(purrr::pluck("H"))
  IFs <-
    inputD$iaSpec |> unlist(recursive = FALSE) |>
    purrr::map_dbl(purrr::pluck("atIF"))
  if (!all(Hs %in% 1:numHyp))
    stop("H index outside of range in inputD$iaSpec")
  if (!all(IFs > 0 &
           IFs <= 1))
    stop("Info fraction must be in (0,1]")
  # TODO
  # 1) check that enrollment and sample sizes for binary EPs are in agreement
  # 2) ensure that spending function is specified if there are IAs

  if (is.null(inputD$tag))
    stop("Missing tag field in inputD")
}

#' Input Processing
#'
#' @param inputD An input dataset with specifications
#'
#' @returns A list with main objects to be outputted: D, ia_details and hyp_testing_dataset
#' @export
#'
#' @examples
#' \dontrun{
#'   # Need to define inputD first
#'   exec_calc(inputD)
#' }
#' @keywords internal
exec_calc <- function(inputD) {
  checkInput(inputD, G)    # check the inputs TODO

  D <- inputD %>% dplyr::mutate(
    # effect on a natural scale, e.g, mean group difference or HR
    delta       = purrr::map_dbl(.data$endpointParam, getDelta),
    deltaStr    = purrr::map_chr(.data$endpointParam, getEffectSizeDetails),
    # get standardization factor to calculate effect sizes 'theta'
    standFactor = purrr::map2_dbl(.data$endpointParam, .data$allocRatio, getStandardizingCoef),
    logDelta  = purrr::map_lgl(.data$endpointParam, function(x)
      class(x) %in% c("tte_exp", "tte_pwe")),
    theta =  dplyr::if_else(.data$logDelta,-log(.data$delta) * .data$standFactor, .data$delta * .data$standFactor),
    hypNames = paste(.data$id, .data$tag, sep = ": ")  # combine ID and tag names into full name
  )

  # derive calendar time of IAs usisng 'iaSpec', enrollment, and endpoint maturity time in 'endpointParam'
  D$iaTime <- deriveCalendarTime(D)
  # derive sample sizes (number of events)
  D$n.I   <- derive_Ns(D)
  # derive information fractions
  D$infoFr  <-
    lapply(
      D$n.I,
      FUN = function(x)
        unique(round(x / x[length(x)], idigits))
    )
  D$hypN    <- sapply(
    D$n.I,
    FUN = function(x)
      x[length(x)]
  )

  # add description columns
  descr_fields <-
    match(c("regiment", "ep", "suffix"), names(D), nomatch = 0)
  if (!any(descr_fields))
    descr_fields <- match(c("tag"), names(D), nomatch = 0)
  D$descr <- dplyr::select(D, descr_fields) %>% purrr::pmap_chr(.f = paste)
  D <-
    tibble::add_column(D, grSeqTestingCh = sapply(D$grSeqTesting, print_sfInfo))

  ia_details <- D %>%
    dplyr::distinct(.data$id, .data$iaSpec, .keep_all = TRUE)     %>%
    dplyr::mutate(id_tag = paste0(.data$id, " (", .data$tag, ")"))   %>%
    dplyr::select("id_tag", "iaSpec", "iaTime", "n.I", "infoFr")        %>%
    dplyr::mutate(iaSpec = purrr::map(.data$iaSpec, function(x) {
      lapply(x, tibble::as_tibble)
    })) %>%
    tidyr::unnest(c(.data$iaSpec, .data$iaTime, .data$n.I, .data$infoFr)) %>%
    tidyr::unnest(.data$iaSpec) %>%
    dplyr::group_by(.data$id_tag) %>%
    dplyr::mutate(ia = dplyr::row_number(),
                  criterion = paste0("H", .data$H, " at information fraction ",
                                     round(.data$atIF, idigits))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$iaTime) %>%
    dplyr::mutate(ia_ind = dplyr::cur_group_id()) %>%
    dplyr::ungroup()

  # the main call to calculate details of MT in group sequential testing
  hyp_testing_dataset <-
    report_MT_grSeq(G, D, pdigits = pdigits, idigits = idigits)
  return(list(
    D                   = D ,
    ia_details          = ia_details,
    hyp_testing_dataset = hyp_testing_dataset
  ))
}


