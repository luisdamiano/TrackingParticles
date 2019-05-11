#' @useDynLib TrackingParticles
NULL

#' Measurements from two passive sensors tracking a moving vehicle.
#'
#' @format A data frame with 11027 rows and 2 numeric variables:
#' \describe{
#'   \item{a1}{Angle in radians of a moving vehicle as measured from one
#'   passive sensor located at -93.2494663765932, 41.5563518606521.}
#'   \item{a2}{Angle in radians of a moving vehicle as measured from another
#'   passive sensor located at -93.2475338232000, 41.5576632356000.}
#' }
#' @source \url{https://github.com/ISU-STRIPS/STRIPSyield/}
"vehicle"

#' Compute posterior mean of the latent state (position and velocity).
#'
#' For a given parameter vector, this function runs a Particle Filter to
#' estimate the posterior mean of the latent state for the bearing-only
#' tracking problem with two passive sensors.
#'
#' @param y A two-column matrix with the measurements.
#' @param dt The time step between observations.
#' @param location1 A two-element vector with the longitude (x) and latitude
#' (y) of the first sensor.
#' @param location2 A two-element vector with the longitude (x) and latitude
#' (y) of the second sensor.
#' @param sr The variance of the measurement model error.
#' @param q1 The first difussion constant of the state model.
#' @param q2 The second difussion constant of the state model.
#' @param statepriorMu A two-element vector with the longitude (x) and latitude
#' (y) of the location where the state prior density should be centered.
#' @param statepriorCholesky A four-element vector with the diagonal of the
#' Cholesky factor corresponding to the variance of the state prior distribution.
#' @param importanceCholesky A four-element vector with the diagonal of the
#' Cholesky factor corresponding to the variance of the importance distribution.
#' @param nParticles An integer with the number of particles.
#'
#' @return A named list with four elements.
#' `noiseless` is a T x 2 matrix with the noiseless approximation of the
#' vehicle position (assumes no noise and velocity equal to zero).
#' `stateMean` is a T x 4 matrix with the posterior mean of the latent state
#' at each time step.
#' `weights` is a T x nParticles matrix with the normalized weights.
#' `ess` is a T-sized vector with the effective sample size at each time step.
#' @note The resampling step is currently not implemented. Expect particle
#' degeneracy (i.e. rapidly decaying ESS).
#' @seealso \code{\link{plot.filtered}{plot}}
#' @export
particle_filter <- function(y, dt, location1, location2, sr, q1, q2,
                            statepriorMu, statepriorCholesky,
                            importanceCholesky, nParticles) {
  # Ready...
  DIM_MEASUREMENT <- 2
  DIM_STATE       <- 4
  y               <- as.matrix(y)
  RT              <- nrow(y)
  location1       <- as.numeric(location1)
  location2       <- as.numeric(location2)

  # Steady...
  if (ncol(y) != DIM_MEASUREMENT)
    stop(sprintf("`y` must be a matrix with %i columns.", DIM_MEASUREMENT))

  if ((length(location1) != DIM_MEASUREMENT) ||
      (length(location2) != DIM_MEASUREMENT))
    stop(sprintf("`location1` and `location2` must be two vectors of size %i.",
                 DIM_MEASUREMENT))

  if (min(sr, q1, q2, statepriorCholesky, importanceCholesky) < 0)
    stop("Variance components may only take positive values.")

  # Go!
  out <- .C(
    "Rfilter",
    Ry1                   = as.double(y[, 1]),
    Ry2                   = as.double(y[, 2]),
    RT                    = as.integer(RT),
    LOCATION_1_X          = as.double(location1[1]),
    LOCATION_1_Y          = as.double(location1[2]),
    LOCATION_2_X          = as.double(location2[1]),
    LOCATION_2_Y          = as.double(location2[2]),
    DT                    = as.double(dt),
    MEASUREMENT_ERROR_1   = as.double(sr),
    STATE_DIFFUSION_1     = as.double(q1),
    STATE_DIFFUSION_2     = as.double(q2),
    STATEPRIOR_MU_X       = as.double(statepriorMu[1]),
    STATEPRIOR_MU_Y       = as.double(statepriorMu[2]),
    STATEPRIOR_L_00       = as.double(statepriorCholesky[1]),
    STATEPRIOR_L_11       = as.double(statepriorCholesky[2]),
    STATEPRIOR_L_22       = as.double(statepriorCholesky[3]),
    STATEPRIOR_L_33       = as.double(statepriorCholesky[4]),
    IMPORTANCE_L_00       = as.double(importanceCholesky[1]),
    IMPORTANCE_L_11       = as.double(importanceCholesky[2]),
    IMPORTANCE_L_22       = as.double(importanceCholesky[3]),
    IMPORTANCE_L_33       = as.double(importanceCholesky[4]),
    NPARTICLES            = as.integer(nParticles),
    RnoiselessOut         = as.double(
      matrix(0, nrow = RT, ncol = DIM_MEASUREMENT)),
    RxMeanOut             = as.double(
      matrix(0, nrow = RT + 1, ncol = DIM_STATE)),
    RwOut                 = as.double(
      matrix(0, nrow = RT + 1, ncol = nParticles)),
    RessOut               = as.double(
      vector("numeric", RT + 1)),
    PACKAGE = "TrackingParticles"
  )

  # Return
  structure(
    list(
      noiseless = matrix(out$RnoiselessOut, RT, DIM_MEASUREMENT),
      stateMean = matrix(out$RxMeanOut, RT + 1, DIM_STATE)[-1, ],
      weights   = matrix(out$RwOut, RT + 1, nParticles)[-1, ],
      ess       = out$RessOut[-1]
    ),
    class = c("filtered")
  )
}

#' Plot the posterior mean of the position as estimated by the Particle Filter.
#'
#' @param x An object returned by the function \code{\link{particle_filter}}.
#' @param ... Currently nothing.
#' @note The function may hide values of the state posterior mean that are
#' estimated as zero (funny estimates).
#' @export
#' @importFrom graphics par plot
plot.filtered <- function(x, ...) {
  ind <- !((x$stateMean[, 1] == 0) | (x$stateMean[, 2] == 0))

  par(mfrow = c(1, 2))

  plot(
    x$noiseless,
    main = "Noiseless approximation",
    xlab = "Latitude",
    ylab = "Longitude",
    ...
  )

  plot(
    x$stateMean[ind, 1:2],
    main = "Particle filter",
    xlab = "Latitude",
    ylab = "Longitude",
    ...
  )
}
