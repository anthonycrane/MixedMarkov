# File: Drugs_Functions.R
# Purpose: Defines convenience functions for use in the mHMM project.

AnySuccess.Z <- function(data, subject){
  # Function: AnySuccess.Z(data, subject)
  # Purpose:  Identifies if there is a sequence of 14 consecutive 0's for a
  #           given subject in the data frame., for 'Z' variable.
  # Input:    data    - A data frame
  #           subject - Subject id in data frame
  # Output:   '1' if the subjects had a 14 day non-using streak
  #           '0' otherwise
  sequences  <- rle(data$Z[data$Subject == subject])
  sequences.nouse.length <- sequences$lengths[which(sequences$values == 0)]
  n.success <- length(sequences.nouse.length[sequences.nouse.length >= 14])
  return(ifelse(n.success > 0,1, 0))
}

AnySuccess.C <- function(data, subject){
  # Function: AnySuccess.C(data, subject)
  # Purpose:  Identifies if there is a sequence of 14 consecutive 0's for a
  #           given subject in the data frame., for 'C' variable.
  # Input:    data    - A data frame
  #           subject - Subject id in data frame
  # Output:   '1' if the subjects had a 14 day non-using streak
  #           '0' otherwise
  sequences  <- rle(data$C[data$Subject == subject])
  sequences.nouse.length <- sequences$lengths[which(sequences$values == 0)]
  n.success <- length(sequences.nouse.length[sequences.nouse.length >= 14])
  return(ifelse(n.success > 0,1, 0))
}

AnySuccess.C.fit <- function(data, subject){
  # Function: AnySuccess.C.fit(data, subject)
  # Purpose:  Identifies if there is a sequence of 14 consecutive 0's for a
  #           given subject in the data frame., for 'C.fit' variable.
  # Input:    data    - A data frame
  #           subject - Subject id in data frame
  # Output:   '1' if the subjects had a 14 day non-using streak
  #           '0' otherwise
  sequences  <- rle(data$C.fit[data$Subject == subject])
  sequences.nouse.length <- sequences$lengths[which(sequences$values == 0)]
  n.success <- length(sequences.nouse.length[sequences.nouse.length >= 14])
  return(ifelse(n.success > 0,1, 0))
}

invlogit <- function (x) {
  # Function: invlogit(x)
  # Purpose:  Compute inverse logistic value of an adds.
  # Input:    An odds
  # Output:   A probability
  return(exp(x)/(1+exp(x)))
}

visits.gen <- function(X.it, C.itm1, Beta, beta.0, beta.1, b.0, lambda.0, lambda.1) {
  # Function: visits.gen(X.it, C.itm1, Beta, beta.0,beta.1, lambda.0, lambda.1)
  # Purpose:  Generate the response at time t for subject i, given covariates
  #           at time t and hidden state at time t - 1.
  # Input:    X.it     - Treatment group time t
  #           C.itm1   - Hidden state at time t - 1
  #           Beta     - Treatment effect of being in state 1 compared to state
  #                      0 (log OR)
  #           beta.0   - Intercept for being in state 1 compared to state 0
  #                      (log odds of being in state 0, given in state 0 at
  #                      previous time point and in control group)
  #           beta.1   - Log OR of being in state 1 given in state 1 at previous
  #                      time point
  #               b.0  - random intercept for subject
  #           lambda.0 - Log odds of observing drug use while in state 0
  #           lambda.1 - Log OR of observing drug use in state 1 compared to
  #                      state 0
  # Output:   A list of length 2 with first element being the state at time t
  #           and the second element being observed use at time t
  p.1 <- invlogit(beta.0 + beta.1*C.itm1 + Beta*X.it + b.0)
  C.it <- rbinom(n = 1, size = 1, p = p.1)
  theta.it <- invlogit(lambda.0 + lambda.1*C.it)
  Z.it <- rbinom(n = 1, size = 1, p = theta.it)
  return(list(C.it,Z.it))
}


Data.gen <- function(X.i, nday, Beta, beta.0, beta.1, b.0, lambda.0, lambda.1) {
  # Function: Data.gen(X.i, nday, Beta, beta.0, beta.1, lambda.0, lambda.1)
  # Purpose:  Given covatiates across all measurements for a subject, generate
  #           the response vector and hidden states for each measurement time.
  # Input:    X.i      - nday x 3 matrix, first column is Subject, Second Column
  #                      is Day, 3rd column is Treatment group
  #           nday     - Number of days of observations
  #           Beta     - Treatment effect of being in state 1 compared to state
  #                      0 (log OR)
  #           beta.0   - Intercept for being in state 1 compared to state 0
  #                      (log odds of being in state 0, given in state 0 at
  #                      previous time point and in control group)
  #           beta.1   - Log OR of being in state 1 given in state 1 at previous
  #                      time point
  #               b.0  - random intercept for subject
  #           lambda.0 - Log odds of observing drug use while in state 0
  #           lambda.1 - Log OR of observing drug use in state 1 compared to
  #                      state 0
  # Output:   nday x 5 matrix; columns are: Subj, Day, Treatment, C (state),
  #           Z (response)
  #
  # D.i is a nday x 5 matrix; columns are: Subj, Day, Treatment, C (state),
  # Z (response)
  D.i <- cbind(X.i, C = rep(NA, nday), Z = rep(NA, nday))
  for (row in 1:nrow(D.i)) {
    if (row == 1) {
      X.it <- D.i[row, 3]
      C.itm1 <- 1 # We assume all subjects were in state 1 prior to entry
      V.gen <- visits.gen(X.it = X.it, C.itm1 = C.itm1, Beta = Beta,
                          beta.0 = beta.0, beta.1 = beta.1, b.0 = b.0,
                          lambda.0 = lambda.0, lambda.1 = lambda.1)
      D.i[row, 4] <- V.gen[[1]]
      D.i[row, 5] <- V.gen[[2]]
    } else {
      X.it <- D.i[row, 3]
      C.itm1 <- D.i[row - 1, 4]
      V.gen <- visits.gen(X.it = X.it, C.itm1 = C.itm1, Beta = Beta,
                          beta.0 = beta.0, beta.1 = beta.1, b.0 = b.0,
                          lambda.0 = lambda.0, lambda.1 = lambda.1)
      D.i[row, 4] <- V.gen[[1]]
      D.i[row, 5] <- V.gen[[2]]
    }
  }
  return(D.i)
}

expected <- function(x){
  (rowSums(x)%*%t(colSums(x)))/sum(x)
}