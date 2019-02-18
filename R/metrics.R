## This file contains the functions that calculate the QC metrics for peaks and peak groups.


#' Compute maximum intensity of each transition peak in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and returns a numeric vector of maximum intensity for individual transition peaks in the peak group.
#'
#' @param peak A peak group of class peakObj
#'
#' @return  A numeric vector of maximum intensity for individual transition peaks in the peak group
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' transition.max.intensities <- CalculatePeakMaxIntensity(peak)

CalculatePeakMaxIntensity <- function(peak, ...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format = simpleError("CalculatePeakMaxIntensity: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)

  # function body  ---------------------------------------

  # max peak intensity is calculated
  r.max.intensity <- sapply(peak@sig,max)

  # return output
  return(r.max.intensity)
}

#' Compute maximum intensity at the peak boundary of each transition peak in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and returns a numeric vector of maximum intensity at peak boundary for individual transition peaks in the peak group.
#'
#' @param peak A peak group of class peakObj
#'
#' @return A numeric vector of maximum intensity at peak boundary for individual transition peaks in the peak group
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' transition.max.at.boundary.intensities <- CalculateMaxBoundaryIntensity(peak)

CalculateMaxBoundaryIntensity <- function(peak, ...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculateMaxBoundaryIntensity: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)

  # function body  ---------------------------------------

  #   max intensity at peak boundary is calculated
  r.max.boundary.intensity <- sapply(peak@sig,function(x) max(head(x,1),tail(x,1)))

  # return output
  return(r.max.boundary.intensity)
}


#' Compute elution shift of a transition peak relative to a reference peak.
#'
#' The function takes two numeric vectors representing transition peaks as input and calculates their relative shift. Additionally, it needs to timepoints corresponding to the input transition peaks. The shift between these transition peaks is defined as the difference between time at max for the two peaks normalized by the peak width at base. This function should only be used for peaks that belong to the same peak group and have the same peak boundary and timepoints. Although this function can be used independently, it is meant to be called by CalculatePeakElutionShift. For high quality peaks, the elution shift is expected to be close to 0.
#'
#' @param sig1 A numeric vector representing a transition peak whose relative elution shift to the reference peak is to be calculated.
#' @param sig2 A numeric vector representing the reference peak.
#' @param time A numeric vector representing the timepoints series of the peaks
#'
#' @return Numeric value of the relative shift (no unit) between sig1 and sig2.
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' transition.shift.y6 <- CalculateElutionShift(sig1 = peak@@sig$y6.1.light, sig2 = peak@@sig$y6.1.heavy,time = peak@@time)


CalculateElutionShift <- function(sig1,sig2,time,...) {

  # error and warning handling ---------------------------------------

  #   input errors: if the inputs are not numeric or  not equal in length
  error.input.format <- simpleError("CalculateElutionShift: sig1, sig2 and time should be numeric vectors of equal length")

  if (!is.numeric(sig1)) stop(error.input.format) # input should be a numeric vector
  if (!is.numeric(sig2)) stop(error.input.format) # input should be a numeric vector
  if (!is.numeric(time)) stop(error.input.format) # input should be a numeric vector
  if (length(sig1) != length(sig2)) stop(error.input.format) # the two numeric vectors should be of equal length
  if (length(sig1) != length(time)) stop(error.input.format) # numeric vectors of intensity should have equal length with the time vector



  # function body  ---------------------------------------

  # if at least one of the inputs is all zeros, return 0 for the peak shift
  if ((abs(sum(sig1)) == 0) || (abs(sum(sig2)) == 0)) {
    return(round(0,digits = 4))
  }

  #   peak shift is the difference between time at max of the two signals. If there are multiple points with max intensity the smallest shift is used.
  peak.shift <- min(abs(time[which(sig1 == max(sig1))] - time[which(sig2 == max(sig2))])/(tail(time,1) - head(time,1)))

  # return output
  return(round(peak.shift,digits = 4))
}

#' Compute pairwise elution shift of transition peaks in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and calculates the pairwise shift between transition peaks in the peak group using the CalculateElutionShift function. The shift between two transition peaks is defined as the difference between time at max for the two peaks normalized by the peak width at base. The shift between an individual transition peak relative to the peak group is defined as the different between time at mas of the transition peak and the median of time at max of all the peaks in the peak group normalized by the peak width at base. For high quality peaks, the elution shift is expected to be close to 0.
#'
#' @param peak A peak group of class peakObj
#'
#' @return A list with the following objects:
#'               r.shift: A numeric matrix where the non-diagonal elements represent the pairwise shift between transition peaks and diagonal  elements represent the shift between each transition peak relative to the median time at max of all the transition peaks in the peak group.
#'
#' @export
#'
#' @importFrom tidyr spread
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' peak.elution.shift <- CalculatePeakElutionShift(peak)

CalculatePeakElutionShift <- function(peak,...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculatePeakElutionShift: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)


  # function body  ---------------------------------------

  #   elution shift is calculated for each pair of transitions in the peak group using the CalculateElutionShift function

  r.shift <- data.frame(shift = mapply(CalculateElutionShift,t(combn(peak@sig,2))[,1],t(combn(peak@sig,2))[,2],data.frame(peak@time)),
                       ion = combn(colnames(peak@sig),2)[1,],
                       ion2 = combn(colnames(peak@sig),2)[2,]) %>%
    spread(key = ion, value = shift)
  rownames(r.shift) <- r.shift$ion2
  r.shift <- r.shift %>% select(-ion2)
  r.shift <- cbind(rbind(r.shift,rep(NA,ncol(r.shift))),rep(NA,ncol(r.shift) + 1))
  colnames(r.shift)[ncol(r.shift)] <- setdiff(colnames(peak@sig),colnames(r.shift))
  rownames(r.shift)[nrow(r.shift)] <- setdiff(colnames(peak@sig),rownames(r.shift))

  # the rows and columns are ordered according to the order of the transitions in the peak object
  r.shift <- r.shift[names(peak@area),names(peak@area)]

  # shift for each transition vs the peak groups is determined by the difference between time at max for each transition and the median of time at max for all the transitions in the peak group.
  max.intensity.times <- peak@time[unlist(sapply(peak@sig,function(sig) min(which(sig == max(sig)))))]
  diag(r.shift) <- round(abs(max.intensity.times - median(max.intensity.times))/(tail(peak@time,1) - head(peak@time,1)),digits = 4)

  # the NA values in this matrix correspond to transition pairs that are ordered differently. NAs fpr (tr1,tr2) transitions pairs are replaced by the peak elution shift calculated for (tr2,tr1) pairs
  ind_na <- which(is.na(r.shift), TRUE)
  if (nrow(ind_na) == 1)
    r.shift[ind_na] <- r.shift[ind_na[1,2],ind_na[1,1]]
  else
    r.shift[ind_na] <- r.shift[ind_na[,2:1]]

  # format output
  r.shift <- as.matrix(r.shift)
  shift <- list(r.shift = r.shift)

  # return output
  return(shift)
}


#' Compute the jaggedness score for a transition peak.
#'
#' The function takes a numeric vector representing a transition peak as input and calculates the jaggedness score for the peak. The jaggedness score is defined as the fraction of timepoints where the signal changes direction, excluding the change of direction at the peak apex. This function can be called independently, however it is mainly meant to be called by CalculatePeakJaggedness to calculate the jaggedness metric for a peak group. For high quality peaks, the jaggedness score is expected to be close to 0.
#'
#' @param sig A numeric vector representing a transition peak
#' @param flatness.factor A numeric parameter between 0 and 1 that determines the sensitivity of the jaggedness score to low levels of noise. To avoid high jaggedness scores due to small levels of noise, near-flat ranges in the peak are artificially flattened before calculating the jaggedness score. A range is defined as near-flat and thus flattened when the difference between intensities of adjacent time points is smaller than flatness.factor times the peak maximum intensity.
#'
#' @return The numeric value of jaggedness score for the transition peak, which is defined as the fraction of timepoints where the signal changes direction
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[196]]
#' transition.jaggedness.y6 <- CalculateJaggedness(peak@@sig$y6.1.heavy, flatness.factor = 0.05)
#' PlotChromPeak(peak,transition.list = c("y6"))

CalculateJaggedness <- function(sig, flatness.factor = 0.05, ...) {

  # error and warning handling ---------------------------------------

  #   input errors: if the input vectors are na or too short
  error.input.format <- simpleError("CalculateJaggedness: input sig should be a numeric vector with at least 3 elements")

  if (length(sig) < 3) stop(error.input.format)

  error.flatness.format <- simpleError("CalculateJaggedness: flatness.factor should be a numeric value between 0 and 1")

  if (flatness.factor < 0 || flatness.factor > 1) stop(error.flatness.format)


  # function body  ---------------------------------------

  #   changes in the sign of differential of the peak is used to quantify the jaggedness of the peak
  diff.sig <- diff(sig)

  #   the near-flat ranges of peak, where the differential is less than flatness.factor x peak max are assumed to be flat in measurements of jaggedness
  diff.sig[which(abs(diff.sig) < flatness.factor*max(abs(sig)))] = 0

  #   jaggeddness is calculated
  jaggedness <- (sum(abs(diff(sign(diff.sig))) > 1) - 1)/(length(diff.sig) - 1)

  #   if jaggedness is negative return zero
  jaggedness <- max(0,jaggedness)

  # return output
  return(round(jaggedness,digits = 4))
}

#' Compute jaggedness scores for transition peaks in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and calculates the jaggedness score for the transition peaks in the peak group using the CalculateJaggedness function. The jaggedness score is defined as the fraction of timepoints where the signal changes direction, excluding the change of direction at the peak apex. For high quality peaks, the jaggedness score is expected to be close to 0.
#'
#' @param peak A peak group of class peakObj
#' @param flatness.factor A numeric parameter between 0 and 1 that determines the sensitivity of the jaggedness score to low levels of noise. To avoid high jaggedness scores due to small levels of noise, near-flat ranges in the peak are artificially flattened before calculating the jaggedness score. A range is defined as near-flat and thus flattened when the difference between intensities of adjacent time points is smaller than flatness.factor times the peak maximum intensity.
#'
#' @return A list with the following objects:
#'               r.jaggedness: A numeric vector of jaggedness scores for each transition in the peak group.
#'               peak.jaggedness: Mean of the r.jaggedness vector. This score represents the overall jaggedness of all the transitions in the peak group.
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[196]]
#' peak.group.jaggedness <- CalculatePeakJaggedness(peak)

CalculatePeakJaggedness <- function(peak, flatness.factor = 0.05, ...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculatePeakJaggedness: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)

  # function body  ---------------------------------------

  #   jaggeddness is calculated using the CalculateJaggedness function
  r.jaggedness <- mapply(CalculateJaggedness,peak@sig,flatness.factor = flatness.factor)
  peak.jaggedness <- round(mean(r.jaggedness),digits = 4)

  # format output
  jaggedness <- list(r.jaggedness = r.jaggedness, peak.jaggedness = peak.jaggedness)

  # return output
  return(jaggedness)
}

#' Compute similarity scores for transition peaks in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and calculates pairwise similarity scores for the transition peaks in the peak group. Similarity score is defined as the Pearson correlation coefficient between the signal intensity vectors of the two peaks. For high quality peaks, the similarity score is expected to be close to 1.
#'
#' @param peak A peak group of class peakObj
#'
#' @return A list with the following objects:
#'               r.similarity: A numeric matrix of pairwise similarity scores for transition peaks in the peak group.
#'               peak.jaggedness: Mean of the r.similarity. vector. This score represents the overall similarity of all the transitions in the peak group.
#' @export
#'
#' @importFrom psych corr.test
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' peak.group.similarity <- CalculatePeakShapeSimilarity(peak)

CalculatePeakShapeSimilarity <- function(peak,...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculatePeakShapeSimilarity: input peak group should be non-empty")

  if (is.na(peak)) stop(error.input.format)


  # function body  ---------------------------------------

  #   pearson correlation coefficient between sig and the ref peak
  r.similarity <- psych::corr.test(peak@sig,method = "pearson",adjust = "holm", alpha = .05,ci = TRUE)$r

  # NA values are imputed to 0. they are a result of all zero signals.
  r.similarity[is.na(r.similarity)] <- 0

  # with the exception of the diagonal NA values which must be imputed to 1.
  diag(r.similarity) <- 1
  peak.similarity <- round(mean(r.similarity),digits = 4)

  # return output
  similarity <- list(r.similarity = r.similarity, peak.similarity = peak.similarity)
  return(similarity)
}

#' Compute symmetry scores for transition peaks in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and calculates the symmetry score for the transition peaks in the peak group. The symmetry score is defined as the Pearson correlation coefficient between the left and right halves of the peak intensity vector divided at the timepoint in the middle. For high quality peaks, the symmetry score is expected to be close to 1.
#'
#' @param peak A peak group of class peakObj
#'
#' @return A list with the following objects:
#'               r.symmetry: A numeric vector of symmetry scores for each transition in the peak group.
#'               peak.symmetry: Mean of the r.symmetry vector. This score represents the overall symmetry of all the transitions in the peak group.
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' peak.group.symmetry <- CalculatePeakSymmetry(peak)

CalculatePeakSymmetry <- function(peak,...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculatePeakSymmetry: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)


  # function body  ---------------------------------------

  #   left and right half of the peak
  left <- peak@sig[1:floor(nrow(peak@sig)/2),]

  right <- peak@sig[seq(nrow(peak@sig),nrow(peak@sig) + 1 - floor(nrow(peak@sig)/2),by = -1),]

  #   pearson correlation coefficient between sig and the ref peak
  r.symmetry <- diag(suppressWarnings(cor(left,right,method = "pearson")))

  # replace NA values with 1. cor returns NA if there is missing data or if right and left half are identical as seen in all zero signals
  r.symmetry[is.na(r.symmetry)] <- 1
  peak.symmetry <- round(mean(r.symmetry),digits = 4)

  symmetry <- list(r.symmetry = r.symmetry, peak.symmetry = peak.symmetry)

  # return output
  return(symmetry)
}


#' Compute FWHM and ratio of FWHM to peak width at base for transition peaks in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and calculates the full-width at half-max (FWHM) and the FWHM to peak width at the base ratios for the transition peaks in the peak group. See \url{https://en.wikipedia.org/wiki/Full_width_at_half_maximum}
#'
#' @param peak A peak group of class peakObj
#'
#' @return A list with the following objects:
#'               r.fwhm: A numeric vector of FWHM for each transition in the peak group.
#'               peak.fwhm: Mean of the r.fwhm vector.
#'               r.fwhm2base: A numeric vector of FWHM to peak width at base ratio for each transition in the peak group.
#'               peak.fwhm2base: Mean of the r.fwhm2base vector.
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' peak.group.fwhm <- CalculateFWHM(peak)

CalculateFWHM <- function(peak, ...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculateFWHM: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)

  # function body  ---------------------------------------

  time <- peak@time
  sig <- peak@sig

  calc.fwhm <- function(sig,time) {

    # find the peak max
    peakmax <- max(sig)

    # determine the first timepoint that crosses the half line
    left.index <- c(which(sig - peakmax/2 > 0)[1] - 1,which(sig - peakmax/2 > 0)[1])
    right.index <- c(tail(which(sig - peakmax/2 > 0),1),tail(which(sig - peakmax/2 > 0),1) + 1)

    # if the leftmost left.index  is 0, which can happen if the peak value on the boundary is high, or if it's NA, which can happen if sig is all zeros, assign the peak boundary value to them:
    if (left.index[1] == 0 || is.na(left.index[1])) {

        t.left <- time[1]

      } else {

      # use linear interpolation to find the timepoint at half max.
      t.left <- (time[left.index[2]] - time[left.index[1]])/(sig[left.index[2]] - sig[left.index[1]])*(peakmax/2 - sig[left.index[1]]) + time[left.index[1]]
    }

    # if the rightmost right.index  is greater than length of time, which can happen if the peak value on the boundary is high, or if it's NA, which can happen if sig is all zeros, assign the peak boundary value to them:

    if (right.index[2] > length(time) || is.na(right.index[2])) {

      t.right <- tail(time,1)}

    else {
      t.right <- (time[right.index[2]] - time[right.index[1]])/(sig[right.index[2]] - sig[right.index[1]])*(peakmax/2 - sig[right.index[1]]) + time[right.index[1]]
    }

    # if t.left or t.right returns nothing, which can happen if the peak value on the boundary is high assign the peak boundary value to them:
    if (length(t.left) == 0) t.left <- time[1]
    if (length(t.right) == 0) t.right <- tail(time,1)

    # fwhm is the difference in time between the two timepoints that are crossed by the half max line
    fwhm <- t.right - t.left
    return(fwhm)
  }

  # calculate fwhm for each transition
  r.fwhm <- mapply(calc.fwhm,sig,data.frame(time))

  peak.fwhm <- round(mean(r.fwhm),digits = 4)

  # calculate fwhm to base ratio for each transition
  r.fwhm2base <- r.fwhm/(tail(time,1) - time[1])

  peak.fwhm2base <- round(mean(r.fwhm2base),digits = 4)

  # format output
  fwhm <- list(r.fwhm = r.fwhm, peak.fwhm = peak.fwhm, r.fwhm2base = r.fwhm2base, peak.fwhm2base = peak.fwhm2base)

  # return output
  return(fwhm)
}

#' Compute modality scores for transition peaks in a peak group of class peakObj.
#'
#' The function takes a peak group of class peakObj as input and calculates the modality score for the transition peaks in the peak group. Modality score is defined as the largest unexpected dip in the peak, normalized by peak height. For high quality peaks, the modality score is expected to be close to 0.
#'
#' @param peak A peak group of class peakObj
#' @param flatness.factor A numeric parameter between 0 and 1 that determines the sensitivity of the modality score to low levels of noise. To avoid high modality scores due to small levels of noise, near-flat ranges in the peak are artificially flattened before calculating the modality score. A range is defined as near-flat and thus flattened when the difference between intensities of adjacent time points is smaller than flatness.factor times the peak maximum intensity.
#'
#' @return A list with the following objects:
#'               r.modality: A numeric vector of modality scores for each transition in the peak group.
#'               peak.modality: Mean of the r.modality vector. This score represents the overall modality of all the transitions in the peak group.
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[196]]
#' peak.group.modality <- CalculateModality(peak)
#' PlotChromPeak(peak,transition.list = c("y6","y8"))

CalculateModality <- function(peak, flatness.factor = 0.05, ...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculateModality: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)

  # function body  ---------------------------------------

  time <- peak@time
  sig <- peak@sig

  calc.modality <- function(sig,time) {

    # find the differential of the peak
    diff.sig <- diff(sig)

    # any differences that are below the flatnessfactor of the maximum peak height are flattened.
    diff.sig[which(abs(diff.sig) < flatness.factor*max(abs(sig)))] <- 0

    # find the first and last timepoint where the differential changes sign
    first.fall <- head(which(diff.sig < 0),1)
    last.rise <- tail(which(diff.sig > 0),1)

    if (length(first.fall) == 0) first.fall <- length(time) + 1
    if (length(last.rise) == 0) last.rise <- -1

    # if first fall is after last rise, peak cannot be bi or multi-modal, so max.dip is set to 0. Otherwise it is set to the largest fall or rise between the first fall and last rise

    max.dip <- 0

    if (!is.na(first.fall) & !is.na(last.rise) & first.fall < last.rise) {
      max.dip <- max(abs(diff.sig[first.fall:last.rise]))
    }

    # The output is the maximum dip normalized by the peak height
    if (max(sig) == 0) {
      modality <- 0
    } else {
      modality <- max.dip/max(sig)
    }

    return(modality)
  }

  # calculate modality for each transition
  r.modality <- mapply(calc.modality,sig,data.frame(time))
  peak.modality <- round(mean(r.modality),digits = 4)

  # format output
  modality <- list(r.modality = r.modality, peak.modality = peak.modality)

  # return output
  return(modality)
}

#' Compute the sum of transition peaks in a peak group of class peakObj
#'
#' The function takes a peak group of class peakObj as input and calculates the sum of the peak intensities of transition peaks with identical labels. For example, all the endogenous transition peaks with labels of "light" will be added up to create the "sum.0.light" peak and all the spiked-in standard transition peaks with labels of "heavy" will be added up to create the "sum.0.heavy" peak.
#'
#' @param peak A peak group of class peakObj
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte: default is "heavy"
#'
#' @return A peak group of class peakObj, which contains two peaks: sum of light transition peaks and sum of heavy transition peaks in the input peak group
#'
#' @export
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' peak.sum <- CalculateTransitionSum(peak)
#' PlotChromPeak(peak.sum)

CalculateTransitionSum <- function(peak,endogenous.label = "light", standard.label = "heavy", ...) {

  # error and warning handling ---------------------------------------

  #   input errors
  error.input.format <- simpleError("CalculateTransitionSum: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)


  # function body  ---------------------------------------

  # separate columns of light and heavy isotopes

  endogenous.cols <- grepl(endogenous.label,colnames(peak@sig))
  standard.cols <- grepl(standard.label,colnames(peak@sig))

  # calculate the signal for the sum transition peak

  sig <- data.frame(rowSums(peak@sig[,endogenous.cols, drop = FALSE]), rowSums(peak@sig[,standard.cols, drop = FALSE]))

  colnames(sig) <- paste0("sum.0.",c(endogenous.label,standard.label))
  # calculate the area for the sum transition peak

  area <- c(sum(peak@area[endogenous.cols]),sum(peak@area[standard.cols]))
  names(area) <- paste0("sum.0.",c(endogenous.label,standard.label))

  # create the peak object

  peak.sum <- peakObj(
    time = peak@time,
    sig = sig,
    area = area
  )

  return(peak.sum)
}
