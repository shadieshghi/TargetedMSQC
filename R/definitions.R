# helper classes

matrixORnumeric <- setClassUnion("matrixORnumeric", c("matrix", "numeric"))

dataframeORnumeric <- setClassUnion("dataframeORnumeric", c("data.frame", "numeric"))


#' A custom class for storing single peaks or multiple peaks that belong to the same peptide and sample in a single object
#'
#' @slot time A numeric vector of the peak timepoints
#' @slot sig A numeric vector or dataframe of the peak signal intensities. If this is a data
#' @slot area A numeric vector of the AUC of peaks in the sig slot
#'
#' @exportClass peakObj
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' peak@@time
#' peak@@sig
#' peak@@area

peakObj <- setClass("peakObj",slots =
                       c(time = "numeric",
                         sig = "dataframeORnumeric",
                         area = "numeric")
                     )

# helper validate peak function

validPeakObject <- function(object) {

  # Number of time points in the time vector
  time.points = length(object@time)

  # Number of data points in each peak vector and number of transitions
  if(is.data.frame(object@sig)) {
    sig.points = nrow(object@sig)
    sig.transitions = ncol(object@sig)
  } else {
    sig.points = length(object@sig)
    sig.transitions = 1
  }


  # Check if the number of time points and data points in the signal are equal
  if(time.points != sig.points) {
      paste("Unequal time and sig lengths: ",time.points, ", ",sig.points,sep = "")
      return(FALSE)

  # Check if the number of transitions and signal columns are equal
  } else if(sig.transitions != length(object@area)) {
      paste("Unequal number of transitions in sig and area: ",sig.transitions, ", ",length(object@area),sep = "")
      return(FALSE)

    # The peak requires at least 4 time points for the peak QC functions to work.
  } else if (time.points < 4) {
      paste("Peak requires at least 4 time points. Time points: ", length(object@time), sep="")
      return(FALSE)
  } else return(TRUE)
}

#
#
setValidity("peakObj",validPeakObject)

#' Default S3 is.na method for peakObj
#' peakObj object is declared na if the time or sig portion of it is of length 0
#'
#
setMethod("is.na", "peakObj", definition =
            function(x) {
              if (length(x@time) * dim(x@sig)[1] * dim(x@sig)[2] * length(x@area) == 0) {
                return(TRUE)
              } else return(FALSE)
            }
)

#' Convert a dataframe with Times and Intensities columns into a peakObj object
#'
#' This is a helper function to format the output of Skyline into peakObj that are compatible with TargetedMSQC functions.
#'
#' @param df the dataframe with Times and Intensities columns for creating the peakObj object
#'
#' @return A peak group of class peakObj
#'
#' @importFrom tidyr unite_
#'

BuildPeakGroup <- function(df,...) {

  # purpose: converts the dataframe with Times and Intensities columns into a peakObj object
  #
  # args:
  #   df: the dataframe with Times and Intensities columns for creating the peakObj object.
  #
  # returns:
  #   the peakObj object


  # error and warning handling ---------------------------------------

  #   input errors: if the time and intensity vectors are na or of unequal length, an error is triggered

  error.input.format <- simpleError("BuildPeakGroup: input time and intensities should be non empty and of equal length")

  # function body  ---------------------------------------

  time <- as.numeric(strsplit(as.character(df$Times[1]),split = ",")[[1]])

  # if there are any peaks that are shorter than the time vector, pad the signal vector with zeros.
  tmp <- lapply(df$Intensities,function(x) as.numeric(strsplit(as.character(x),split = ",")[[1]]))

  # if time is longer than x, pad x with extra 0s. If x is longer than time, cut it to match the length of the time vector
  tmp <- lapply(tmp, function(x)
      {
        d <- length(time) - length(x)
        if (d > 0) {x <- c(x,rep(0,length(time) - length(x)))}
        if (d < 0) {x <- x[1:length(time)]}
        return(x)
      }
    )

  peak.sig <- data.frame(tmp)


  # input format error
  if (length(time) != nrow(peak.sig) || length(time)*dim(peak.sig)[1]*dim(peak.sig)[2] == 0) stop(error.input.format)
  # modifying the peak names
  identifying.cols <- colnames(df)[colnames(df) %in% c("FragmentIon","IsotopeLabelType","ProductCharge")]
  if (length(identifying.cols) > 0) {
    columns <- unite_(df,"PeakName",identifying.cols,sep = ".")$PeakName
    colnames(peak.sig) <- columns
  }

 # Total area
#  area <- as.numeric(df$TotalArea)
#  names(area) <- columns

  # calculate peak area using the trapozoidal approximation
  area <- sapply(peak.sig,function(x) pracma::trapz(time,x))


 # identifying missing peaks
 # peaks that do not have their pair

 single.peaks <- columns[which(mapply(function(x,columns) sum(grepl(substr(x,1,regexpr("\\.[^\\.]*$", x)),columns)),columns,list(columns)) != 2)]

 all.isotopes <- unique(df$IsotopeLabelType)

 if (length(single.peaks) > 0) {

   # isotope label of the missing peaks
   missing.isotopes <- mapply(function(x,all.isotopes) all.isotopes[which(all.isotopes != substr(x,regexpr("\\.[^\\.]*$", x) + 1,nchar(x)))],single.peaks,list(all.isotopes))

   # names of the missing peaks
   missing.peaks <- sapply(single.peaks,function(x) substr(x,1,regexpr("\\.[^\\.]*$", x)))

   # colnames of the missing peaks that are to be imputated by zero
   imputed.peaks <- paste(missing.peaks,missing.isotopes,sep = "")


   # imputing the peak.sig and total area with vectors of 0s for missing peaks
   for (column in imputed.peaks) {
     peak.sig[,column] <- rep(0,nrow(peak.sig))
     area[column] <- 0
   }
 }


 # create the output peakObj object
  peakgroup <- peakObj(
    time = time,
    sig = peak.sig,
    area = area
  )
  return(peakgroup)
}

