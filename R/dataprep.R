## This file contains the functions that clean up the input data and calculate the features associated with peaks and peak groups.

#' Convert the exported skyline chromatogram and peak boundary reports to a dataframe format compatible with TargetedMSQC functions.
#'
#' The function takes the directories that contain chromatogram and peak boundary files of Skyline documents and merges them. The rows with NA values for peak boundaries, rows corresponding to peptides with fewer than 3 transitions and transitions missing a standard isotope pair are separated into a dataframe that is returned by the function in output$removed and excluded from subsequence QC steps. The remainder of the rows are returned by the function in output$data for subsequence analysis. Transition peaks of each peptide in each run are formed into peak groups of custom class peakObj and stored in the ChromGroup column. The ApplyPeakBoundary function is applied to chromatograms to restrict them to the provided peak boundaries and the results are stored in the PeakGroup column. Additionally, rows are added for the "sum" transition, which are is the sum of peptides transitions with the same isotope label in each sample.
#'
#' @param chromatogram.path Path to the directory containing the .tsv files of the peak chromatograms. For each Skyline document, this file is exported from Skyline through File > Export > Chromatograms. Here, check runs of interest and include Precursors, Products, Base Peaks and TICs. Each chromatogram .tsv file corresponds to a single Skyline document, which may contain any number of runs. Multiple chromatogram files, corresponding to multiple Skyline documents can be copied into the chromatogram.path directory. For each chromatogram file in this folder, there should be a peak boundary file with an identical name in peak.boundary.path directory.
#' @param peak.boundary.path  Path to the directory containing the .csv files of the peak boundaries. For each Skyline document, this file is exported from Skyline through File > Export > Report. Here, select Peak Boundaries. Each peak boundary .csv file corresponds to a single Skyline document, which may contain any number of runs. Multiple peak boundary files, corresponding to multiple Skyline documents can be copied into the peak.boundary.path directory. For each peak boundary file in this folder, there should be a peak chromatogram file with an identical name in chromatogram.path directory.
#' @param labkey.url.base URL of the labkey server, if data is being imported from labkey (labkey = TRUE). This feature has not been implemented yet.
#' @param labkey.url.path Path to the directory containing the Skyline documents on the labkey server, if data is being imported from labkey (labkey = TRUE). This feature has not been implemented yet.
#' @param parallel Logical parameter to determine whether the function should run in parallel. This feature has not been implemented yet. Set parallel = FALSE.
#' @param labkey Logical parameter to indicate whether data will be imported from labkey. This feature has not been implemented yet. Set labkey = FALSE.
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte: default is "heavy"
#' @param iRT.list List of iRT standards used in the experiment. These peptides will be removed from the training set.
#'
#' @return A list with the following objects:
#'               data: A dataframe that contains all rows of the input data with assigned peak boundaries. This dataframe also includes columns for peakObj objects created for each peak group.
#'               removed: A dataframe that contains the rows with missing peak boundary values, too few transitions and missing isotope pairs. These rows are removed from downstream feature extraction and QC analysis.
#'
#' @export
#'
#' @importFrom tidyr spread
#' @importFrom data.table setDT
#' @importFrom data.table rbindlist
#' @importFrom data.table rbindlist
#'
#' @examples
#'
#' extdata.path <- system.file("extdata",package = "TargetedMSQC")
#' project.folder.name <- "CSF_Panel"
#' project.path <- file.path(extdata.path,project.folder.name)
#' chromatogram.path <- file.path(project.path,"Chromatograms")
#' peak.boundary.path <- file.path(project.path,"Peak_boundary")
#' data <- CleanUpChromatograms(chromatogram.path = chromatogram.path,
#'                             peak.boundary.path = peak.boundary.path,
#'                             endogenous.label = "light",
#'                             standard.label = "heavy",
#'                             iRT.list = c("LGGNETQVR","AGGSSEPVTGLADK","VEATFGVDESANK","YILAGVESNK","TPVISGGPYYER","TPVITGAPYYER","GDLDAASYYAPVR","DAVTPADFSEWSK","TGFIIDPGGVIR","GTFIIDPAAIVR","FLLQFGAQGSPLFK"))


CleanUpChromatograms <- function(chromatogram.path = NULL, peak.boundary.path = NULL, labkey.url.base = NULL, labkey.url.path = NULL, parallel = FALSE, labkey = FALSE, endogenous.label = "light", standard.label = "heavy" , iRT.list = c("LGGNETQVR","AGGSSEPVTGLADK","VEATFGVDESANK","YILAGVESNK","TPVISGGPYYER","TPVITGAPYYER","GDLDAASYYAPVR","DAVTPADFSEWSK","TGFIIDPGGVIR","GTFIIDPAAIVR","FLLQFGAQGSPLFK","LGGNEQVTR","GAGSSEPVTGLDAK","VEATFGVDESNAK","YILAGVENSK","TPVISGGPYEYR","TPVITGAPYEYR","DGLDAASYYAPVR","ADVTPADFSEWSK","GTFIIDPGGVIR","GTFIIDPAAVIR","LFLQFGAQGSPFLK"), ...) {

  # parallel = TRUE indicates that function should run in parallel. Currently, the code to use parallel is not implemened/tested and parallel is automatically set to FALSE.
  if (parallel) {
    warning('The code to handle parallel = TRUE has not been fully implemented and tested yet. parallel is set to FALSE.')
    parallel = FALSE
  }


  # labkey = TRUE indicates that data is going to be read from labkey using the Rlabkey API. Currently, the code to read data from labkey has not been developed and labkey is automatically set to FALSE.
  if (labkey) {
    warning('The code to handle Labkey = TRUE has not been implemented yet. Labkey is set to FALSE.')
    labkey = FALSE
  }
  if (!labkey) {

    # error and warning handling ---------------------------------------

    #   input errors: if the input are empty

    error.input.format1 = simpleError("CleanUpChromatograms: Input chromatogram file should not be empty")
    error.input.format2 = simpleError("CleanUpChromatograms: Input peak boundary file should not be empty")
    error.input.format3 = simpleError("CleanUpChromatograms: For each chromatograph file in chromatogram.path, there should be a peak boundary file in peak.boundary.path with an identical file name and vice versa")

    # chromgram.path and peak.boundary.path can only be NULL if the code is being run on labkey.

    if (is.null(chromatogram.path)) stop(error.input.format1)
    if (is.null(peak.boundary.path)) stop(error.input.format2)


    # check the names of files in chromatogram.path and peak.boundary.path

    chromatogram.files <- sort(list.files(chromatogram.path,pattern = "\\.tsv"))
    chromatogram.files <- unlist(lapply(chromatogram.files,function(x) substr(x,1,regexpr("\\.[^\\.]*$", x)[[1]] - 1)))

    peak.boundary.files <- sort(list.files(peak.boundary.path,pattern = "\\.csv"))
    peak.boundary.files <- unlist(lapply(peak.boundary.files,function(x) substr(x,1,regexpr("\\.[^\\.]*$", x)[[1]] - 1)))

    if (!identical(chromatogram.files,peak.boundary.files)) stop(error.input.format3)

    # read the chromatrogram files

    chrom <- setDT(rbindlist(lapply(chromatogram.files,
                                    function(chromatogram.file) {
                                      chromatogram.file.path <- file.path(chromatogram.path,paste0(chromatogram.file,".tsv"))
                                      tmp <- unique(read.table(chromatogram.file.path,sep = "\t", header = TRUE))
                                      tmp$File <- chromatogram.file
                                      return(tmp)
                                }
                        )
      )
    )


    # convert isotopelabeltypes to light and heavy. Only these two values are acceptable:

    chrom$IsotopeLabelType <-  tolower(chrom$IsotopeLabelType)
    chrom$IsotopeLabelType <- factor(chrom$IsotopeLabelType,levels = c(endogenous.label,standard.label))

    # check if each FragmentIon/Peptide pair is represented by both endogenous and standard labels in each sample

    tmp <- chrom %>%
      select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType, TotalArea) %>%
      spread(IsotopeLabelType,TotalArea)

    tmp <- as.data.frame(tmp)

    tmp$RemoveIsotopePair <- FALSE
    tmp[is.na(tmp[,endogenous.label]) | is.na(tmp[,standard.label]),"RemoveIsotopePair"] <- TRUE

    chrom <- chrom %>%
      left_join(tmp[,!colnames(tmp) %in% c(endogenous.label,standard.label,"<NA>")],by = c("File","FileName","PeptideModifiedSequence","PrecursorCharge","FragmentIon","ProductCharge"))

    # mark the rows with any isotopelabeltype except endogenous. and standard.label are marked for removal
    chrom[is.na(chrom$IsotopeLabelType),"RemoveIsotopePair"] <- TRUE

    # read the peak boundary file

    peak.boundary <- setDT(rbindlist(lapply(peak.boundary.files,
                                            function(peak.boundary.file) {
                                              peak.boundary.file.path <- file.path(peak.boundary.path,paste0(peak.boundary.file,".csv"))
                                              tmp <- unique(read.csv(peak.boundary.file.path,sep = ",", header = TRUE))
                                              tmp$File <- peak.boundary.file
                                              return(tmp)
                                    }
                        )
      )
    )

    peak.boundary <- as.data.frame(peak.boundary)

    if (nrow(chrom) == 0) stop(error.input.format1)
    if (nrow(peak.boundary) == 0) stop(error.input.format2)


    # change the column names of the peak boundary
    colnames(peak.boundary) = gsub("\\.","",colnames(peak.boundary))

    # change the column names of the chrom
    colnames(chrom) = gsub("\\.","",colnames(chrom))


    # change the class of time columns to numeric
    numeric.cols = c("MinStartTime","MaxEndTime")

    peak.boundary[,numeric.cols] = apply(peak.boundary[,numeric.cols],2,function(x) as.numeric(as.character(x)))


    # mark the peaks with missing peak boundaries or the ones with multiple peak boundaries for removal. Note that sometimes skyline exports multiple lines for a peak, where one line is NA. In these cases, I keep the non-NA line and remove the NA line as redundant.
    peak.boundary <- unique(peak.boundary)
    peak.boundary$RemovePeakBoundary <- FALSE
    peak.boundary$Redundant <- FALSE

    peak.boundary <- peak.boundary %>% group_by(FileName,File,PeptideModifiedSequence) %>% mutate(n = n())
    peak.boundary[(is.na(peak.boundary$MinStartTime) | is.na(peak.boundary$MaxEndTime)) & (peak.boundary$n > 1),"Redundant"] <- TRUE
    peak.boundary <- peak.boundary[!peak.boundary$Redundant,] %>% select(-Redundant)

    peak.boundary <- peak.boundary %>% group_by(FileName,File,PeptideModifiedSequence) %>% mutate(n = n())
    peak.boundary[is.na(peak.boundary$MinStartTime) | is.na(peak.boundary$MaxEndTime) | (peak.boundary$n > 1),"RemovePeakBoundary"] <- TRUE

    peak.boundary <- peak.boundary %>% select(-n)

    # merge peak boundary and chromatrogram
    data <- chrom %>% left_join(peak.boundary,by = c("File","FileName","PeptideModifiedSequence","PrecursorCharge"))

    data$Remove <- data$RemoveIsotopePair | data$RemovePeakBoundary

  }

  # Rows where peak boundaries are NA are separated into a data frame and returned as output.

  removed <- data %>%
    filter(Remove | PeptideModifiedSequence %in% iRT.list) %>%
    select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,MinStartTime,MaxEndTime,RemoveIsotopePair,RemovePeakBoundary,Remove)

  data <- data %>% filter(!Remove & !(PeptideModifiedSequence %in% iRT.list)) %>%
    select(FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,TotalArea,Times,Intensities,MinStartTime,MaxEndTime,File)

  #
  if (nrow(data) == 0) {
    warnings("data is not appropriately formatted for TargetedMSQC. Please check the RemoveIsotopePair and RemovePeakBoundary columns in output$removed.")
  } else {
    # substitute NA values in the TotalArea column with 0.

    if (sum(is.na(data$TotalArea)) > 0) {
      warnings("NA values are detected in the peak area column of the input and replaced with 0")
      data$TotalArea[is.na(data$TotalArea)] = 0
    }


    # build the chromatogram group for each (File,FileName,peptide,precursorcharge) trio

    peak.groups = data %>%
      group_by(File,FileName,PeptideModifiedSequence,PrecursorCharge) %>%
      do(ChromGroup = BuildPeakGroup(.))

    data = unique(merge(data,peak.groups))

    # build the peak group for each (File,FileName,peptide,precursorcharge) trio: the only difference with the chromatogram group is that the peak boundaries are applied.

    if (parallel) {
      data$PeakGroup =  mcmapply(ApplyPeakBoundary,data$ChromGroup,boundary = as.list(data.frame(t(cbind(data$MinStartTime,data$MaxEndTime)))),
                                 mc.silent = FALSE, mc.cores = detectCores() - 1)

    } else {
      data$PeakGroup =  mapply(ApplyPeakBoundary,data$ChromGroup,boundary = as.list(data.frame(t(cbind(data$MinStartTime,data$MaxEndTime)))))
    }

    # remove those rows where the peak is not a proper peak object. This can happen if the peak boundaries are too narrow and the peak has only 3 or fewer points across it.

    removed.na.peaks <- data %>%
      filter(is.na(PeakGroup)) %>%
      select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,MinStartTime,MaxEndTime)

    if (nrow(removed.na.peaks) > 0) {

      removed.na.peaks$RemoveIsotopePair <- FALSE
      removed.na.peaks$RemovePeakBoundary <- TRUE
      removed.na.peaks$Remove <- TRUE

    }

    removed <- rbind(removed,removed.na.peaks)


    data <- data %>% filter(!is.na(PeakGroup))

    # add the peak area column to the data frame
    data$PeakArea <- NULL
    for (i in 1:nrow(data)) {
      tmp <- paste(data$FragmentIon[[i]],data$ProductCharge[[i]],data$IsotopeLabelType[[i]],sep = ".")
      data$PeakArea[i] <- data$PeakGroup[[i]]@area[[tmp]]
    }

    # calculate the sum of AUCs of individual transitions and assign to SumArea
    data = data %>%
      select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,PeakArea) %>%
      group_by(PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType, FileName,File) %>%
      summarise(SumArea = sum(PeakArea)) %>%
      right_join(data,by = c("PeptideModifiedSequence","PrecursorCharge", "IsotopeLabelType","FileName","File"))

    #
    # calculate sum of transitions with the FragmentIon name of "sum". Add the new "sum" transition rows to the dataframe
    tmp = data %>%
      select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,PeakGroup,MinStartTime,MaxEndTime) %>%
      group_by(File,FileName,PeptideModifiedSequence,PrecursorCharge) %>%
      summarise(Peak = PeakGroup[1],
                MinStartTime = MinStartTime[1],
                MaxEndTime = MaxEndTime[1])

    if (parallel) {
      tmp$PeakGroup = mcmapply(function(x,endogenous.label,standard.label) CalculateTransitionSum(x,endogenous.label,standard.label),
                               tmp$Peak,endogenous.label = endogenous.label,standard.label = standard.label,
                               mc.silent = FALSE, mc.cores = detectCores() - 1)

    } else {
      tmp$PeakGroup = mapply(function(x,endogenous.label,standard.label)
        CalculateTransitionSum(x,endogenous.label,standard.label),
        tmp$Peak,endogenous.label = endogenous.label,standard.label = standard.label)
    }


    tmp = tmp %>% select(-Peak)

    # for "sum" transitions assign default values to appropriate columns
    tmp = tmp %>%
      mutate(FragmentIon = as.factor("sum")) %>%
      mutate(ProductCharge = as.integer(0)) %>%
      mutate(IsotopeLabelType = as.factor(endogenous.label)) %>%
      mutate(Times = as.factor(NA)) %>%
      mutate(Intensities = as.factor(NA)) %>%
      mutate(ChromGroup = as.factor(NA))

    # peak area for sum transitions is the sum of the PeakArea of all transitions. This value is calculated for endogenous and standard sums separately

    if (parallel) {
      tmp$PeakArea = parSapply(cl,tmp$PeakGroup, function(peak) {
        area = peak@area
        return(area[paste0("sum.0.",endogenous.label)])
      })
    } else {
        tmp$PeakArea = sapply(tmp$PeakGroup, function(peak) {
          area = peak@area
          return(area[paste0("sum.0.",endogenous.label)])
        })
    }

    # SumArea is the sum of indivodual transitions for each peak and sample and is equal to TotalAre of "sum" transitions
    tmp$SumArea = tmp$PeakArea

    # append the new rows to data frame
    data <- data %>% select(-TotalArea)
    tmp = data.frame(tmp[,colnames(data)])
    data = rbind.data.frame(data,tmp)

    # do the same for heavy sum isotopes. The only difference is in the IsotopeLabelType and PeakArea
    tmp$IsotopeLabelType <- standard.label

    if (parallel) {
      tmp$PeakArea <- parSapply(cl,tmp$PeakGroup, function(peak) {
        area <- peak@area
        return(area[paste0("sum.0.",standard.label)])
      })
    } else {
      tmp$PeakArea <- sapply(tmp$PeakGroup, function(peak) {
        area = peak@area
        return(area[paste0("sum.0.",standard.label)])
      })
    }

    # SumArea is the sum of indivodual transitions for each peak and sample and is equal to TotalAre of "sum" transitions
    tmp$SumArea <- tmp$PeakArea

    # append the sum transitions for heavy isotopes

    data <- rbind.data.frame(data,tmp) %>% ungroup()
    data <- unique(data)
  }
  # deactivate the cluster
  if (parallel) stopCluster(cl)

  # return output
  return(list(data = data,removed = removed))
}


#' Calculate an ensemble of QC features/metrics for a dataset that has been cleaned up by CleanUpChromatograms.
#'
#' The function takes the output of CleanUpChromatograms as input and for each transition pair of each peptide, calculates a number of QC metrics that are used as features for training a predictive peak QC model.
#'
#' @param data A dataframe that contains identifier columns and peak groups of class peakObj for each transition peak. This dataframe is the output of CleanUpChromatograms (output$data).
#' @param parallel Logical parameter to determine whether the function should run in parallel. This feature has not been implemented yet. Set parallel = FALSE.
#' @param blanks A dataframe with two columns of File and FileName. The FileName column contains the names of blank runs for each Skyline File in the File column. The values in File and FileName should be consistent with corresponding columns in input data. Example: If Skyline document "SkylineFile1" includes three blank runs of "Blank1", "Blank2" and "Blank3", blanks <- data.frame(File = "SkylineFile1", FileName = c("Blank1","Blank2","Blank3")). The function uses the blank runs to estimate the peak intensity at LLOQ for each pepide and transition. If blanks = NA, the function uses the intensity.threshold by default as approximation of intensity at LLOQ for all peaks. This threshold can be adjusted.
#' @param intensity.threshold The threshold that the function uses as an approximation of the intensity at LLOQ for all peptides. If blank samples are not provided to estimate the intensity at LLOQ, this value is used as an approximation to determine transition peaks that are below limit of quantitation.
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte: default is "heavy"
#' @param export.features Logical parameter to indicate whether calculated features should be exported as a .csv file.
#' @param feature.path Path to the directory, where the features should be saved if export.features = TRUE.
#'
#'
#' @return A list with the following objects:
#'               features: A dataframe with columns that contain the calculated QC features for each transition pair in the input data.
#'
#' @export
#'
#' @importFrom tidyr spread
#' @importFrom data.table dcast
#' @importFrom data.table setDT
#'
#' @examples
#'
#' data.features <- ExtractFeatures(data = data.CSF$data,
#'                                 export.features = FALSE,
#'                                 intensity.threshold = 1000)

ExtractFeatures <- function(data, parallel = FALSE, blanks = NA, intensity.threshold = 100000, endogenous.label = "light", standard.label = "heavy", export.features = FALSE, feature.path = "", ...) {

  # parallel = TRUE indicates that function should run in parallel. Currently, the code to use parallel is not implemened/tested and parallel is automatically set to FALSE.
  if (parallel) {
    warning('The code to handle parallel = TRUE has not been fully implemented and tested yet. parallel is set to FALSE.')
    parallel = FALSE
  }


  # if the dataset includes blank sample, use those to estimate the intensity at lloq by mean(intensity) + 10*sd(intensity) and mark what measurements are above the lloq level. If blanks are not available, just use the default value of 1e6 as your threshold. Based on looking at a dataset ~90% of peptides have an lloq.intesnity below 1e6 for individual fragment ions.
  if (!is.na(blanks)) {
    blank = data %>%
      filter(interaction(File,FileName) %in% interaction(blanks$File,blanks$FileName)) %>%
      group_by(File,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,FragmentIon,ProductCharge) %>%
      summarise(lloq.intensity = mean(PeakArea) + 10*sd(PeakArea))
    data = data %>%
      left_join(blank,by = c("File","PeptideModifiedSequence","PrecursorCharge","IsotopeLabelType","FragmentIon","ProductCharge"))
   data$aboveThreshold = data$PeakArea > data$lloq.intensity
   data = data %>% select(-lloq.intensity)
  }
  else {
    data$aboveThreshold = data$PeakArea > intensity.threshold
  }

  if (parallel) {
    # Calculate the number of cores
    no_cores <- detectCores() - 1

    # Initiate cluster
    cl <- makeCluster(no_cores)

    # export the unknown functions to individual processors
   # clusterExport(cl,"corr.test")

  }

  if (parallel) {
    # calculate jaggedness for each peak group
    tmp = parSapply(cl,data$PeakGroup,CalculatePeakJaggedness)

    data$PeakGroupJaggedness.r = tmp[1,]
    data$PeakGroupJaggedness.m = as.numeric(tmp[2,])

    # calculate symmetry for each peak group
    tmp = parSapply(cl,data$PeakGroup,CalculatePeakSymmetry)

    data$PeakGroupSymmetry.r = tmp[1,]
    data$PeakGroupSymmetry.m = as.numeric(tmp[2,])


    # export the unknown functions to individual processors
   # clusterExport(cl,"corr.test")

    # calculate shape similarity for each peak group
    tmp = parSapply(cl,data$PeakGroup,CalculatePeakShapeSimilarity)

    data$PeakGroupSimilarity.r = tmp[1,]
    data$PeakGroupSimilarity.m = as.numeric(tmp[2,])

    # export the unknown functions to individual processors
    #  clusterExport(cl,"CalculateElutionShift") # toghiess: if any error, uncomment this, if not remove the lie

    # calculate the elution shift for each peak group
    tmp = parSapply(cl,data$PeakGroup,CalculatePeakElutionShift)

    data$PeakGroupShift.r = tmp

    # calculate fwhm and fwhm2base ratio for each peak group
    tmp = parSapply(cl,data$PeakGroup,CalculateFWHM)

    data$PeakGroupFWHM.r = tmp[1,]
    data$PeakGroupFWHM.m = as.numeric(tmp[2,])
    data$PeakGroupFWHM2base.r = tmp[3,]
    data$PeakGroupFWHM2base.m = as.numeric(tmp[4,])

    # calculate modality ratio for each peak group
    tmp = parSapply(cl,data$PeakGroup,CalculateModality)

    data$PeakGroupModality.r = tmp[1,]
    data$PeakGroupModality.m = as.numeric(tmp[2,])

    # calculate max intensity for each transition
    tmp = parSapply(cl,data$PeakGroup,CalculatePeakMaxIntensity)
    data$TransitionMaxIntensity = mcmapply(function(r,trn) r[trn],
                                         tmp,
                                         trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                         mc.silent = FALSE, mc.cores = detectCores() - 1)

    # calculate max intensity at boundary for each transition
    tmp = parSapply(cl,data$PeakGroup,CalculateMaxBoundaryIntensity)
    data$TransitionMaxBoundaryIntensity = mcmapply(function(r,trn) r[trn],
                                           tmp,
                                           trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                           mc.silent = FALSE, mc.cores = detectCores() - 1)



    # extract jaggedness for each transition from the jaggedness.r matrix
    data$TransitionJaggedness = mcmapply(function(r,trn) r[trn],
                                         data$PeakGroupJaggedness.r,
                                         trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                         mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract symmetry for each transition from the symmetry.r matrix
    data$TransitionSymmetry = mcmapply(function(r,trn) r[trn],
                                       data$PeakGroupSymmetry.r,
                                       trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                       mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract shift for each transition from shift.r matrix using
    data$TransitionShift = mcmapply(function(r,trn) r[trn,trn],
                                  data$PeakGroupShift.r,
                                  trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                  mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract similarity between light and heavy isotope for each transition from the similarity.r matrix
    data$PairSimilarity = mcmapply(function(r,trn.endg,trn.std) r[trn.endg,trn.std],
                                   data$PeakGroupSimilarity.r,
                                   trn.endg = paste(data$FragmentIon,data$ProductCharge,endogenous.label,sep = "."),
                                   trn.std = paste(data$FragmentIon,data$ProductCharge,standard.label,sep = "."),
                                   mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract elution shift between light and heavy isotope for each transition from the shift.r matrix
    data$PairShift = mcmapply(function(r,trn.endg,trn.std) r[trn.endg,trn.std],
                              data$PeakGroupShift.r,
                              trn.endg = paste(data$FragmentIon,data$ProductCharge,endogenous.label,sep = "."),
                              trn.std = paste(data$FragmentIon,data$ProductCharge,standard.label,sep = "."),
                              mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract fwhm2base ratio for each transition from the fwhm2base.r matrix
    data$TransitionFWHM2base = mcmapply(function(r,trn) r[trn],
                                        data$PeakGroupFWHM2base.r,
                                        trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                        mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract fwhm  for each transition from the fwhm.r matrix
    data$TransitionFWHM = mcmapply(function(r,trn) r[trn],
                                   data$PeakGroupFWHM.r,
                                   trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                   mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract modality  for each transition from the modality.r matrix
    data$TransitionModality = mcmapply(function(r,trn) r[trn],
                                       data$PeakGroupModality.r,
                                       trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."),
                                       mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average jaggedness  for all light (heavy) transitions in each peak group from the jaggedness.r matrix
    data$IsotopeJaggedness = mcmapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                      data$PeakGroupJaggedness.r,
                                      isotope = data$IsotopeLabelType,
                                      mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average symmetry  for all light (heavy) transitions in each peak group from the symmetry.r matrix
    data$IsotopeSymmetry = mcmapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                    data$PeakGroupSymmetry.r,
                                    isotope = data$IsotopeLabelType,
                                    mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average similarity  for all light (heavy) transitions in each peak group from the similarity.r matrix
    data$IsotopeSimilarity = mcmapply(function(r,isotope) round(mean(abs(r[grep(isotope,colnames(r)),grep(isotope,rownames(r))])),digits = 4),
                                      data$PeakGroupSimilarity.r,
                                      isotope = data$IsotopeLabelType,
                                      mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average elution shift  for all light (heavy) transitions in each peak group from the shift.r matrix
    data$IsotopeShift = mcmapply(function(r,isotope) round(mean(abs(r[grep(isotope,colnames(r)),grep(isotope,rownames(r))])),digits = 4),
                                 data$PeakGroupShift.r,
                                 isotope = data$IsotopeLabelType,
                                 mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average fwhm2base ratio  for all light (heavy) transitions in each peak group from the fwhm2base.r matrix
    data$IsotopeFWHM2base = mcmapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                     data$PeakGroupFWHM2base.r,
                                     isotope = data$IsotopeLabelType,
                                     mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average fwhm  for all light (heavy) transitions in each peak group from the fwhm.r matrix
    data$IsotopeFWHM = mcmapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                data$PeakGroupFWHM.r,
                                isotope = data$IsotopeLabelType,
                                mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average modality  for all light (heavy) transitions in each peak group from the modality.r matrix
    data$IsotopeModality = mcmapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                    data$PeakGroupModality.r,
                                    isotope = data$IsotopeLabelType,
                                    mc.silent = FALSE, mc.cores = detectCores() - 1)

    # extract average co-elution  for all transitions in each peak group from the shift.r matrix
    data$PeakGroupShift = parSapply(cl,data$PeakGroupShift.r, function(r) round(mean(diag(abs(r))),digits = 4))


  }
  else {


    # calculate jaggedness for each peak group
    tmp = sapply(data$PeakGroup,CalculatePeakJaggedness)

    data$PeakGroupJaggedness.r = tmp[1,]
    data$PeakGroupJaggedness.m = as.numeric(tmp[2,])

    # calculate symmetry for each peak group
    tmp = sapply(data$PeakGroup,CalculatePeakSymmetry)

    data$PeakGroupSymmetry.r = tmp[1,]
    data$PeakGroupSymmetry.m = as.numeric(tmp[2,])

    # calculate shape similarity for each peak group
    tmp = sapply(data$PeakGroup,CalculatePeakShapeSimilarity)

    data$PeakGroupSimilarity.r = tmp[1,]
    data$PeakGroupSimilarity.m = as.numeric(tmp[2,])

    # calculate the elution shift for each peak group using the CalculatePeakElutionShift method
    tmp = sapply(data$PeakGroup,CalculatePeakElutionShift)

    data$PeakGroupShift.r = tmp

    # calculate fwhm and fwhm2base ratio for each peak group
    tmp = sapply(data$PeakGroup,CalculateFWHM)

    data$PeakGroupFWHM.r = tmp[1,]
    data$PeakGroupFWHM.m = as.numeric(tmp[2,])
    data$PeakGroupFWHM2base.r = tmp[3,]
    data$PeakGroupFWHM2base.m = as.numeric(tmp[4,])

    # calculate modality for each peak group
    tmp = sapply(data$PeakGroup,CalculateModality)

    data$PeakGroupModality.r = tmp[1,]
    data$PeakGroupModality.m = as.numeric(tmp[2,])

    # calculate max intensity for each transition
    tmp = sapply(data$PeakGroup,CalculatePeakMaxIntensity)
    data$TransitionMaxIntensity = mapply(function(r,trn) r[trn],
                                         tmp,
                                         trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))


    # calculate max intensity at boundary for each transition
    tmp = sapply(data$PeakGroup,CalculateMaxBoundaryIntensity)
    data$TransitionMaxBoundaryIntensity = mapply(function(r,trn) r[trn],
                                         tmp,
                                         trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))

    # calculate max intensity at boundary normalized by max intensity for each transition
    data$TransitionMaxBoundaryIntensityNormalized = data$TransitionMaxBoundaryIntensity/data$TransitionMaxIntensity

    # impute TransitionMaxBoundaryIntensityNormalized values of NaN to 0. This happens transition is all zeros.
    data$TransitionMaxBoundaryIntensityNormalized[!is.finite(data$TransitionMaxBoundaryIntensityNormalized)] <- 0

    # extract jaggedness for each transition from the jaggedness.r matrix
    data$TransitionJaggedness = mapply(function(r,trn) r[trn],
                                       data$PeakGroupJaggedness.r,
                                       trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))

    # extract symmetry  for each transition from the symmetry.r matrix
    data$TransitionSymmetry = mapply(function(r,trn) r[trn],
                                     data$PeakGroupSymmetry.r,
                                     trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))

    # extract shift for each transition from shift.r matrix using
    data$TransitionShift = mapply(function(r,trn) r[trn,trn],
                                     data$PeakGroupShift.r,
                                     trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))


    # extract similarity between light and heavy isotope for each transition from the similarity.r matrix
    data$PairSimilarity = mapply(function(r,trn.endg,trn.std) r[trn.endg,trn.std],
                                 data$PeakGroupSimilarity.r,
                                 trn.endg = paste(data$FragmentIon,data$ProductCharge,endogenous.label,sep = "."),
                                 trn.std = paste(data$FragmentIon,data$ProductCharge,standard.label,sep = "."))


    # extract elution shift between light and heavy isotope for each transition from the shift.r matrix
    data$PairShift = mapply(function(r,trn.endg,trn.std) r[trn.endg,trn.std],
                            data$PeakGroupShift.r,
                            trn.endg = paste(data$FragmentIon,data$ProductCharge,endogenous.label,sep = "."),
                            trn.std = paste(data$FragmentIon,data$ProductCharge,standard.label,sep = "."))


    # extract fwhm2base  for each transition from the fwhm2base.r matrix
    data$TransitionFWHM2base = mapply(function(r,trn) r[trn],
                                      data$PeakGroupFWHM2base.r,
                                      trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))

    # extract fwhm  for each transition from the fwhm.r matrix
    data$TransitionFWHM = mapply(function(r,trn) r[trn],
                                 data$PeakGroupFWHM.r,
                                 trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))

    # extract modality  for each transition from the modality.r matrix
    data$TransitionModality = mapply(function(r,trn) r[trn],
                                     data$PeakGroupModality.r,
                                     trn = paste(data$FragmentIon,data$ProductCharge,data$IsotopeLabelType,sep = "."))

    # extract average jaggedness  for all light (heavy) transitions in each peak group from the jaggedness.r matrix
    data$IsotopeJaggedness = mapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                    data$PeakGroupJaggedness.r,
                                    isotope = data$IsotopeLabelType)

    # extract average symmetry  for all light (heavy) transitions in each peak group from the symmetry.r matrix
    data$IsotopeSymmetry = mapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                  data$PeakGroupSymmetry.r,
                                  isotope = data$IsotopeLabelType)

    # extract average similarity  for all light (heavy) transitions in each peak group from the similarity.r matrix
    data$IsotopeSimilarity = mapply(function(r,isotope) round(mean(abs(r[grep(isotope,colnames(r)),grep(isotope,rownames(r))])),digits = 4),
                                    data$PeakGroupSimilarity.r,
                                    isotope = data$IsotopeLabelType)

    # extract average co-elution  for all light (heavy) transitions in each peak group from the shift.r matrix
    data$IsotopeShift = mapply(function(r,isotope) round(mean(diag(matrix(abs(r[grep(isotope,colnames(r)),grep(isotope,rownames(r))])))),digits = 4),
                               data$PeakGroupShift.r,
                               isotope = data$IsotopeLabelType)

    # extract average fwhm2base  for all light (heavy) transitions in each peak group from the fwhm2base.r matrix
    data$IsotopeFWHM2base = mapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                   data$PeakGroupFWHM2base.r,
                                   isotope = data$IsotopeLabelType)

    # extract average fwhm  for all light (heavy) transitions in each peak group from the fwhm.r matrix
    data$IsotopeFWHM = mapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                              data$PeakGroupFWHM.r,
                              isotope = data$IsotopeLabelType)

    # extract average modality  for all light (heavy) transitions in each peak group from the modality.r matrix
    data$IsotopeModality = mapply(function(r,isotope) round(mean(r[grep(isotope,names(r))]),digits = 4),
                                  data$PeakGroupModality.r,
                                  isotope = data$IsotopeLabelType)
    # extract average co-elution  for all transitions in each peak group from the shift.r matrix
    data$PeakGroupShift = sapply(data$PeakGroupShift.r, function(r) round(mean(diag(abs(r))),digits = 4))

  }


  # calculate the center for each peak which is equivalent to RT for the peak
  data$PeakCenter = (data$MaxEndTime + data$MinStartTime)/2



  # calculate the ratio consistency between the endogenous and standard for each sample and peptide and tranistion:
  # PairRatioConsistency =  abs(endogenous - standard)/standard
  # Area2SumRatio is the ratio of eahc transition to sum of all transitions in each sample for endogenous and standard isotopes

  data = data %>% ungroup()

  data = data %>%
    mutate(Area2SumRatio = PeakArea/SumArea)

  # impute Area2SumRaio values of NaN to 0. This happens when all the peaks are zero resulting in a SumArea of 0.
  data$Area2SumRatio[!is.finite(data$Area2SumRatio)] <- 0

  tmp = data %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,Area2SumRatio,IsotopeLabelType) %>%
    spread(key = IsotopeLabelType,value = Area2SumRatio)

  colnames(tmp)[grepl(endogenous.label,colnames(tmp))] = "endogenous"
  colnames(tmp)[grepl(standard.label,colnames(tmp))] = "standard"

  data = tmp %>%
    mutate(PairRatioConsistency = abs(endogenous - standard)/standard) %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,PairRatioConsistency) %>%
    right_join(data,by = c("PeptideModifiedSequence", "FragmentIon","PrecursorCharge","ProductCharge","FileName","File"))

  # impute the values that are not finite (NA,NaN (0/0) and Inf (number/0)) to the maximum observed for PairRatioConsistency
  data$PairRatioConsistency[!is.finite(data$PairRatioConsistency)] <- max(data$PairRatioConsistency[is.finite(data$PairRatioConsistency)],na.rm = TRUE)


  # calculate the peak width consistency between the light and heavy for each sample and peptide and tranistion
  # PairFWHMConsistency: abs(endogenous - standard)/standard
  tmp = data %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,TransitionFWHM,IsotopeLabelType) %>%
    spread(key = IsotopeLabelType,value = TransitionFWHM)

  colnames(tmp)[grepl(endogenous.label,colnames(tmp))] = "endogenous"
  colnames(tmp)[grepl(standard.label,colnames(tmp))] = "standard"

  data = tmp %>%
    mutate(PairFWHMConsistency = abs(endogenous - standard)/standard) %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,PairFWHMConsistency) %>%
    right_join(data,by = c("PeptideModifiedSequence", "FragmentIon","PrecursorCharge","ProductCharge","FileName","File"))

  # impute the values that are not finite (NA,NaN (0/0) and Inf (number/0)) to the maximum observed for PairFWHMConsistency
  data$PairFWHMConsistency[!is.finite(data$PairFWHMConsistency)] <- max(data$PairFWHMConsistency[is.finite(data$PairFWHMConsistency)],na.rm = TRUE)


  # calculate the area ratio consistency of each transition with average of area ratios in all samples
  # MeanArea2SumRatio: mean of Area2SumRatio across all samples is calculated for each transition and isotope label.
  # MeanIsotopeRatioConsistency: abs(Area2SumRatio - MeanArea2SumRatio)/MeanArea2SumRatio)
  data = data %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,Area2SumRatio,IsotopeLabelType,aboveThreshold) %>%
    filter(aboveThreshold) %>%
    group_by(PeptideModifiedSequence,PrecursorCharge,File,FragmentIon,ProductCharge,IsotopeLabelType) %>%
    summarise(MeanArea2SumRatio = mean(Area2SumRatio)) %>%
    right_join(data,by = c("PeptideModifiedSequence","FragmentIon","PrecursorCharge","ProductCharge","IsotopeLabelType","File"))

  data = data %>%
    mutate(MeanIsotopeRatioConsistency = abs(Area2SumRatio - MeanArea2SumRatio)/MeanArea2SumRatio)

  # impute the NA MeanIsotopeRatioConsistency values to max of this column. This happens when for a certain fragment ion of a peptide, nothing is detected above a certain threshold (either determined by intensity at lloq or chosen arbitrarily by the user). In these cases we just assume that the peak is too low to be used for determining the quality of the peak. Also impute values that are not finite (NA,NaN (0/0) and Inf (number/0))
  data$MeanIsotopeRatioConsistency[!is.finite(data$MeanIsotopeRatioConsistency)] <- max(data$MeanIsotopeRatioConsistency[is.finite(data$MeanIsotopeRatioConsistency)],na.rm = TRUE)


  # calculate the peak width consistency of each transition with average of peak widths in all samples
  # MeanTransitionFWHM: mean of TransitionFWHM across all samples is calculated for each transition and isotope label.
  # MeanIsotopeFWHMConsistency: abs(TransitionFWHM - MeanTransitionFWHM)/MeanTransitionFWHM)
  data = data %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,TransitionFWHM,IsotopeLabelType,aboveThreshold) %>%
    filter(aboveThreshold) %>%
    group_by(PeptideModifiedSequence,PrecursorCharge,File,FragmentIon,ProductCharge,IsotopeLabelType) %>%
    summarise(MeanTransitionFWHM = mean(TransitionFWHM)) %>%
    right_join(data,by = c("PeptideModifiedSequence","FragmentIon","PrecursorCharge","ProductCharge","IsotopeLabelType","File"))

  data = data %>%
    mutate(MeanIsotopeFWHMConsistency = abs(TransitionFWHM - MeanTransitionFWHM)/MeanTransitionFWHM)


  # impute the NA MeanIsotopeFWHMConsistency values to max of this column. This happens when for a certain fragment ion of a peptide, nothing is detected above a certain threshold (either determined by intensity at lloq or chosen arbitrarily by the user). In these cases we just assume that the peak is to low to be used for determining the quality of the peak. Also impute values that are not finite (NA,NaN (0/0) and Inf (number/0))
  data$MeanIsotopeFWHMConsistency[!is.finite(data$MeanIsotopeFWHMConsistency)] <- max(data$MeanIsotopeFWHMConsistency[is.finite(data$MeanIsotopeFWHMConsistency)],na.rm = TRUE)


  # calculate the peak center (RT) consistency of each transition with average of RT in all samples
  # MeanPeakCenter: mean of PeakCenter across all samples is calculated for each transition and isotope label.
  # MeanIsotopeRTConsistency: abs(PeakCenter - MeanPeakCenter)/MeanPeakCenter

  data = data %>%
    select(PeptideModifiedSequence,PrecursorCharge,File,FileName,FragmentIon,ProductCharge,PeakCenter,IsotopeLabelType,aboveThreshold) %>%
    filter(aboveThreshold) %>%
    group_by(PeptideModifiedSequence,PrecursorCharge,File,FragmentIon,ProductCharge,IsotopeLabelType) %>%
    summarise(MeanPeakCenter = mean(PeakCenter)) %>%
    right_join(data,by = c("PeptideModifiedSequence","FragmentIon","PrecursorCharge","ProductCharge","IsotopeLabelType","File"))

  data = data %>%
    mutate(MeanIsotopeRTConsistency = abs(PeakCenter - MeanPeakCenter)/MeanPeakCenter)

  # impute the NA MeanIsotopeRTConsistency values to max of this column. This happens when for a certain fragment ion of a peptide, nothing is detected above a certain threshold (either determined by intensity at lloq or chosen arbitrarily by the user). In these cases we just assume that the peak is to low to be used for determining the quality of the peak. Also impute values that are not finite (NA,NaN (0/0) and Inf (number/0))
  data$MeanIsotopeRTConsistency[!is.finite(data$MeanIsotopeRTConsistency)] <- max(data$MeanIsotopeRTConsistency[is.finite(data$MeanIsotopeRTConsistency)],na.rm = TRUE)

  # calculate the correlation between the light and heavy Area2SumRatios
  # the sum transitions are first filtered out, and the correlation coefficiant between the Area2SumRatio of light and heavy isotopes are calculated for each peptide and sample
  # also calculate the absolute value of log of Area2SumRatio between endogenous and standard. This helps identify transitions that have interference or different area ratios between the isotope pairs
  tmp = data %>%
    filter(FragmentIon != "sum") %>%
    select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Area2SumRatio) %>%
    spread(key = IsotopeLabelType,value = Area2SumRatio)

  colnames(tmp)[grepl(endogenous.label,colnames(tmp))] = "endogenous"
  colnames(tmp)[grepl(standard.label,colnames(tmp))] = "standard"


  data = suppressWarnings(tmp %>%
    group_by(PeptideModifiedSequence,PrecursorCharge, FileName,File) %>%
    summarise(PeakGroupRatioCorr = cor(endogenous,standard,method = "pearson")) %>%
    right_join(data,by = c("PeptideModifiedSequence", "PrecursorCharge","FileName","File")))

  # impute NA PeakGroupRatioCorr values to 0. The NA values are a result of all zero signals or peptides with only one transition
  data$PeakGroupRatioCorr[is.na(data$PeakGroupRatioCorr)] <- 0

  # calculate the endogenous to standard ratio for each transition

  tmp = data %>%
    select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,PeakArea) %>%
    spread(key = IsotopeLabelType,value = PeakArea)

  colnames(tmp)[grepl(endogenous.label,colnames(tmp))] = "endogenous"
  colnames(tmp)[grepl(standard.label,colnames(tmp))] = "standard"

  data = tmp %>%
    mutate(Endogenous2StandardRatio = endogenous/standard) %>%
    select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,Endogenous2StandardRatio) %>%
    right_join(data,by = c("File","FileName","PeptideModifiedSequence", "PrecursorCharge","FragmentIon","ProductCharge"))



  # for each transition, calculate the CV% of the area to sum ration across all samples
  data = data %>%
    select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Area2SumRatio) %>%
    group_by(File,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType) %>%
    mutate(Area2SumRatioCV = sd(Area2SumRatio,na.rm = TRUE)/mean(Area2SumRatio,na.rm = TRUE)) %>%
    select(-Area2SumRatio) %>%
    right_join(data,by = c("File","FileName","PeptideModifiedSequence", "PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType"))


  # deactivate the cluster
  if (parallel) stopCluster(cl)

  # extract the feature columns
  features = data %>%
    select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,IsotopeLabelType,Area2SumRatioCV,PeakGroupRatioCorr,PairFWHMConsistency,PairRatioConsistency,PeakGroupJaggedness.m,PeakGroupSymmetry.m,PeakGroupSimilarity.m,PeakGroupShift,PeakGroupFWHM.m,PeakGroupFWHM2base.m,PeakGroupModality.m,TransitionJaggedness,TransitionSymmetry,PairSimilarity,PairShift,TransitionFWHM2base,TransitionFWHM,TransitionModality,IsotopeJaggedness,IsotopeSymmetry,IsotopeSimilarity,IsotopeShift,IsotopeFWHM2base,IsotopeFWHM,IsotopeModality,MeanIsotopeRatioConsistency,MeanIsotopeFWHMConsistency,MeanIsotopeRTConsistency,TransitionMaxIntensity,TransitionMaxBoundaryIntensity,TransitionMaxBoundaryIntensityNormalized,TransitionShift)

  # change the IsotopeLabelType to endogenous and standard

  features$IsotopeLabelType <- as.character(features$IsotopeLabelType)
  features$IsotopeLabelType[features$IsotopeLabelType == endogenous.label] = "endogenous"
  features$IsotopeLabelType[features$IsotopeLabelType == standard.label] = "standard"
  features$IsotopeLabelType <- factor(features$IsotopeLabelType,levels = c("endogenous","standard"))

  # convert the feature.data to wide format
  features.to.cast <- colnames(features)[!(colnames(features) %in% c("File","FileName","PeptideModifiedSequence","PrecursorCharge","FragmentIon","ProductCharge","IsotopeLabelType"))]

  features <- dcast(setDT(features),File + FileName + PeptideModifiedSequence + PrecursorCharge + FragmentIon + ProductCharge ~ IsotopeLabelType,value.var = features.to.cast)

  # removing duplicated columns (this happens for the features where the values of standard and endogenous are identical e.g. the peakgroup features)
  features <- features[ ,which(!duplicated(t(features))),with = FALSE]


  # if export.features is true save the features in the csv file specified by feature.path
  if (export.features == TRUE) {


    # template file name
    feature.file <- file.path(feature.path,"features.csv")

    # if the template path does not exist, create it
    if (!dir.exists(feature.path)) {
      dir.create(feature.path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }

    write.table(features, file = feature.file, sep = ",",row.names = FALSE)
  }

  # return output
  return(list(features = features))

}

#' Extract the peak from a chromatogram by restricting the chromatogram to peak boundaries.
#'
#' The function takes a peak group of class peakObj and the peak boundary as input. The input peak group is usually a chromatogram with a time range that is wider than the peak of interest. This function extracts peak of interest from the chromatogram by limiting the time and signal intensities of the peak to a range that is within the provided peak boundaries.
#'
#' @param peak A peak group of class peakObj
#' @param boundary a numeric vector of size two with min start time and max end time of the peak (peak boundaries)
#'
#' @return A peak group of class peakObj with time and intensity vectors that are within the boundaries of the peak
#'
#' @export
#'
#' @examples
#'
#' chrom <- data.CSF$data$ChromGroup[[1]]
#' PlotChromPeak(chrom)
#' peak <- ApplyPeakBoundary(chrom,c(data.CSF$data$MinStartTime[[1]],data.CSF$data$MaxEndTime[[1]]))
#' PlotChromPeak(peak)

ApplyPeakBoundary <- function(peak, boundary,...) {

  # error and warning handling ---------------------------------------

  #   input errors: if the input peak is empty
  error.input.format <- simpleError("ApplyPeakBoundary: input peak should be non-empty and of class peakObj")

  if (is.na(peak)) stop(error.input.format)

  #   input errors: if the boundary is not a numeric vector of length 2
  error.boundary.format <- simpleError("ApplyPeakBoundary: boudary should be a numeric vector of length 2")

  if (length(boundary) != 2) stop(error.boundary.format)

  #   if the peak boundaries are NA
  warning.na.boundary <- simpleWarning("ApplyPeakBoundary: the peak boundaries are NA")

  if (sum(is.na(boundary)) > 0) {
    warnings(warning.na.boundary)

    # return the original peak as input
    return(peak)
  }

  # function body  ---------------------------------------

  # determine the timepoints and the signal within the boundaries

  # if there are only 3 time points within the peak boundary pad it with one additional datapoint. If there are fewer that 3, the peak is removed from qc analysis.
  time.index <- which(peak@time > boundary[1] & peak@time < boundary[2])
  if (length(time.index) == 3) {
    if (time.index[1] > 1) time.index <- c(time.index[1] - 1,time.index)
    else time.index <- c(time.index,tail(time.index,1) + 1)
  }

  time <- peak@time[time.index]

  sig <- peak@sig %>%
    slice(time.index)

  # calculate peak area using the trapozoidal approximation
  area <- sapply(sig,function(x) pracma::trapz(time,x))

  # create a new peak with the time and sig data points within the boundaries
  peak <- tryCatch({
    peakObj(time = time, sig = sig, area = area)
  }, error = function(e) {
    return(NA)
  }
  )

  return(peak)
}
