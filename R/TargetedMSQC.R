# TargetedMSQC documentation
#
# @author     Shadi Eshghi
# @copyright  Copyright 2017 OMNI-BD, Genentech. All rights reserved
# @license    Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
# @version    0.99.0
# @date       2018.05.11


#' TargetedMSQC: a tool for quality assessment and interference detection in targeted mass spectrometry data using machine learning
#'
#' TargetedMSQC provides a semi-automated workflow for quality control (QC) of chromatographic peaks in targeted proteomics experiments, with the aim of improving the efficiency and reducing the subjectivity of data QC. The package offers a toolkit to build and apply statistical models for predicting peak qualities in proteomics datasets using supervised learning methods. The package contains functions to calculate an ensemble of >30 well-established and newly introduced peak quality metrics, such as jaggedness, FWHM, modality, shift, coefficient of variation, consistency of transition peak area ratios, etc., to quantify the quality of chromatographic peaks. These quality control metrics calculated in a training dataset of peaks with pre-annotated quality status labels are used as the feature set in supervised learning algorithms to flag peaks with poor chromatography or interference in other targeted proteomics experiments.
#'
#' @author Shadi Eshghi, \email{toghiess@@gene.com}
#'
#' @seealso \url{https://github.com/shadieshghi/TargetedMSQC}
#'
#' @name        TargetedMSQC
#' @docType     package
#' @keywords    targeted proteomics, quantitative, mass spectrometry, bioinformatics, automated analysis, interference detection, quality control, machine learning


NULL
