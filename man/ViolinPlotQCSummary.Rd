% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{ViolinPlotQCSummary}
\alias{ViolinPlotQCSummary}
\title{Graph violin plots to visualize distribution of QC features for a dataset}
\usage{
ViolinPlotQCSummary(data, runs = "all", features = "all", labels = NULL,
  response.var = c("Status"), description.columns = c("Notes"),
  font.size = 50, ...)
}
\arguments{
\item{data}{A dataframe that contains columns for features of interest and FileName where FileName refers to names of individual mass spectrometry runs included in the input. This function is mainly used in conjunction with the MakeDataSet function. data.set$data.merged is used as the input where data.set is the output of MakeDataSet.}

\item{runs}{List of indivodual mass spectrometry runs to be plotted. Default is "all".}

\item{features}{List of QC features of interest to be plotted. Default is "all".}

\item{labels}{Name of the column that holds the annotated labels for each peak. When labels are provided, the violin plots will be overlaid with dot plots, where the dots are colored based on their corresponding label.}

\item{response.var}{If the input dataframe contains columns corresponding to response variables (such as Status or Status.prediction), it should be indicated here. Response and description columns as well as identifier columns will be removed from the data before creating the violin plots.}

\item{description.columns}{If the input dataframe contains columns corresponding to description variables (Such as Notes), it should be indicated here. Response and description columns as well as identifier columns will be removed from the data before creating the violin plots.}

\item{font.size}{Font size of text in the plot}
}
\value{
List of plots of class ggplot for each run in the input data
}
\description{
The function takes a dataframe and list of QC features of interest as input and graphs violin plots for distribution of these features. The results are returned as a list of ggplot objects. Each element in the list corresponds to a mass spectrometry run in the input dataframe.
}
\examples{

feature.set <- c("Area2SumRatioCV_standard","TransitionJaggedness_standard","TransitionSymmetry_standard","TransitionFWHM2base_standard","TransitionFWHM_standard","TransitionModality_standard","TransitionShift_standard")
violin.plots = ViolinPlotQCSummary(data.set.CSF$data.merged,runs = "all",features = feature.set,font.size = 15, labels = "Status")
violin.plots[[1]]
}
