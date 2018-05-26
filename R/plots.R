## This file contains the customized plotting functions for targeted MS data analysis

#' Plot the chromatographic peaks of class peakObj
#'
#' The function takes a peak group of class peakObj as input and plots its chromatogram. It can generate a static (ggplot) or interactive (plotly) plot.
#'
#' @param peak A peak group of class peakObj
#' @param plottype This parameter can take values of "ggplot" for static and "plotly" for interactive plots
#' @param split.label A logical parameter. If split.label = TRUE, the peaks corresponding to different isotope label types will be split into two panels
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte: default is "heavy"
#' @param transition.list List of transitions to be plotted. If transition.list = NULL, all transitions will be plotted. Example: transition.list = c("b3","y4")
#' @param label.list List of isotope label types to be plotted. If label.list = NULL, all isotope label types will be plotted. Example: label.list = c("light")
#' @param font.size Font size of text in the plot
#'
#' @return A plot object of class ggplot or plotly
#'
#' @export
#'
#' @rawNamespace import(ggplot2, except = margin)
#' @import grid
#' @importFrom plotly ggplotly
#' @importFrom tidyr gather
#'
#' @examples
#'
#' peak <- data.CSF$data$PeakGroup[[1]]
#' PlotChromPeak(peak,transition.list = c("b5","y6"), label.list = NULL,split.label = FALSE)
#' PlotChromPeak(peak,split.label = TRUE)

PlotChromPeak <- function(peak,plottype = c("ggplot","plotly"), split.label = TRUE, endogenous.label = "light", standard.label = "heavy", transition.list = NULL, label.list = NULL, font.size = 20, ...) {

  # error and warning handling ---------------------------------------

  #   input errors: if the input vectors are na or not of equal length
  error.input.format = simpleError("PlotChromPeak: incompatible input")

  if (is.na(peak)) stop(error.input.format)

  #   if the specified plottype is not acceptable
  error.plottype <- simpleError("PlotChromPeak: 'plottype' should be either 'ggplot' or 'plotly'")

  # function body  ---------------------------------------

  #   match inputs with options
  plottype <- tryCatch(match.arg(plottype),
                   error = function(e) stop(error.plottype))

  #   convert peak to dataframe
  d <- data.frame(peak@time,peak@sig) %>%
    gather(transition,intensity,-peak.time)

  d$transition.name <- as.factor(substr(d$transition,1,regexpr("\\.[^\\.]*$",d$transition) - 1))
  d$label <- factor(substr(d$transition,regexpr("\\.[^\\.]*$",d$transition) + 1,100),levels = c(endogenous.label,standard.label))

  # limit the plot to transitions in transition.list
  d.transition <- data.frame()

  if (!is.null(transition.list)) {
    for (t in transition.list) {
      d.transition <- rbind(d.transition,
                            d %>% filter(grepl(t,transition.name)))
    }
    d <- d.transition
  }

  # limit the plot to isotope labels in label.list
  d.label <- data.frame()

  if (!is.null(label.list)) {
    for (l in label.list) {
      d.label <- rbind(d.label,
                            d %>% filter(grepl(l,label)))
    }
    d <- d.label
  }

  #   ggplot
  p <- tryCatch(
    ggplot(d,aes(x = peak.time,y = intensity,colour = transition.name,linetype = label)) +
      geom_line(size = 2) +
      scale_x_continuous(name = "Retention Time (min)") +
      scale_y_continuous(name = "Intensity") +
      scale_linetype_manual(values = c("solid", "twodash")) +
      theme_classic() +
      theme(text = element_text(size = font.size))
  )

  if (split.label) p <- p + facet_grid(label ~ ., scales = "free")

  # return an interactive plot
  if (plottype == "plotly") p <- plotly::ggplotly(p)

  return(p)
}


#' Graph violin plots to visualize distribution of QC features for a dataset
#'
#' The function takes a dataframe and list of QC features of interest as input and graphs violin plots for distribution of these features. The results are returned as a list of ggplot objects. Each element in the list corresponds to a mass spectrometry run in the input dataframe.
#'
#' @param data A dataframe that contains columns for features of interest and FileName where FileName refers to names of individual mass spectrometry runs included in the input. This function is mainly used in conjunction with the MakeDataSet function. data.set$data.merged is used as the input where data.set is the output of MakeDataSet.
#' @param runs List of indivodual mass spectrometry runs to be plotted. Default is "all".
#' @param features List of QC features of interest to be plotted. Default is "all".
#' @param labels Name of the column that holds the annotated labels for each peak. When labels are provided, the violin plots will be overlaid with dot plots, where the dots are colored based on their corresponding label.
#' @param response.var If the input dataframe contains columns corresponding to response variables (such as Status or Status.prediction), it should be indicated here. Response and description columns as well as identifier columns will be removed from the data before creating the violin plots.
#' @param description.columns If the input dataframe contains columns corresponding to description variables (Such as Notes), it should be indicated here. Response and description columns as well as identifier columns will be removed from the data before creating the violin plots.
#' @param font.size Font size of text in the plot
#'
#' @return  List of plots of class ggplot for each run in the input data
#'
#' @export
#'
#' @rawNamespace import(ggplot2, except = margin)
#' @import grid
#' @importFrom tidyr gather_
#'
#' @examples
#'
#' feature.set <- c("Area2SumRatioCV_standard","TransitionJaggedness_standard","TransitionSymmetry_standard","TransitionFWHM2base_standard","TransitionFWHM_standard","TransitionModality_standard","TransitionShift_standard")
#' violin.plots = ViolinPlotQCSummary(data.set.CSF$data.merged,runs = "all",features = feature.set,font.size = 15, labels = "Status")
#' violin.plots[[1]]

ViolinPlotQCSummary <- function(data, runs = "all", features = "all", labels = NULL, response.var = c("Status"), description.columns = c("Notes"), font.size = 50, ...){

  # error and warning handling ---------------------------------------

  if (length(runs) == 1 && tolower(runs) == "all") runs <- unique(data$FileName)

  # identifier columns should be removed from the list of features
  identifier.columns = c("File","FileName","PeptideModifiedSequence","FragmentIon",
                         "IsotopeLabelType","PrecursorCharge","ProductCharge")

  # response.var and description columns should be removed from the list of features
  if (length(features) == 1 && tolower(features) == "all")
    features <- colnames(data)[!colnames(data) %in% identifier.columns & !colnames(data) %in% description.columns & !colnames(data) %in% response.var]

  # initiating the list of output plots
  count <- 0
  v.plots <- list()

  # create a violin plot for each run
  for (run in runs) {
    count <- count + 1
    tmp <- data %>% filter(FileName == run) %>%
      gather_(key_col = "feature",value_col = "score",gather_cols = features)

    v.plots[[count]] <- ggplot(tmp,aes(feature,score)) +
      geom_violin(fill = "#F2F2F2",scale = "width") +
      coord_flip() +
      ggtitle(run) +
      theme_classic() +
      theme(text = element_text(size = font.size))

    # if a vector of labels is provided, overlay the violin plot with a jittered dot plot
    if (!is.null(labels)) {
      v.plots[[count]] <- v.plots[[count]] + geom_jitter(aes_string(color = labels), position = position_jitter(0.25),size = 0.5)
    }
  }

  # create a violing plot for all the runs
  count <- count + 1
  tmp <- data %>%
      gather_(key_col = "feature",value_col = "score",gather_cols = features)

  v.plots[[count]] <- ggplot(tmp,aes(feature,score)) +
    geom_violin(fill = "#F2F2F2",scale = "width") +
    coord_flip() +
    ggtitle("All Runs") +
    theme_classic() +
    theme(text = element_text(size = font.size))

  # if a vector of labels is provided, overlay the violin plot with a jittered dot plot
  if (!is.null(labels)) {
    v.plots[[count]] <- v.plots[[count]] + geom_jitter(aes_string(color = labels), position = position_jitter(0.25),size = 0.5)
  }

  # return output
  return(v.plots)
}


#' Draws violin plots to  show distribution of QC features for a dataset
#'
#' A longer description of how the function works.
#'
#' @param response.data A dataframe with identifier columns (File, FileName, PeptideModifiedSequence, FragmentIon, PrecursorCharge ,ProductCharge) and the response variable column (Status.prediction by default). Usually, the input dataframe is the output of the ApplyQCModel function.
#' @param response.var The name of the column in the input data that contains the response variable to be used for creating the dataset QC report. This variable is usually the predicted status for each peak generated by the ApplyQCModel function.
#' @param report.path The path to directory where the pdf file of the QC report should be saved.
#' @param min.ok.transitions.threshold Minimum number of transitions that should pass QC for a peptide to pass QC in each sample.
#' @param plot.prob A logical value to indicate whether the class probabilities should be used in heatmaps. If plot.prob = FALSE,  of the class probabilities, the assigned binary class will be used.
#' @param flag.class.prob.var The name of the response.data column that holds class probabilities for the "flag" class. This parameter is used only if plot.prob = TRUE.
#'
#'
#' @return Saves a .pdf file of the QC report in report.path
#'
#' @export
#'
#' @rawNamespace import(ggplot2, except = margin)
#' @import grid
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom ComplexHeatmap draw
#' @importFrom circlize colorRamp2
#'
#'
#' @examples
#'
#' response.data <- ApplyQCModel(data.set.CSF$feature.data,
#'                               model.rrf.CSF,
#'                               response.var = c("Status"),
#'                               description.columns = c("Notes"),
#'                               flag.prob.threshold = 0.5,
#'                               type = "prob")
#' PlotQCReport(response.data,report.path = "",response.var = c("Status.prediction"),plot.prob = TRUE)


PlotQCReport <- function(response.data,response.var = c("Status.prediction"), report.path = "", min.ok.transitions.threshold = 3, plot.prob = FALSE, flag.class.prob.var = c("flag.prob.prediction"), ...){

  # function body  ---------------------------------------

  # determine OS
  if (grepl("unix",.Platform$OS.type)) {
    mark =  "/"
  } else {
    mark = "\\"
  }

  # identifier columns are used to clean up the input response data
  identifier.columns = c("File","FileName","PeptideModifiedSequence","FragmentIon",
                         "PrecursorCharge","ProductCharge")

  if (!plot.prob) {
    # only identifier columns and response columns are kept
    response.data <- response.data[,c(identifier.columns,response.var)]

    # replace the column name of the response.var
    colnames.vec <- colnames(response.data)
    colnames.vec[colnames.vec == response.var] <- "Status"
    colnames(response.data) <- colnames.vec

  } else {
    # only identifier columns and response columns are kept
    response.data <- response.data[,c(identifier.columns,response.var,flag.class.prob.var)]

    # replace the column name of the response.var
    colnames.vec <- colnames(response.data)
    colnames.vec[colnames.vec == response.var] <- "Status"
    colnames.vec[colnames.vec == flag.class.prob.var] <- "Flag.Class.Prob"
    colnames(response.data) <- colnames.vec

  }

  # if any of the "sum" transitions are imported, remove them from the QC report.
  response.data <- response.data %>% filter(FragmentIon != "sum")

  # for each skyline doc in the data:
  for (File_ in unique(response.data$File)) {

    # the report will be saved as a pdf file in the report path directory
    file.name <- paste0(File_,"_TargetedMSQC_report",format(Sys.time(), "_%Y%m%d_%H%M"),".pdf")
    report.file <- file.path(report.path,file.name)
    pdf(report.file)

    # filter rows for the skyline file
    response.data.file <- response.data %>% filter(File == File_)


    # calculate the number of "ok" transitions for each peptide and precursor charge combination
    response.data.file.plot <- response.data.file %>%
      group_by(FileName,PeptideModifiedSequence,PrecursorCharge) %>%
      summarize(Ok.Transition.No = sum(Status == "ok")) %>%
      ungroup()

    ## generating the global QC bargraphs
    response.data.file.bargraph.plot <- response.data.file.plot

    # create a single identifier for peptide and charge combinations
    response.data.file.bargraph.plot$Peptide.Charge <- paste(response.data.file.bargraph.plot$PeptideModifiedSequence,response.data.file.bargraph.plot$PrecursorCharge,sep = ".")

    # for peptides whose number of ok transitions does not meet the min.ok.transitions.threshold criteria, the QC.Status is changed to Fail
    response.data.file.bargraph.plot$QC.Status <- "Pass"
    response.data.file.bargraph.plot$QC.Status[response.data.file.bargraph.plot$Ok.Transition.No < min.ok.transitions.threshold] <- "Fail"

    response.data.file.bargraph.plot$QC.Status <- factor(response.data.file.bargraph.plot$QC.Status, levels = c("Fail","Pass"))

    # shorten the filename by removing the prefix and suffix
    prefix <- Biobase::lcPrefix(response.data.file.bargraph.plot$FileName,ignore.case = TRUE)
    suffix <- Biobase::lcSuffix(response.data.file.bargraph.plot$FileName,ignore.case = TRUE)

    response.data.file.bargraph.plot$FileName <- gsub(prefix,"",response.data.file.bargraph.plot$FileName,ignore.case = TRUE)
    response.data.file.bargraph.plot$FileName <- gsub(suffix,"",response.data.file.bargraph.plot$FileName,ignore.case = TRUE)


    # plot the number of peptides that passed QC for each sample
    sample.summary.plot <- ggplot(response.data.file.bargraph.plot, aes(FileName, fill = QC.Status)) +
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5)) +
      xlab("Sample") +
      ylab("No. of Peptides")

    print(sample.summary.plot)

    # plot the number of samples that passed QC for each peptides
    peptide.summary.plot <- ggplot(response.data.file.bargraph.plot, aes(Peptide.Charge, fill = QC.Status)) +
      geom_bar() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)) +
      xlab("Peptide.Charge") +
      ylab("No. of Samples")

    print(peptide.summary.plot)


    ## generating the global heatmap

    response.data.file.heatmap.plot <- response.data.file.plot

    # create a single identifier for peptide and charge combinations
    response.data.file.heatmap.plot$Peptide.Charge <- paste(response.data.file.heatmap.plot$PeptideModifiedSequence,response.data.file.heatmap.plot$PrecursorCharge,sep = ".")

    # remove extra columns
    response.data.file.heatmap.plot <- response.data.file.heatmap.plot %>% select(-PeptideModifiedSequence) %>% select(-PrecursorCharge)

    # change to wide format for plotting with complexheatmap
    response.data.file.heatmap.plot <- response.data.file.heatmap.plot %>% spread(key = Peptide.Charge,value = Ok.Transition.No)

    # assign the file name to rows before changing the df to mat
    rownames(response.data.file.heatmap.plot) <- response.data.file.heatmap.plot$FileName

    # format to matrix
    response.data.file.heatmap.plot <- as.matrix(response.data.file.heatmap.plot %>% select(-FileName))

    # transpose the matrix
    response.data.file.heatmap.plot <- t(response.data.file.heatmap.plot)

    # shorten the filename by removing the prefix and suffix
    prefix <- Biobase::lcPrefix(colnames(response.data.file.heatmap.plot),ignore.case = TRUE)
    suffix <- Biobase::lcSuffix(colnames(response.data.file.heatmap.plot),ignore.case = TRUE)

    colnames(response.data.file.heatmap.plot) <- gsub(prefix,"",colnames(response.data.file.heatmap.plot),ignore.case = TRUE)
    colnames(response.data.file.heatmap.plot) <- gsub(suffix,"",colnames(response.data.file.heatmap.plot),ignore.case = TRUE)

    # missing values correspond to filename and peptides with no good transition. should be imputed to zero.
    response.data.file.heatmap.plot[is.na(response.data.file.heatmap.plot)] <- 0


    # calculate the  number of transition monitored for each peptides
    tmp <- response.data.file %>%
      group_by(FileName,PeptideModifiedSequence,PrecursorCharge) %>%
      summarize(Total.Transition.No = n())

    # calcualte the maximum number of transitions monitored among all peptides
    max.total.transitions <- max(tmp$Total.Transition.No)


    # heatmap plot of the skyline QC summary
    h.file <- ComplexHeatmap::Heatmap(response.data.file.heatmap.plot,
                 col = colorRamp2(c(0,min.ok.transitions.threshold,max(c(max.total.transitions,min.ok.transitions.threshold + 1))), c("#EC7063","white","#48C9B0")),
                 row_title = "",
                 column_title = File_,
                 cluster_rows = F,
                 cluster_columns = F,
                 show_column_names = T,
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp = gpar(fontsize = 5),
                 column_names_max_height = unit(10, "cm"),
                 name = "No. of High Quality Transitions",
                 heatmap_legend_param = list(legend_direction = "horizontal",title_position = "lefttop"))


    draw(h.file,heatmap_legend_side = "top")

    # Generating QC plots for each peptide
    peptide.charge.vec <- unique(response.data.file[,c("PeptideModifiedSequence","PrecursorCharge")])

    # for each peptide and charge combination:
    for (counter in 1:nrow(peptide.charge.vec)) {

      Peptide.Charge_ <- peptide.charge.vec[counter,]

      # first filter the rows corresponding to the file, peptide and charge combination
      response.data.file.peptide <- response.data.file %>%
        filter(PeptideModifiedSequence == Peptide.Charge_$PeptideModifiedSequence & PrecursorCharge == Peptide.Charge_$PrecursorCharge)

      if (!plot.prob) {
        # determine which transitions passed QC in each run
        response.data.file.peptide.plot <- response.data.file.peptide %>%
          group_by(FileName,FragmentIon,ProductCharge) %>%
          mutate(Ok.Transition.No = sum(Status == "ok"))

        # remove the extra columns
        response.data.file.peptide.plot <- unique(response.data.file.peptide.plot %>%
                                                    select(FileName,FragmentIon,ProductCharge,Ok.Transition.No)) %>%
          ungroup()

        # create an identifier for each peptide and charge
        response.data.file.peptide.plot$Fragment.Charge <- paste(response.data.file.peptide.plot$FragmentIon,response.data.file.peptide.plot$ProductCharge,sep = ".")

        # remove extra columns
        response.data.file.peptide.plot <- response.data.file.peptide.plot %>% select(-FragmentIon) %>% select(-ProductCharge)

        # change to wide format for the plotting the heatmap with complexHeatmap
        response.data.file.peptide.plot <- response.data.file.peptide.plot %>% spread(key = Fragment.Charge,value = Ok.Transition.No)

        # assign the filename (run) to rownames before converting to matrix
        rownames(response.data.file.peptide.plot) <- response.data.file.peptide.plot$FileName

        # convert to matrix for use with the complexHeatmap package
        response.data.file.peptide.plot <- as.matrix(response.data.file.peptide.plot %>% select(-FileName))

        # transpose the matrix
        response.data.file.peptide.plot <- t(response.data.file.peptide.plot)

        # shorten the filename by removing the prefix and suffix
        prefix <- Biobase::lcPrefix(colnames(response.data.file.peptide.plot),ignore.case = TRUE)
        suffix <- Biobase::lcSuffix(colnames(response.data.file.peptide.plot),ignore.case = TRUE)

        colnames(response.data.file.peptide.plot) <- gsub(prefix,"",colnames(response.data.file.peptide.plot),ignore.case = TRUE)
        colnames(response.data.file.peptide.plot) <- gsub(suffix,"",colnames(response.data.file.peptide.plot),ignore.case = TRUE)

        # missing values correspond to filename and peptides with no good transition. should be imputed to zero.
        response.data.file.peptide.plot[is.na(response.data.file.peptide.plot)] <- 0

        # just change the 0 and 1 labels to "Fail" and "Pass"
        response.data.file.peptide.plot[response.data.file.peptide.plot == 1] = "Pass"
        response.data.file.peptide.plot[response.data.file.peptide.plot == 0] = "Fail"

        # set columns for "Pass" and "Fail" in the heatmap
        disc.color <- as.character(c("#48C9B0","#EC7063"))
        names(disc.color) <- c("Pass","Fail")

        # create the QC summary heatmap for each peptide
        h.file.peptide <- ComplexHeatmap::Heatmap(response.data.file.peptide.plot,
                                  col = disc.color,
                                  row_title = "",
                                  column_title = paste(Peptide.Charge_$PeptideModifiedSequence,Peptide.Charge_$PrecursorCharge,sep = "."),
                                  column_title_gp = gpar(fontsize = 8),
                                  cluster_rows = F,
                                  cluster_columns = F,
                                  show_column_names = T,
                                  row_names_gp = gpar(fontsize = 8),
                                  column_names_gp = gpar(fontsize = 3),
                                  column_names_max_height = unit(10, "cm"),
                                  name = "QC Status",
                                  heatmap_legend_param = list(legend_direction = "horizontal",title_position = "lefttop"),
                                  rect_gp = gpar(col = "white"))


      } else {
        # remove the extra columns
        response.data.file.peptide.plot <- unique(response.data.file.peptide %>%
                                                    select(FileName,FragmentIon,ProductCharge,Flag.Class.Prob))

        # create an identifier for each peptide and charge
        response.data.file.peptide.plot$Fragment.Charge <- paste(response.data.file.peptide.plot$FragmentIon,response.data.file.peptide.plot$ProductCharge,sep = ".")

        # remove extra columns
        response.data.file.peptide.plot <- response.data.file.peptide.plot %>% select(-FragmentIon) %>% select(-ProductCharge)

        # change to wide format for the plotting the heatmap with complexHeatmap
        response.data.file.peptide.plot <- response.data.file.peptide.plot %>% spread(key = Fragment.Charge,value = Flag.Class.Prob)

        # assign the filename (run) to rownames before converting to matrix
        rownames(response.data.file.peptide.plot) <- response.data.file.peptide.plot$FileName

        # convert to matrix for use with the complexHeatmap package
        response.data.file.peptide.plot <- as.matrix(response.data.file.peptide.plot %>% select(-FileName))

        # transpose the matrix
        response.data.file.peptide.plot <- t(response.data.file.peptide.plot)

        # shorten the filename by removing the prefix and suffix
        prefix <- Biobase::lcPrefix(colnames(response.data.file.peptide.plot),ignore.case = TRUE)
        suffix <- Biobase::lcSuffix(colnames(response.data.file.peptide.plot),ignore.case = TRUE)

        colnames(response.data.file.peptide.plot) <- gsub(prefix,"",colnames(response.data.file.peptide.plot),ignore.case = TRUE)
        colnames(response.data.file.peptide.plot) <- gsub(suffix,"",colnames(response.data.file.peptide.plot),ignore.case = TRUE)

        # create the QC summary heatmap for each peptide
        h.file.peptide <- ComplexHeatmap::Heatmap(response.data.file.peptide.plot,
                                  col = colorRamp2(c(0, 0.5, 1), c("#48C9B0", "white", "#EC7063")),
                                  row_title = "",
                                  column_title = paste(Peptide.Charge_$PeptideModifiedSequence,Peptide.Charge_$PrecursorCharge,sep = "."),
                                  column_title_gp = gpar(fontsize = 8),
                                  cluster_rows = F,
                                  cluster_columns = F,
                                  show_column_names = T,
                                  row_names_gp = gpar(fontsize = 8),
                                  column_names_gp = gpar(fontsize = 3),
                                  column_names_max_height = unit(10, "cm"),
                                  name = "Flag Class Probability",
                                  heatmap_legend_param = list(legend_direction = "horizontal",title_position = "lefttop"),
                                  rect_gp = gpar(col = "white"))


      }

      # create a unique identifier for fragment and charge combinations
      response.data.file.peptide$Fragment.Charge <- paste(response.data.file.peptide$FragmentIon,response.data.file.peptide$ProductCharge,sep = ".")

      # unique transitions for each peptide
      transition.charge.vec <- unique(response.data.file.peptide$Fragment.Charge)

      # Transition.Profile holds the data that is used to find the best transition combination for quantitation
      Transition.Profile <- data.frame(Transition.Pattern = as.character(),Transition.Pattern.Length = as.numeric(),FileName = as.character(),QC.Fail = as.logical())

      # calculate the number of samples that pass QC for each transition combination
      for (transition.comb.n in 1:length(transition.charge.vec)) {

        # transition.comb.n determines the number of transitions in the transition pattern
        transition.comb <- combn(transition.charge.vec,transition.comb.n)


        for (iter in 1:ncol(transition.comb)) {

          # filter rows that correspond to our specific transition pattern
          tmp <- response.data.file.peptide[response.data.file.peptide$Fragment.Charge %in% transition.comb[,iter],]

          # count the number of samples that failed QC if they were quantified based on the specified transition pattern
          tmp <- tmp %>%
            group_by(FileName) %>%
            summarise(QC.Fail = sum(Status == "flag") > 0)

          # create an identifier for the transition pattern
          tmp$Transition.Pattern <- paste(transition.comb[,iter],collapse = "/")

          # include the length of the transition pattern in the final results
          tmp$Transition.Pattern.Length = transition.comb.n

          # merge the results of this transition pattern with the previous result matrix
          Transition.Profile <- rbind(Transition.Profile,tmp)

        }
      }

      # clean up the QC.Status columns for the result matrix
      Transition.Profile$QC.Status <- "Pass"
      Transition.Profile$QC.Status[Transition.Profile$QC.Fail] <- "Fail"
      Transition.Profile$QC.Status <- factor(Transition.Profile$QC.Status,levels = c("Fail","Pass"))

      # determine the order by which the transition patterns should be plotted (order by decrease in the number of successfullly QC'd samples and the length of the transition pattern)
      tmp <- Transition.Profile %>%
        group_by(Transition.Pattern,Transition.Pattern.Length) %>%
        summarise(n = sum(QC.Status == "Pass")) %>%
        arrange(desc(n),desc(Transition.Pattern.Length))

      # apply the order to the result matrix
      Transition.Profile$Transition.Pattern <- factor(Transition.Profile$Transition.Pattern,levels = tmp$Transition.Pattern)

      tr.profile.bar <- ggplot(Transition.Profile) +
        geom_bar(aes(fill = QC.Status,x = Transition.Pattern)) +
        coord_flip() +
        theme_classic() +
        ylab("No. of Samples") +
        xlab("Transitions for Quantitation") +
        guides(fill = guide_legend(title = "QC Status")) +
        theme(legend.position = "bottom") +
        scale_fill_manual("legend", values = c("Pass" = "#48C9B0", "Fail" = "#EC7063"))

      # move to a new page in the pdf report
      grid.newpage()

      # plot the figures for each peptide side by side
      # heatmap
      pushViewport(viewport(x = 0, width = 0.5, just = "left"))
      draw(h.file.peptide,heatmap_legend_side = "top",newpage = FALSE)
      popViewport()

      # bar graph
      print(tr.profile.bar,vp = viewport(x = 0.5, width = 0.5,just = "left"))

    }
    dev.off()
  }
}
