## This file contains the functions used to build and evaluate the training model.


#' Create a training set template that can be shared with an analyst for manual annotation.
#'
#' The function takes the directory that contains chromatogram of Skyline documents and creates a template for manual annotation of the peaks. The template is saved as a .csv file in the directory specified by template.path. This templae should be shared with an analyst for manual annotation and the "Status" column should be filled with "flag" or "ok" for low and high quality peaks, respectively. Any comments about the manual annotation can be entered in the "Notes" column.
#'
#' @param chromatogram.path Path to the directory containing the .tsv files of the peak chromatograms. For each Skyline document, this file is exported from Skyline through File > Export > Chromatograms. Here, check runs of interest and include Precursors, Products, Base Peaks and TICs. Each chromatogram .tsv file corresponds to a single Skyline document, which may contain any number of runs. Multiple chromatogram files, corresponding to multiple Skyline documents can be copied into the chromatogram.path directory. For each chromatogram file in this folder, there should be a peak boundary file with an identical name in peak.boundary.path directory.
#' @param template.path Path to the directory, where the template file will be saved.
#' @param training.filename.list List of the runs that are going to be used for training. if set to defaul of "all", all the runs in the chromatogram file will be used for training.
#' @param endogenous.label Label of the endogenous analyte: default is "light"
#' @param standard.label Label of the spiked-in isotopically labeled analyte: default is "heavy"
#' @param iRT.list List of iRT standards used in the experiment. These peptides will be removed from the training set.
#'
#' @return Saves the template .csv file in template.path
#'
#' @export
#'
#' @rawNamespace import(dplyr, except = combine)
#'
#' @examples
#'
#' extdata.path <- system.file("extdata",package = "TargetedMSQC")
#' project.folder.name <- "CSF_Panel"
#' project.path <- file.path(extdata.path,project.folder.name)
#' chromatogram.path <- file.path(project.path,"Chromatograms")
#' template.path <- file.path(project.path,"Templates")
#' MakeTemplate(chromatogram.path = chromatogram.path,
#'               template.path = template.path,
#'               endogenous.label = "light",standard.label = "heavy")

MakeTemplate <- function(chromatogram.path,template.path, training.filename.list = "all", endogenous.label = "light", standard.label = "heavy" , iRT.list = c("LGGNETQVR","AGGSSEPVTGLADK","VEATFGVDESANK","YILAGVESNK","TPVISGGPYYER","TPVITGAPYYER","GDLDAASYYAPVR","DAVTPADFSEWSK","TGFIIDPGGVIR","GTFIIDPAAIVR","FLLQFGAQGSPLFK","LGGNEQVTR","GAGSSEPVTGLDAK","VEATFGVDESNAK","YILAGVENSK","TPVISGGPYEYR","TPVITGAPYEYR","DGLDAASYYAPVR","ADVTPADFSEWSK","GTFIIDPGGVIR","GTFIIDPAAVIR","LFLQFGAQGSPFLK"), ...) {


  # read all the tsv files in the chromatogram directory
  chromatogram.files <- sort(list.files(chromatogram.path,pattern = "\\.tsv"))
  chromatogram.files <- unlist(lapply(chromatogram.files,function(x) substr(x,1,regexpr("\\.",x)[[1]] - 1)))

  tmp <- lapply(chromatogram.files, function(chromatogram.file) {

    chromatogram.file.path <- file.path(chromatogram.path,paste0(chromatogram.file,".tsv"))

    data <- unique(read.table(chromatogram.file.path,sep = "\t", header = TRUE))

    # assign the skyline file name as a column
    data$File <- chromatogram.file

    # remove the rows with an isotope label that is not included in the pair
    data$IsotopeLabelType <-  tolower(data$IsotopeLabelType)
    data$IsotopeLabelType <- factor(data$IsotopeLabelType,levels = c(endogenous.label,standard.label))

    data <- data %>% filter(!is.na(IsotopeLabelType))

    # remove the rows that correspond to iRTs

    data <- data %>% filter(!PeptideModifiedSequence %in% iRT.list)

    # QC.data holds File,FileName, Peptide, Transition and Isotopelabel Info

    QC.data <- unique(data %>%
                         select(File,FileName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge) %>%
                         arrange(File,PeptideModifiedSequence,PrecursorCharge,FileName,FragmentIon,ProductCharge))

    # add a column for status
    QC.data$Status = NA

    # add a column for additional notes
    QC.data$Notes = NA

    # template file name
    template.file <- file.path(template.path,paste0(chromatogram.file,"_training_template.csv"))

    # if the template path does not exist, create it
    if (!dir.exists(template.path)) {
      dir.create(template.path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }

    # filter to keep the rows in the training.filename.list

    if (!identical(training.filename.list,"all")) {
      QC.data <- QC.data %>% filter(FileName %in% training.filename.list)
    }


    # writing the template
    write.table(QC.data, file = template.file, sep = ",",row.names = FALSE)
  }
  )
}

#' Create the dataset used as input to the TrainQCModel function.
#'
#' The function takes the feature dataframe (containing the calculated ensemble of QC features) and the annotated training dataframe and merges them to create the input of the TrainQCModel function.
#'
#' @param feature.data A dataframe that contains peak identifiers (File,FileName,PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge) as well as  QC metrics calcualted for each transition pair. This dataframe is the output of ExtractFeatures function (output$features).
#' @param feature.path Alternative to providing the feature.data, a path to the directory that contains the features.csv file exported by ExtractFeatures function can be provided. This file can be generated by by setting export.features = TRUE and specifying the feature.path in ExtractFeatures function. If both feature.data and feature.path are provided, feature.path is ignored.
#' @param training.path Path to the directory containing the annotated template file. The template file is generated using the MakeTemplate function should be manually annotated by an expert analyst and saved as a separate .csv file. The path to this annotated file should be provided. Please note that the directory should contain only annotated .csv files that are meant to be in the study.
#'
#'
#' @return A list with the following objects:
#'               data.merged: A dataframe that is product of merging and cleaning up feature.data and training.data.
#'               feature.data: The input feature.data
#'               training.data: The annotated input training.data
#'               data.training.feature: A dataframe that is product of merging and cleaning up feature.data and training.data. This dataframe contain only peaks that have been manually annotated and are included in training.data This data can be used by the TrainQCModel for training a predictive peak QC model.
#'
#' @export
#'
#' @importFrom data.table setDT
#'
#' @examples
#'
#' extdata.path <- system.file("extdata",package = "TargetedMSQC")
#' project.folder.name <- "CSF_Panel"
#' project.path <- file.path(extdata.path,project.folder.name)
#' training.path <- file.path(project.path,"Training")
#' data.set <- MakeDataSet(feature.data = data.features.CSF$features,training.path = training.path)
#' feature.path <- file.path(project.path,"Features")
#' data.set <- MakeDataSet(feature.path = feature.path,training.path = training.path)

MakeDataSet <- function(feature.data = NULL, feature.path = NULL, training.path = NULL, ...) {


  # error and warning handling ---------------------------------------
  # In the future, MakeTemplate will be populated by interactive selections made by user on Skyline.
  # input errors: if feature.data and feature.file are empty

  error.input = simpleError("MakeDataSet: Please provide the feature data frame or path to the feature spreadsheet file")
  warning.input = simpleWarning("MakeDataSet: feature.data is provided. Therefore, feature.path is ignored.")


  if (is.null(feature.data) && is.null(feature.path)) stop(error.input)
  if (!is.null(feature.data) && !is.null(feature.path)) warnings(warning.input)


  # reading and cleaning up the feature data set (this has been generated by the ExtractFeatures function)

  if (is.null(feature.data) && !is.na(feature.path)) {

    feature.files <- sort(list.files(feature.path,pattern = "\\.csv"))
    feature.files <- unlist(lapply(feature.files,function(x) substr(x,1,regexpr("\\.",x)[[1]] - 1)))

    feature.data <- setDT(rbindlist(lapply(feature.files,
                                            function(feature.file) {
                                              feature.file.path <- file.path(feature.path,paste0(feature.file,".csv"))
                                              tmp <- unique(read.csv(feature.file.path,sep = ",", header = TRUE))
                                              return(tmp)
                                            }
    )
    )
    )


    numeric.cols = colnames(feature.data)[!(colnames(feature.data) %in% c("File","FileName",
                                                                          "PeptideModifiedSequence","FragmentIon","IsotopeLabelType","PrecursorCharge","ProductCharge"))]


    feature.data <- as.data.frame(feature.data)

    # when reading the excel file using read.xlsx2, the numeric values are read as levels. so they need to be converted to numeric values.
    for (col in numeric.cols) {
      feature.data[,col] = as.numeric(as.character(feature.data[,col]))
    }

  }



  # removing the rows of FragmentIon = "sum" from the feature.data

  feature.data = feature.data %>% filter(FragmentIon != "sum")
  feature.data$FragmentIon = factor(feature.data$FragmentIon)
  feature.data$File = factor(feature.data$File,levels = sort(as.character(unique(feature.data$File))))


  # reading and cleaning up the training data set (this was the template file generated by MakeTemplate that has been populated by the analyst)

  if (!is.null(training.path))  {

    training.files <- sort(list.files(training.path,pattern = "\\.csv"))
    training.files <- unlist(lapply(training.files,function(x) substr(x,1,regexpr("\\.",x)[[1]] - 1)))

    training.data <- setDT(rbindlist(lapply(training.files,
                                           function(training.file) {
                                             training.file.path <- file.path(training.path,paste0(training.file,".csv"))
                                             tmp <- unique(read.csv(training.file.path,sep = ",", header = TRUE))
                                             return(tmp)
                                           }
    )
    )
    )

    training.data$Status = tolower(as.character(training.data$Status))
    training.data$Notes = tolower(as.character(training.data$Notes))

    training.data = training.data %>%
      filter(Status != "#N/A") %>%
      filter(!is.na(Status))


    training.data[training.data$Notes == "#N/A" | is.na(training.data$Notes),c("Notes")] = "ok"

    training.data$Status = as.factor(training.data$Status)
    training.data$Notes = as.factor(training.data$Notes)

    training.data$File = factor(training.data$File,levels = sort(as.character(unique(training.data$File))))
  }


  # to synchronize the peptides in feature.data and training.data and avoid casting them from factors to characters, let's remove the training set peptides missing from features data frame

  training.data.tmp <- training.data %>%
      filter(PeptideModifiedSequence %in% feature.data$PeptideModifiedSequence)

  training.data.tmp$PeptideModifiedSequence <- droplevels(training.data.tmp$PeptideModifiedSequence)

  feature.data$PeptideModifiedSequence <- droplevels(feature.data$PeptideModifiedSequence)

  # Merge features data with training data
  data.merged = feature.data %>%
    left_join(training.data.tmp,
      by = c("File","FileName","PeptideModifiedSequence","PrecursorCharge",
      "FragmentIon","ProductCharge"))

  data.training.feature <- data.merged %>% filter(!is.na(Status))

  return(list(data.merged = data.merged,feature.data = feature.data,training.data = training.data,data.training.feature = data.training.feature))
}

#' Train a binary classification model to flag peaks with poor chromatography or interference.
#'
#' The function acts as a wrapper to several functions from the caret package to train and optimize a binary predictive peak QC model for the provided training data. Twenty percent of the training dataset is randomly selected as validation set and left out from the training process to estimate the performance of the models on unseen data. The features are mean centered and scaled by diving by the standard deviation before being used for training. Repeated 10-fold cross validation (3 repeats) is applied to the remainder of the training set to minimize over-fitting. The model offering the highest accuracy is used and returned by the function.
#'
#' @param data.merged A dataframe that contains peak identifiers (File,FileName,PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge), the calculated QC metrics as well as the Status assigned by the expert analyst to each transition pair. data.merged is usually output of MakeDataSet function (output$data.merged or output$data.training.feature).
#' @param response.var This variable indicates the name of the column that stores the "ok" and "flag" labels for the transition pairs in the training data.
#' @param description.columns If the input dataframe contains columns corresponding to description variables (such as Notes), it should be indicated here. Description and identifier columns will be removed from the data before training the model.
#' @param method The machine learning algorithm for training the classifier. The algorithm can be chosen from the list of available packages in caret \url{https://topepo.github.io/caret/available-models.html}. The following have been tested: RRF, regLogistic, svmLinear3, svmPoly, kknn. Before using TrainQCModel with any of these packages, you will need to first install the machine learning package using the install.packages command.
#' @param tuneGrid Use this parameter of you want to specify  tuneGrid for the caret train method. Otherwise, set tuneGrid to NULL. See the caret package help for more details: \url{https://topepo.github.io/caret/model-training-and-tuning.html}.
#' @param random.seed To fix the random seed for splitting the dataset into training and validation and the data splitting for cross validation, provide a vector of length 2 e.g. random.seed = c(1000,2000). This is particularly useful if you want to compare multiple models with the same data split.
#' @param export.model A Logical parameter to indicate whether the model should be saved. If export.model = TRUE the model will be saved in model.path.
#' @param model.path Path to the directory where the model will be saved if export.model = TRUE.
#'
#' @return A list with the following objects:
#'               model: Trained model to flag peaks with poor chromatography or interference.
#'               performance.testing: Confusion matrix of applying the model on the unseen validation data (20% of the input data). This parameter can help evaluate the performance of the model on unseen data and identify potential overfitting issues.
#'               model.file.path: If export.model = TRUE and the model is saved, the path and file name for the model is stored in this field.
#'
#' @export
#'
#' @import caret
#' @import RRF
#'
#' @examples
#'
#' rrf.grid <-  expand.grid(mtry = c(2,10),
#'                          coefReg = c(0.5,1),
#'                          coefImp = c(0))
#' model.rrf <- TrainQCModel(data.set.CSF$data.training.feature,
#'                           response.var = c("Status"),
#'                           description.columns = c("Notes"),
#'                           method = "RRF",
#'                           tuneGrid = rrf.grid,
#'                           random.seed = c(100,200))
#'

TrainQCModel <- function(data.merged, response.var = c("Status"), description.columns = c("Notes"), method = "RRF", metric = c("Accuracy","ROC"),tuneGrid = NULL, random.seed = NULL, export.model = FALSE, model.path = "", ...) {

  identifier.columns = c("File","FileName","PeptideModifiedSequence","FragmentIon",
                         "IsotopeLabelType","PrecursorCharge","ProductCharge")

  # change data table to data frame. Otherwise you need to add with = FALSE to slice subsets.
  data.merged <- data.frame(data.merged)

  # The identifier columns are removed.
  data.merged.feature.only = data.merged[,colnames(data.merged)
                                         [!(colnames(data.merged) %in% identifier.columns)]]

  # The description columns are removed
  data.merged.feature.only = data.merged.feature.only[,colnames(data.merged.feature.only)
                                         [!(colnames(data.merged.feature.only) %in% description.columns)]]


  # The response variable columns are removed to leave only the features
  data.merged.feature.only = data.merged.feature.only[,colnames(data.merged.feature.only)
                                                      [!(colnames(data.merged.feature.only) %in% response.var)]]

  # response vector
  resp_vector <- as.matrix(data.merged[,response.var])

  # a 10-fold repeated cross validation (3 repeats) is used for model optimization.
  # if metric = ROC, the trainControl
  if (metric == "ROC") {
    train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,verboseIter = FALSE,sampling = "up",grid,classProbs = TRUE,summaryFunction = twoClassSummary)
  } else {
    train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,verboseIter = FALSE,sampling = "up",grid)
  }

  # split the data into training and testing sets. Training set is used to optimize model parameters and train the model. Testing set is used as unseen data to estimate model performance and evaluate overfitting
  if (!is.null(random.seed)) set.seed(random.seed[1])
  trainIndex <- createDataPartition(as.matrix(data.merged[,response.var]),
                                    p = .8, list = FALSE, times = 1)

#  datasetTrain <- data.merged.feature.only.transformed[trainIndex,]
  datasetTrain <- data.merged.feature.only[trainIndex,]

  # if a seed is provided for controlling the randomness of the algorithm, apply here
  if (!is.null(random.seed)) set.seed(random.seed[2])

  # Train the model using the training dataset
  if (!is.null(tuneGrid)) {
    model <- train(as.matrix(datasetTrain),
                   resp_vector[trainIndex],
                   method = method,
                   preProcess = c("center","scale"),
                   trControl = train_control,
                   importance = TRUE,
                   metric = metric,
                   tuneGrid = tuneGrid)

  } else {
      model <- train(as.matrix(datasetTrain),
                     resp_vector[trainIndex],
                     method = method,
                     preProcess = c("center","scale"),
                     trControl = train_control,
                     importance = TRUE,
                     metric = metric)

  }

  # predicting the response on inputs
  response.prediction <- predict(model, newdata = data.merged.feature.only)

  # model performance based on confusion matrix for unseen data (testing set)
  performance.testing <- confusionMatrix(as.factor(response.prediction[-trainIndex]),as.factor(resp_vector[-trainIndex]))

  QC.model = list(model = model, performance.testing = performance.testing, model.file.path = NA)

  # if export.features is true save the features in the csv file specified by feature.path
  if (export.model == TRUE) {


    # template file name
    model.file <- file.path(model.path,paste0("model_",method,format(Sys.time(), "_%Y%m%d_%H%M"),".rda"))

    # if the template path does not exist, create it
    if (!dir.exists(model.path)) {
      dir.create(model.path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }

    QC.model$model.file.path <- model.file

    save(QC.model, file = model.file)
  }
  return(QC.model)

}

#' Apply a binary predictive peak QC model to flag peaks that suffer from poor chromatography or interference.
#'
#' The function takes the feature dataframe (output of ExtractFeatures output$features) and the trained predictive QC model (output of TrainQCModel) and applies the model to the input feature data using the function provided in the caret package.
#'
#' @param data.feature A dataframe that contains peak identifiers (File,FileName,PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge) as well as  QC metrics calcualted for each transition pair. This dataframe is the output of ExtractFeatures function (output$features).
#' @param model  The predictive model of peak QC. The model is the output of TrainQCModel.
#' @param response.var If the input dataframe contains columns corresponding to response variables, it should be indicated here. it should be indicated here. Response and description columns as well as identifier columns will be removed from the data before applying the model.
#' @param description.columns If the input dataframe contains columns corresponding to description variables (Such as Notes), it should be indicated here. Response and description columns as well identifier columns will be removed from the data before applying the model.
#' @param flag.prob.threshold A numeric value between 0 and 1 which determines the cut-off threshold for assigning classes to each peak based on corresponding class probabilitis. By default, the caret package uses a probability threshold of 0.5. This parameter can be used to overrise the default probability threshold.
#' @param standard.intensity.threshold This parameter can be used to set an intensity threshold to identify and flag transitions where the spiked-in standard is too low. If the numerical value of desired intensity threshold is provided, it is used to flag any transition whose standard signal intensity is below this threshold. For such transitions, this will override the model output.
#' @param type If type = "prob", the function will return class probabilities for the binary classification. This feature can be used only if the model supports classification probabilities e.g. logistic regression and random forest.
#'
#'
#' @return A dataframe of the predicted response (final class and/or class probabilities)  appended to the input data.feature.
#'
#' @export
#'
#' @import caret
#'
#' @examples
#'
#' response.data <- ApplyQCModel(data.set.CSF$feature.data,
#'                               model.rrf.CSF,
#'                               response.var = c("Status"),
#'                               description.columns = c("Notes"),
#'                               flag.prob.threshold = 0.5,
#'                               type = "prob")


ApplyQCModel <- function(data.feature, model, response.var = c("Status"), description.columns = c("Notes"), flag.prob.threshold = 0.5, standard.intensity.threshold = NULL, type = NULL, ...) {

  # purpose: Apply a binary classification model to flag peaks that require manual inspection
  #
  # args:
  #   data.feature: A dataframe that contains the peak identifiers (File,FileName,PeptideModifiedSequence,FragmentIon,IsotopeLabelType,PrecursorCharge and ProductCharge) as well as the QC metrics calcualted for each. data.feature is the output of ExtractFeatures$features
  #   model: the peak binary classifier. The model is the output of TrainQCModel. We recommend to use the default model provided with the package.
  #   response.var: If the data.feature includes any response variables, it should be indicated here. Response and description columns as well identifier columns will be removed from the data before applying the model
  #   description.columns: If the data.feature includes any description columns, it should be indicated here. Response and description columns as well identifier columns will be removed from the data before applying the model
  #   standard.intensity.threshold: If not NULL, it should be a numerical value that is used to flag any transition whose signal intensity is below this threshold. For such transitions, this will override the model output.
  #
  # returns:
  #   data.feature.pred: the predicted response is appended to the input data.feature by a column


  # error and warning handling ---------------------------------------

  #   input errors: if the input standard.intensity.threshold is not numeric
  error.input.format = simpleError("ApplyQCModel: standard.intensity.threshold must be numeric")

  # error handling
  if (!is.null(standard.intensity.threshold) & !is.numeric(standard.intensity.threshold)) stop(error.input.format)

  # function body  ---------------------------------------

  # default identifier.columns, which will be removed from data.merged
  identifier.columns = c("File","FileName","PeptideModifiedSequence","FragmentIon",
                         "IsotopeLabelType","PrecursorCharge","ProductCharge")

  # change data table to data frame. Otherwise you need to add with = FALSE to slice subsets.
  data.feature <- data.frame(data.feature)


  # The identifier columns are removed.
  data.merged.feature.only = data.feature[,colnames(data.feature)
                                         [!(colnames(data.feature) %in% identifier.columns)]]

  # The description columns are removed
  data.merged.feature.only = data.merged.feature.only[,colnames(data.merged.feature.only)
                                                      [!(colnames(data.merged.feature.only) %in% description.columns)]]

  # The response variable columns are removed to leave only the features
  data.merged.feature.only = data.merged.feature.only[,colnames(data.merged.feature.only)
                                                      [!(colnames(data.merged.feature.only) %in% response.var)]]

  # predicting the response on inputs
  response.prediction.prob <- predict(model$model, newdata = data.merged.feature.only,type = "prob")

  # apply the class probability threshold for flag provided by the user
  response.prediction <- data.frame(Status.prediction = rep("flag",nrow(data.merged.feature.only)))
  levels(response.prediction$Status.prediction) <- c("flag","ok")
  response.prediction$Status.prediction[response.prediction.prob[,1] <= flag.prob.threshold] <- "ok"

  # merge QC predictions with original data
  data.feature.pred <- cbind(data.feature,response.prediction)

  if (!is.null(type) && type == "prob") {

    response.prediction.prob <- data.frame(flag.prob.prediction = response.prediction.prob[,1],
                                           ok.prob.prediction = response.prediction.prob[,2])
    # merge QC predictions with original data
    data.feature.pred <- cbind(data.feature.pred,response.prediction.prob)

  }


  # if a threshold for standard signal intensity is provided, flag any transition whose max intensity is below the threshold:
  if (!is.null(standard.intensity.threshold)) {
    data.feature.pred[data.feature.pred$TransitionMaxIntensity_standard < standard.intensity.threshold,"Status.prediction"] <- "flag"
  }

  return(data.feature.pred)

}
