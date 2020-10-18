

###################### I realized the permutation/ agnostic interpretation (NULL_model) in the above example models cross platforms (C, R, Python) might do not work for some machines (maninly becasue the returned optimal model saved in R and then, permute variables using a C language then transform the model again to Python using a parallel C function, it usually speeds up the computation and saves time, however, if you use a lower version of tensorflow, as I stated in the begaining, it doesn't work for parallel computations). An alternative choice if your machine and envrionmnet can't work with the data transformation between C, R and Python.
#### You can choice to use keras for Neural Network(NN) API & lime for LIME(Local Interpretable Model-agnostic Explanations).
#### An example shows below (without parallel computation)

## create your recipes
library(recipes)
rec_obj <- 
  recipe (Polarity ~ ., data = train_tbl) %>%
  # (1) Log-Transform: 'Topics'-----------------------------------------
step_log (Dining_1, Smoking_2, Breakfast_3, Spa_4, Style_Decoration_5,
          `Room View_6`, Location_7, `Service Encounter_8`, Homeliness_9,
          Bathroom_10, Pool_11, Attraction_12, Communication_13, Wedding_14,
          Amenity_15, Casino_16, Parking_17) %>%
  # (2) One-Hot Encoding: 'ID_HOTEL' (Categorical Var.)-----------------
#step_dummy (all_nominal(), -all_outcomes()) %>%
# (3) Normalizing Vars------------------------------------------------
step_center (all_predictors(), -all_outcomes()) %>%
  step_scale (all_predictors(), -all_outcomes()) %>%
  # (4) Last Step: prep()-----------------------------------------------
prep (data = train_tbl)

save (rec_obj, file = 'object/rec_obj.Rda')

### baking your recipes
Predictors

x_train_tbl <- 
  bake (rec_obj, newdata = train_tbl) %>% 
  select (-Polarity)

x_test_tbl <- 
  bake (rec_obj, newdata = test_tbl) %>% 
  select (-Polarity)

glimpse (x_train_tbl)

#We need to store the actual values (truth) as y_train_vec & y_test_vec, for modeling our ANN. We convert to a series of numeric ones and zeros which can be accepted by the Keras modeling functions: Add vec to the name so we can easily remember the class of the object.

#Response variables for training and testing sets
Convert Positive => 1, and Negative => 0

y_train_vec <- ifelse (pull (train_tbl, Polarity) == "Positive", 1, 0)

y_test_vec  <- ifelse (pull (test_tbl, Polarity) == "Positive", 1, 0)

Building our Artificial Neural Network in R

model_keras <- keras_model_sequential()

model_keras %>% 
  # (1) 1st Hidden Layer-------------------------------------------------
layer_dense (units              = 16, #=> Num Of Nodes
             kernel_initializer = "uniform", 
             activation         = "relu",    
             input_shape        = ncol(x_train_tbl)) %>% 
  layer_dropout (rate = 0.1) %>%  #=> Dropout Below 10%: Prevent overfitting
  # (2) 2nd Hidden Layer-------------------------------------------------
layer_dense (units              = 16,
             kernel_initializer = "uniform", 
             activation         = "relu") %>% 
  layer_dropout (rate = 0.1) %>%  
  # (3) Output Layer-----------------------------------------------------
layer_dense (units              = 1, #=> Binary/Multi?=>That Number
             kernel_initializer = "uniform", 
             activation         = "sigmoid") %>% #=> Common for Binary
  # (4) Compile Model-----------------------------------------------------
compile (optimizer = 'adam', #=> Most Popular for Optimization Algo.
         loss      = 'binary_crossentropy', #=> Binary Classification
         metrics   = c('accuracy') ) #=> Train/Test Evaluation

# Check
model_keras

#Fit ‘keras_model’ to the ‘Training’ Data

system.time ( 
  history <- fit (
    object           = model_keras,             # => Our Model
    x                = as.matrix (x_train_tbl), #=> Matrix
    y                = y_train_vec,             #=> Numeric Vector 
    batch_size       = 50,     #=> #OfSamples/gradient update in each epoch
    epochs           = 35,     #=> Control Training cycles
    validation_split = 0.30) ) #=> Include 30% data for 'Validation' Model

(1) Print a summary of the [Training] history

#We want MINIMAL difference between the Validation & Training accuracy.

print (history)



#(1) Predicted Class

yhat_keras_class_vec <- 
  predict_classes (object = model_keras, 
                   x = as.matrix(x_test_tbl)) %>%
  as.vector()

#(2) Predicted Class “Probability”

yhat_keras_prob_vec <- 
  predict_proba (object = model_keras, 
                 x = as.matrix(x_test_tbl)) %>%
  as.vector()

#Format test data and predictions for yardstick metrics

estimates_keras_tbl <- tibble(
  truth      = as.factor(y_test_vec) %>% 
    fct_recode (Positive = "1", Negative = "0"),
  estimate   = as.factor(yhat_keras_class_vec) %>% 
    fct_recode (Positive = "1", Negative = "0"),
  class_prob = yhat_keras_prob_vec )

options(scipen = 999)

head(estimates_keras_tbl, 10)


###### now connect with lime
Setup

#LIME is not setup out-of-the-box to work with keras. We need to make two custom functions to work properly.

#model_type(): Used to tell lime what type of model we are dealing with. It could be classification, regression, survival, etc.

#predict_model(): Used to allow lime to perform predictions that its algorithm can interpret.

#Identify the class of our model object with class() function.

class (model_keras)



#(1) Setup lime::model_type() for keras

Next we create our model_type() function. The only input is x, the keras model. The function simply returns classification, which tells LIME we are classifying.

model_type.keras.models.Sequential <- function(x, ...) {
  "classification"}

#(2) Setup lime::predict_model() for keras

Now we can create our predict_model(): Wraps keras::predict_proba().

The trick here is to realize that it’s inputs must be

x a model,
newdata a dataframe object (this is important), and,
type which is not used but can be use to switch the output type.

The output is also a little tricky because it must be in the format of probabilities by classification (this is important; shown next).

predict_model.keras.models.Sequential <- function (x, newdata, type, ...) {
  pred <- predict_proba (object = x, x = as.matrix(newdata))
  data.frame (Positive = pred, Negative = 1 - pred) }

Test predict_model()

#Run this next script to show you what the output looks like and to test our predict_model() function. See how it’s the probabilities by classification. It must be in this form for model_type = "classification"

predict_model (x       = model_keras, 
               newdata = x_test_tbl, 
               type    = 'raw') %>%
  tibble::as_tibble()


#Run lime() on Training Set

#Now the fun part, we create an explainer using the lime(). Just pass the Training Data set without the Attribution column. The form must be a data frame, which is OK since our predict_model() will switch it to an keras object.

#Set model = automl_leader our leader model, and
bin_continuous = FALSE.

#We could tell the algorithm to bin continuous variables, but this may not make sense for categorical numeric data that we didn’t change to factors.

explainer <- lime::lime (
  x              = x_train_tbl, 
  model          = model_keras, 
  bin_continuous = FALSE)

#Run explain() on explainer

#Now we run the explain(), which returns our explanation.

system.time (
  explanation <- lime::explain (
    x_test_tbl[1:10, ], # Just to show first 10 cases
    explainer    = explainer, 
    n_labels     = 1, # explaining a `single class`(Polarity)
    n_features   = 4, # returns top four features critical to each case
    kernel_width = 0.5) ) # allows us to increase model_r2 value by shrinking the localized evaluation.

##    user  system elapsed 
##  10.972   1.148  11.781

#Feature Importance Visualization

#The payoff for LIME is feature importance plot. This allows us to visualize each of the first ten cases (observations) from the test data. The top four features for each case are shown. Note that they are not the same for each case.

#Green Bars: The feature supports the model conclusion
#Red Bars contradict.

plot_features (explanation) +
  labs (title = "LIME: Feature Importance Visualization",
        subtitle = "Hold Out (Test) Set, First 10 Cases Shown")


##Another excellent visualization is plot_explanations(): Facetted heatmap of all case/label/feature combinations.

#It’s a more condensed version of plot_features(), but we need to be careful because it does not provide exact statistics and it makes it less easy to investigate binned features (Notice that “tenure” would not be identified as a contributor even though it shows up as a top feature in 7 of 10 cases).

plot_explanations (explanation) +
  labs (title = "LIME Feature Importance Heatmap",
        subtitle = "Hold Out (Test) Set, First 10 Cases Shown")

#Check Explanations With Correlation Analysis

#One thing we need to be careful with the LIME visualization is that we are only doing a sample of the data, in our case the first 10 test observations. Therefore, we are gaining a very localized understanding of how the ANN works.
#However, we also want to know on from a global perspective what drives feature importance.

#Feature Correlations to Polarity

#We can perform a correlation analysis on the Training set as well to help glean What Features Correlate Globally to Polarity. We’ll use the corrr package, which performs tidy correlations with the function correlate().

( corrr_analysis <- x_train_tbl %>%
    mutate (Polarity = y_train_vec) %>%
    correlate () %>%
    focus (Polarity) %>%
    rename (feature = rowname) %>%
    arrange (abs(Polarity)) %>%
    mutate (feature = as_factor(feature)) )

#Correlation Visualization

#The correlation visualization helps in Distinguishing Which Features are Relavant to Polarity.

corrr_analysis %>%
  
  ggplot (aes (x = Polarity, y = fct_reorder(feature, desc(Polarity)))) +
  geom_point () +
  # Positive Correlations - Contribute to Polarity--------------------------------------------
geom_segment (aes(xend = 0, yend = feature), 
              color = palette_light()[[2]], 
              data = corrr_analysis %>% filter(Polarity > 0)) +
  geom_point (color = palette_light()[[2]], 
              data = corrr_analysis %>% filter(Polarity > 0)) +
  # Negative Correlations - Prevent Polarity--------------------------------------------------
geom_segment (aes(xend = 0, yend = feature), 
              color = palette_light()[[1]], 
              data = corrr_analysis %>% filter(Polarity < 0)) +
  geom_point (color = palette_light()[[1]], 
              data = corrr_analysis %>% filter(Polarity < 0)) +
  # Vertical lines-------------------------------------------------------------------------
geom_vline (xintercept = 0, color = palette_light()[[5]], size = 1, linetype = 2) +
  geom_vline (xintercept = -0.25, color = palette_light()[[5]], size = 1, linetype = 2) +
  geom_vline (xintercept = 0.25, color = palette_light()[[5]], size = 1, linetype = 2) +
  # Aesthetics-----------------------------------------------------------------------------
theme_tq () +
  labs (title = "Polarity Correlation Analysis",
        subtitle = "Positive Correlations vs. Negative Correlations", 
        y = "Feature Importance")





### typical example using some ready default models, and other cv method

### grid search

layer1 = c(6,12,18,24,30)
layer2 = c(6,12,18,24,30)
layer3 = c(6,12,18,24,30)

cv.folds = 5

# In order to make models fully reproducible when using parallel processing, we need to pass seeds as a parameter
# https://stackoverflow.com/questions/13403427/fully-reproducible-parallel-models-using-caret

total.param.permutations = length(layer1) * length(layer2) * length(layer3)

seeds <- vector(mode = "list", length = cv.folds + 1)
set.seed(1)  
for(i in 1:cv.folds) seeds[[i]]<- sample.int(n=1, total.param.permutations, replace = TRUE)
seeds[[cv.folds + 1]]<-sample.int(1, 1, replace = TRUE) #for the last model

nn.grid <- expand.grid(layer1 = layer1, layer2 = layer2, layer3 = layer3)

cl <- makeCluster(detectCores()*0.5) # use 50% of cores only, leave rest for other tasks
registerDoParallel(cl)

train_control <- DeepGenomeScanControl(method = "cv" 
                                       ,number=cv.folds 
                                       ,seeds = seeds # user defined seeds for parallel processing
                                       ,verboseIter = TRUE
                                       ,allowParallel = TRUE
)

stopCluster(cl)
registerDoSEQ()

tic("Total Time to NN Training: ")
set.seed(1)
model.nn.caret = DeepGenomeScan(form = formula,
                                data = scaled.train.data,
                                method = 'neuralnet',
                                tuneGrid = nn.grid,
                                trControl = train_control
)
toc()
















######## future work#####################################

##### last but not least, there is a good way to train and compare more models, and chose one best model, not only constraining on neural networks,  called unified automatic machine learning framework.
#### I will incroporate it into the framework and test them in the near future.
##Run AutoML for 100 models,
####The Automatic Machine Learning (AutoML) function automates the supervised machine learning model training process. The current version of AutoML trains and cross-validates a Random Forest, an Extremely-Randomized Forest, a random grid of Gradient Boosting Machines (GBMs), a random grid of Deep Neural Nets, and then trains a Stacked Ensemble using all of the models.
aml <- h2o.automl(x = x, y = y,
                  training_frame = train,
                  max_models = 100,
                  seed = 1)

# View the AutoML Leaderboard
lb <- aml@leaderboard
print(lb, n = nrow(lb))

# Get model ids for all models in the AutoML Leaderboard
model_ids <- as.data.frame(lb$model_id)[,1]

# View variable importance for all the models (besides Stacked Ensemble)
for (model_id in model_ids) {
  print(model_id)
  m <- h2o.getModel(model_id)
  h2o.varimp(m)
  h2o.varimp_plot(m)
}