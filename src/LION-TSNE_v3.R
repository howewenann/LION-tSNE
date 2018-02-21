library(ProjectTemplate)
load.project()

set.seed(100)

data_Y <- tsne.3D 
data_X <- train
data_X <- data_X[, -1]

s_ratio = 0.9
train_rows <- sample(c(1:nrow(data_Y)), 
                     s_ratio * nrow(data_Y), 
                     replace = F)

data_Y_train <- data_Y[train_rows, ]
data_Y_test <- data_Y[-train_rows, ]

data_X_train <- data_X[train_rows, ]
data_X_test <- data_X[-train_rows, ]


Y = LION_tSNE(X_new = data_X_test, Y_train = data_Y_train, X_train = data_X_train, subsample = 0.1)











