library(data.table)
library(Rtsne)
library(quantmod)
library(ggplot2)
library(plotly)
library(plot3D)

set.seed(100)

setwd("D:/NCS/Experiments/TSNE")

train_data <- fread('train.csv')

train_label <- train_data[, 1]
train_vars <- train_data[, -1]
# train_vars <- 
#     data.frame(
#         do.call(cbind, 
#                 lapply(train_vars, function(x) ifelse(x > 0, 1, 0)
#                 )
#         )
#     )


sample_data <- data.frame(train_label, train_vars)
sample_data <- unique(sample_data)
# sample_data <- sample_data[sample(c(1:nrow(train_vars)), 10000), ]

train_tsne <- Rtsne(sample_data[, -1],
                    initial_dims = 50,
                    check_duplicates = F, 
                    dims = 2, 
                    perplexity = 30,
                    max_iter = 1000,
                    verbose = T)


## for plotting
colors = rainbow(length(unique(sample_data[,1])))
names(colors) = unique(sample_data[,1])


# plot1 <- 
#     ggplot(data = data.frame(train_tsne$Y),
#            mapping = aes(x = X1,
#                          y = X2,
#                          color = as.character(sample_data[,1]))) + 
#     geom_point() + 
#     labs(x = "X1", y = 'X2', color = 'Number')

# plot1 <- ggplotly(plot1)

plot1 <-  plot_ly(data.frame(train_tsne$Y),
                  x = ~X1, y = ~X2, 
                  color = as.character(sample_data[,1]),
                  colors = "Paired",
                  marker = list(size = 2)) %>%
    add_markers()

plot1

fwrite(data.frame(train_tsne$Y), 'tsne_2D.csv', row.names = F)

# Try 3D plot
train_tsne2 <- Rtsne(sample_data[, -1],
                     initial_dims = 50,
                     check_duplicates = F, 
                     dims = 3, 
                     perplexity = 30,
                     max_iter = 1000,
                     verbose = T)


plot2 <- plot_ly(data.frame(train_tsne2$Y),
                 x = ~X1, y = ~X2, z = ~X3, 
                 color = as.character(sample_data[,1]),
                 colors = "Paired",
                 marker = list(size = 2)) %>%
    add_markers()

plot2

fwrite(data.frame(train_tsne2$Y), 'tsne_3D.csv', row.names = F)



# Try PCA
# train_pca <- princomp(sample_data[, -1])
# train_pca_data <- data.frame(train_pca$scores)
# 
# plot(train_pca_data[, c(1,2)], 
#      pch = as.character(sample_data[,1]),
#      col = colors[as.character(sample_data[,1])])
# 
# 
# train_pca$sdev[1:5]

