#frist analysis with normalizated values
library(scales)
# scientific annotation function
scientific_annot_fun <- function(l) {
  l <- format(l, scientific = T)
  l <- gsub("^(.*)e.", "10^", l)
  parse(text = l)
}

# log2(x+1) function
log2p <- function(x){
  log2(x + 1)
}
log2p_inv <- function(x){
  (2^x) - 1
}
log2p_trans <- function() trans_new("log2p", log2p, log2p_inv)

# Pre-processing
library(phyloseq)

rrna <- F
if(rrna){
  data <- readRDS("profiled.rds")
  data <- tax_glom(data, "Genus")
  f <- filterfun_sample(topk(100))
  data <- genefilter_sample(data, f, A = .9*nsamples(data)) %>% prune_taxa(data)
  data <- prune_samples(sample_sums(data) > 1000, data)
}else{
  data <- readRDS("profiled.rds")
  data <- tax_glom(data, "s")
  data <- filter_taxa(data, function(x) sum(x > 0) >= .5*nsamples(data), T)
  data <- filter_taxa(data, function(x) sum(x > 300) > 1, TRUE)
}

# Norm
scaled <- scaleData(otu_table(data), zeroTolerant = T)
scaled <- data.frame(t(scaled))

train <- sample(1:nrow(scaled), size = .8*nrow(scaled))

valid <- scaled[-train,]
train <- scaled[train,]

#ga for all vars
predictions <- function(train) {
   pred <- lapply(seq_along(train), function(i) {
     data_x <- train[,-i]
     data_y <- train[[i]]

     ga1 <- gadualModel(data_x = data_x, data_y = data_y, train)
     return(ga1)
   })
   names(pred) <- names(train)
   pred
 }

ga1 <- predictions(train)

# obs and pred values for all vars
df <- function(ga1) {
  df <- lapply(seq_along(ga1), function(i) {
    data <- data.frame(ga1[[i]]$pred, valid[[i]])
    names(data)[1] <- "pred"
    names(data)[2] <- "obs"
    return(data)
  })
  names(df) <- names(train)
  df
}
df <- df(ga1)

library(dplyr)
datas <- bind_rows(df, .id = "sp")

ggplot(data = datas, aes(x = obs, y = pred)) +
  geom_point(alpha = .3) +
  ggplot2::scale_y_continuous(trans = "log2p",
                     breaks = base::c(0, 10, 100, 1000, 10000, 100000),
                     labels = scientific_annot_fun) +
  ggplot2::scale_x_continuous(trans = "log2p",
                     breaks = base::c(0, 10, 100, 1000, 10000, 100000),
                     labels = scientific_annot_fun) +
  ggplot2::facet_wrap(~ sp) +
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ggplot2:: geom_abline(intercept = 0, slope = 1, colour = "red", size = .1) +
  coord_flip() +
  ggsave("prediction_norm.pdf", dpi = 600, height = 6, width = 9)

library(foreach)
library(doMC)

old <- getDoParWorkers()
registerDoMC(2)
toPredict <- colnames(train)
res <- foreach(v=toPredict) %dopar% {
  # Computing variable importance for each
  # variable in the dataset
  y <- train[,v]
  x <- train[,colnames(train) != v]

  model <- trainmod(x = x, y = y)

  # Get test data
  obs <- validation[,v]
  newdata <- validation[,colnames(validation) != v]

  # Get variable importance based on test data
  imp <- variableImp(model, obs, newdata, shuffles = 100)
  r2 <- getR(model, obs, newdata)
  list(imp=imp, r2 = r2)
}
registerDoMC(old)

# Extracting variable importance
# and building a matrix. The importance of
# predictors is reported in each column.
# The name of the column is the name of
# the predicted variables whereas the values
# in each row correspond to the importance of
# the variable. If I want to know the importance
# of the 5th variable in predicting the
# 1st i have to extract the value in the 1st column
# at the 5th row.


cost <- sapply(seq_along(res), function(i) {
  r <- res[[i]]$imp[colnames(train)]
  names(r) <- colnames(train)
  r
})
colnames(cost) <- rownames(cost)

# Overall importance as predictors
suggestion <- order(rowMeans(cost, na.rm = T), decreasing = T)

sort(rowMeans(cost, na.rm = T))

# Getting R-squared values
r2 <- sapply(res, "[[", "r2")
names(r2) <- colnames(cost)

simple_fitness <- function(i, cost){
  c <- cost[i,i]
  gain <- sum(c[upper.tri(c, diag = F)])
  loss <- sum(c[lower.tri(c, diag = F)])
  gain - loss
}

fitness_custom <- function(i, cost, r2){
  ordered <- data.frame(t(combn(colnames(cost)[i], 2)), stringsAsFactors = F)
  keep <- apply(ordered, 1, function(x) cost[x[1], x[2]] > 0)
  ordered <- ordered[keep,]

  from <- unique(ordered[[1]])
  to <- unique(ordered[[2]])

  roots <- from[!from %in% to]

  penalty <- sum(log(r2[!names(r2) %in% roots]))
  c <- cost[i,i]

  gain <- sum(c[upper.tri(c, diag = F)])
  loss <- sum(c[lower.tri(c, diag = F)])
  (gain - loss) + penalty
}

set.seed(11)
ga_1 <- ga(fitness = simple_fitness, # custom fitness function
           type = "permutation", # optimization data type
           lower = 1,
           upper = ncol(train),
           cost = cost,
           # r2 = r2,
           monitor = plot,
           suggestions = matrix(suggestion, nrow=1),
           popSize = 10000,
           maxiter = 100,
           elitism = 1000,
           keepBest = T,
           parallel = T,
           pmutation = 0.2,
           run = 20,
         seed = 11)

best <- data.frame(taxa=colnames(train)[ga_1@solution[1,]],
                   position=1:ncol(train))

lapply(ga_1@bestSol, function(x){
  data.frame(taxa = colnames(train)[x],
             position = 1:ncol(train))
}) %>%
  bind_rows() %>%
  group_by(taxa, position) %>%
  tally() %>%
  ungroup() %>%
  group_by(taxa) %>%
  mutate(perc = n/sum(n),
         ord = median(rep(position, n))) %>%
  ungroup() %>%
  mutate(taxa = reorder(taxa, ord, unique)) %>%
  ggplot(aes(x = position, y = taxa, fill = perc)) +
  geom_tile() +
  geom_text(data = best, aes(x = position, y = taxa, label = "+"),
            inherit.aes = F) +
  scale_fill_distiller(palette = "Spectral") +
  coord_cartesian(expand = F) +
  theme(panel.background = element_rect(fill = "black"),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = 1:ncol(train))

# do.call(rbind, positions) %>%
#   group_by(taxa, position) %>%
#   tally() %>%
#   ungroup() %>%
#   group_by(taxa) %>%
#   mutate(perc = n/sum(n),
#          ord = mean(rep(position, n))) %>%
#   ungroup() %>%
#   mutate(taxa = reorder(taxa, ord, unique)) %>%
# ggplot(aes(x = position, y = taxa, fill = perc)) +
#   geom_tile() +
#   scale_fill_distiller(palette = "Spectral") +
#   coord_cartesian(expand = F) +
#   theme(panel.background = element_rect(fill = "black"),
#         panel.grid = element_blank()) +
#   scale_x_continuous(breaks = 1:ncol(train))

ga_2 <- ga(fitness = function(i, costs) {
  c <- costs[i,i]
  gain <- sum(c[upper.tri(c, diag = F)])
  loss <- sum(c[lower.tri(c, diag = F)])
  gain - loss
}, # custom fitness function
type = "permutation", # optimization data type
lower = 1,
upper = ncol(train),
costs = cost,
monitor = plot,
popSize = 1000,
maxiter = 100,
elitism = 100,
keepBest = T,
parallel = T,
pmutation = 0.2,
run = 10)



data.frame(taxa=colnames(train)[ga_2@solution[1,]],
           position=1:ncol(train))

ordered <- colnames(train)[ga_2@solution[1,]]
ordered <- setNames(data.frame(t(combn(ordered, 2)), stringsAsFactors = F),
                    base::c("from", "to"))

colnames(cost) <- rownames(cost)
names(r2) <- colnames(cost)
apply(ordered, 1, function(x) {
  penalty <- log(r2[x[2]])
  c <- cost[x[1], x[2]]
  c
  # c + penalty
  })
ordered <- ordered[apply(ordered, 1, function(x) cost[x[1] ,x[2]]) > 0,]

graph <- empty.graph(colnames(train))
arcs(graph) <- ordered

graphviz.plot(graph)



#obs vs pred without normalization
nonorm <- data.frame(t(otu_table(data)))
train <- sample(1:nrow(nonorm), size = .8*nrow(nonorm))

validation <- prova[-train,]
train <- prova[train,]

ga2 <- predictions(train)

df <- function(ga1) {
  df <- lapply(seq_along(ga1), function(i) {
    data <- data.frame(ga1[[i]]$pred, valid[[i]])
    names(data)[1] <- "pred"
    names(data)[2] <- "obs"
    return(data)
  })
  names(df) <- names(train)
  df
}
df_2 <- df(ga2)
dt <- bind_rows(df_2, .id = "sp")

ggplot(data = dt, aes(x = obs, y = pred)) +
  geom_point(alpha = .3) +
  ggplot2::scale_y_continuous(trans = "log2p",
                              breaks = base::c(0, 10, 100, 1000, 10000, 100000),
                              labels = scientific_annot_fun) +
  ggplot2::scale_x_continuous(trans = "log2p",
                              breaks = base::c(0, 10, 100, 1000, 10000, 100000),
                              labels = scientific_annot_fun) +
  ggplot2::facet_wrap(~ sp) +
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  ggplot2:: geom_abline(intercept = 0, slope = 1, colour = "red", size = .1) +
  coord_flip() +
  ggsave("prediction.pdf", dpi = 600, height = 6, width = 9)

