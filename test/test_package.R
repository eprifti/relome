data("iris")
results <- runRelome(data = iris, interest = c("Species","Petal.Width"), threshold=1,
                       adjust = "BH", zoom.p = TRUE, verbose=TRUE, plot=TRUE,
                       mfrow = c(4,6), width = 15, height = 10, save.all=TRUE, rerun.all=TRUE,
                       col = list(scatter="gold2",category=gray.colors(2)))

