detach("package:relome", unload=TRUE)
data("iris")
results <- runRelome(data = iris, interest = c("Species","Petal.Width"), threshold=1,
                     adjust = "BH", zoom.p = TRUE, verbose=TRUE, plot=TRUE,
                     mfrow = c(2,3), width = 15, height = 10, save.all=TRUE, rerun.all=TRUE,
                     bg = list(scatter="yellow",category=gray.colors(2)), 
                     col = list(scatter="red",category=gray.colors(2)),
                     pch=21, inpdf = FALSE)

# results <- runRelome(data = dat_bc, interest = interest.all, threshold=0.5,
#                      adjust = "BH", adjust.by.var = FALSE,
#                      zoom.p = TRUE, verbose=TRUE, plot=TRUE,
#                      mfrow = c(4,6), width = 15, height = 10, save.all=TRUE, rerun.all=TRUE,
#                      col = list(scatter="red",category=gray.colors(2)))
