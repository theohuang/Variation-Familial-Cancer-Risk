## OE Simulation for Frailty Analysis (large samples)
## Last updated: July 6, 2020

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


nboot <- 20
w.list <- seq(-2, 2, 0.1)


cnames <- c("W", "O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC",
            paste(c("O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC"), ".cp", sep = ""),
            paste(c("O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC"), ".nocens", sep = ""),
            paste(c("O.C", "E.C", "OE.C", "O.NC", "E.NC", "OE.NC"), ".cp.nocens", sep = ""))

## frailty on both carriers and non-carriers
res.oe.cnc <- array(NA, c(50 * nboot, length(w.list), length(cnames)))
ct <- 0
for(i in 1:50){
  load(paste(getwd(), "/OE Simulation/Both/oe_sim_cnc_", i, ".RData", sep = ""))
  for(j in 1:20){
    ct <- ct + 1
    res.oe.cnc[ct, , ] <- as.matrix(res.oe[[j]])
  }
}

## frailty on carriers only
res.oe.c <- array(NA, c(50 * nboot, length(w.list), length(cnames)))
ct <- 0
for(i in 1:50){
  load(paste(getwd(), "/OE Simulation/Carriers/oe_sim_c_", i, ".RData", sep = ""))
  for(j in 1:20){
    ct <- ct + 1
    res.oe.c[ct, , ] <- as.matrix(res.oe[[j]])
  }
}

## frailty on non-carriers only
res.oe.nc <- array(NA, c(50 * nboot, length(w.list), length(cnames)))
ct <- 0
for(i in 1:50){
  load(paste(getwd(), "/OE Simulation/NonCarriers/oe_sim_nc_", i, ".RData", sep = ""))
  for(j in 1:20){
    ct <- ct + 1
    res.oe.nc[ct, , ] <- as.matrix(res.oe[[j]])
  }
}



# res.oe.mean <- apply(res.oe.tot, 1, mean)

##### Mean of O/E

## Frailty on both carriers and non-carriers
res.oe.mean.cnc <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.mean.cnc[i, j] <- mean(res.oe.cnc[, i, j])
  }
}
res.oe.mean.cnc <- data.frame(res.oe.mean.cnc)
names(res.oe.mean.cnc) <- cnames

## Frailty on carriers only
res.oe.mean.c <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.mean.c[i, j] <- mean(res.oe.c[, i, j], na.rm = TRUE)
  }
}
res.oe.mean.c <- data.frame(res.oe.mean.c)
names(res.oe.mean.c) <- cnames

## Frailty on non-carriers only
res.oe.mean.nc <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.mean.nc[i, j] <- mean(res.oe.nc[, i, j])
  }
}
res.oe.mean.nc <- data.frame(res.oe.mean.nc)
names(res.oe.mean.nc) <- cnames


## median
res.oe.med.cnc <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.med.cnc[i, j] <- median(res.oe.cnc[, i, j])
  }
}
res.oe.med.cnc <- data.frame(res.oe.med.cnc)
names(res.oe.med.cnc) <- cnames

res.oe.med.c <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.med.c[i, j] <- median(res.oe.c[, i, j])
  }
}
res.oe.med.c <- data.frame(res.oe.med.c)
names(res.oe.med.c) <- cnames

res.oe.med.nc <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.med.nc[i, j] <- median(res.oe.nc[, i, j])
  }
}
res.oe.med.nc <- data.frame(res.oe.med.nc)
names(res.oe.med.nc) <- cnames


## variances
res.oe.var <- matrix(NA, length(w.list), length(cnames))
for(i in 1:length(w.list)){
  for(j in 1:length(cnames)){
    res.oe.var[i, j] <- var(res.oe.cnc[, i, j])
  }
}
res.oe.var <- data.frame(res.oe.var)
names(res.oe.var) <- cnames
res.oe.var$W <- w.list




get_plot_limits <- function(plot) {
  gb = ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  c(xmin, xmax, ymin, ymax)
}


#### Frailty on both carriers and non-carriers

p <- ggplot(res.oe.mean.cnc, aes(OE.C.cp, OE.NC.cp)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "A") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_continuous(name = "Frailty")
p.cnc.cp <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_cnc_cp.png", sep = ""))

p <- ggplot(res.oe.mean.cnc, aes(OE.C.cp.nocens, OE.NC.cp.nocens)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "B") +
  theme(plot.title = element_text(hjust = 0.5))
p.cnc.cp.nocens <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_cnc_cp_nocens.png", sep = ""))

p <- ggplot(res.oe.mean.cnc, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "C") +
  theme(plot.title = element_text(hjust = 0.5))
p.cnc <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_cnc.png", sep = ""))

p <- ggplot(res.oe.mean.cnc, aes(OE.C.nocens, OE.NC.nocens)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "D") +
  theme(plot.title = element_text(hjust = 0.5))
p.cnc.nocens <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_cnc_nocens.png", sep = ""))


#### Frailty on carriers only

p <- ggplot(res.oe.mean.c, aes(OE.C.cp, OE.NC.cp)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "E") +
  theme(plot.title = element_text(hjust = 0.5))
p.c.cp <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_c_cp.png", sep = ""))

p <- ggplot(res.oe.mean.c, aes(OE.C.cp.nocens, OE.NC.cp.nocens)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "F") +
  theme(plot.title = element_text(hjust = 0.5))
p.c.cp.nocens <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_c_cp_nocens.png", sep = ""))

p <- ggplot(res.oe.mean.c, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "G") +
  theme(plot.title = element_text(hjust = 0.5))
p.c <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_c.png", sep = ""))

p <- ggplot(res.oe.mean.c, aes(OE.C.nocens, OE.NC.nocens)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "H") +
  theme(plot.title = element_text(hjust = 0.5))
p.c.nocens <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_c_nocens.png", sep = ""))


#### Frailty on non-carriers only

p <- ggplot(res.oe.mean.nc, aes(OE.C.cp, OE.NC.cp)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "I") +
  theme(plot.title = element_text(hjust = 0.5))
p.nc.cp <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_nc_cp.png", sep = ""))

p <- ggplot(res.oe.mean.nc, aes(OE.C.cp.nocens, OE.NC.cp.nocens)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "J") +
  theme(plot.title = element_text(hjust = 0.5))
p.nc.cp.nocens <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_nc_cp_nocens.png", sep = ""))

p <- ggplot(res.oe.mean.nc, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "K") +
  theme(plot.title = element_text(hjust = 0.5))
p.nc <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_nc.png", sep = ""))

p <- ggplot(res.oe.mean.nc, aes(OE.C.nocens, OE.NC.nocens)) +
  geom_point(aes(color = W)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "", title = "L") +
  theme(plot.title = element_text(hjust = 0.5))
p.nc.nocens <- p + scale_x_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p)))) +
  scale_y_continuous(limits = c(min(get_plot_limits(p)), max(get_plot_limits(p))))
# ggsave(paste(, "/plot_mean_nc_nocens.png", sep = ""))

grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

grid_arrange_shared_legend(p.cnc.cp, p.cnc.cp.nocens, p.cnc, p.cnc.nocens,
             p.c.cp, p.c.cp.nocens, p.c, p.c.nocens,
             p.nc.cp, p.nc.cp.nocens, p.nc, p.nc.nocens,
             nrow = 3, ncol = 4)
