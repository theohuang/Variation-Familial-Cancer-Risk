## Plotting all 4 frailty model results
## Last updated: July 6, 2020


library(ggplot2)
library(gridExtra)
library(grid)

setwd("/Users/thuang/Dropbox (Partners HealthCare)/BayesMendel Dropbox folders/Variation Familial Cancer Risk")

load("OE_Sim_CNC_Main.RData")
res.oe.sim.cnc <- res.oe.sim
load("OE_Sim_C_Main.RData")
res.oe.sim.c <- res.oe.sim
load("OE_Sim_NC_Main.RData")
res.oe.sim.nc <- res.oe.sim
load("OE_Sim_none_Main.RData")
res.oe.sim.none <- res.oe.sim


## O/E ratios
p.cnc <- ggplot(res.oe.sim.cnc, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.005, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.005, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "Frailty on all") +
  scale_color_continuous(name = "Frailty") +
  theme(plot.title = element_text(hjust = 0.5))

p.c <- ggplot(res.oe.sim.c, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.005, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.005, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "Frailty on carriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.nc <- ggplot(res.oe.sim.nc, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.005, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.005, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "Frailty on noncarriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.none <- ggplot(res.oe.sim.none, aes(OE.C, OE.NC)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.005, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.005, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "No frailty") +
  theme(plot.title = element_text(hjust = 0.5))


## O-E differences
p.cnc.diff <- ggplot(res.oe.sim.cnc, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.002, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.001, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O-E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "Frailty on all") +
  scale_color_continuous(name = "Frailty") +
  theme(plot.title = element_text(hjust = 0.5))

p.c.diff <- ggplot(res.oe.sim.c, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.002, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.0008, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O/E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "Frailty on carriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.nc.diff <- ggplot(res.oe.sim.nc, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.001, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.001, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O/E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "Frailty on noncarriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.none.diff <- ggplot(res.oe.sim.none, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.0013, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.0007, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O/E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "No frailty") +
  theme(plot.title = element_text(hjust = 0.5))


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






grid_arrange_shared_legend(p.cnc, p.c, p.nc, p.none, nrow = 2, ncol = 2)
grid_arrange_shared_legend(p.cnc.diff, p.c.diff, p.nc.diff, p.none.diff, nrow = 2, ncol = 2)



## square plots, same scale
pl1.lo <- min(c(res.oe.sim.cnc$OE.C.lo, res.oe.sim.cnc$OE.NC.lo))
pl1.hi <- max(c(res.oe.sim.cnc$OE.C.hi, res.oe.sim.cnc$OE.NC.hi))
pl2.lo <- min(c(res.oe.sim.c$OE.C.lo, res.oe.sim.c$OE.NC.lo))
pl2.hi <- max(c(res.oe.sim.c$OE.C.hi, res.oe.sim.c$OE.NC.hi))
pl3.lo <- min(c(res.oe.sim.nc$OE.C.lo, res.oe.sim.nc$OE.NC.lo))
pl3.hi <- max(c(res.oe.sim.nc$OE.C.hi, res.oe.sim.nc$OE.NC.hi))
pl4.lo <- min(c(res.oe.sim.none$OE.C.lo, res.oe.sim.none$OE.NC.lo))
pl4.hi <- max(c(res.oe.sim.none$OE.C.hi, res.oe.sim.none$OE.NC.hi))

grid_arrange_shared_legend(p.cnc + scale_x_continuous(limits = c(pl1.lo, pl1.hi)) +
                             scale_y_continuous(limits = c(pl1.lo, pl1.hi)),
                           p.c + scale_x_continuous(limits = c(pl2.lo, pl2.hi)) +
                             scale_y_continuous(limits = c(pl2.lo, pl2.hi)),
                           p.nc + scale_x_continuous(limits = c(pl3.lo, pl3.hi)) +
                             scale_y_continuous(limits = c(pl3.lo, pl3.hi)),
                           p.none + scale_x_continuous(limits = c(pl4.lo, pl4.hi)) +
                             scale_y_continuous(limits = c(pl4.lo, pl4.hi)),
                           nrow = 2, ncol = 2)


pl1.diff.lo <- min(c(res.oe.sim.cnc$OEdiff.C.lo, res.oe.sim.cnc$OEdiff.NC.lo))
pl1.diff.hi <- max(c(res.oe.sim.cnc$OEdiff.C.hi, res.oe.sim.cnc$OEdiff.NC.hi))
pl2.diff.lo <- min(c(res.oe.sim.c$OEdiff.C.lo, res.oe.sim.c$OEdiff.NC.lo))
pl2.diff.hi <- max(c(res.oe.sim.c$OEdiff.C.hi, res.oe.sim.c$OEdiff.NC.hi))
pl3.diff.lo <- min(c(res.oe.sim.nc$OEdiff.C.lo, res.oe.sim.nc$OEdiff.NC.lo))
pl3.diff.hi <- max(c(res.oe.sim.nc$OEdiff.C.hi, res.oe.sim.nc$OEdiff.NC.hi))
pl4.diff.lo <- min(c(res.oe.sim.none$OEdiff.C.lo, res.oe.sim.none$OEdiff.NC.lo))
pl4.diff.hi <- max(c(res.oe.sim.none$OEdiff.C.hi, res.oe.sim.none$OEdiff.NC.hi))

grid_arrange_shared_legend(p.cnc.diff + scale_x_continuous(limits = c(pl1.diff.lo, pl1.diff.hi)) +
                             scale_y_continuous(limits = c(pl1.diff.lo, pl1.diff.hi)),
                           p.c.diff + scale_x_continuous(limits = c(pl2.diff.lo, pl2.diff.hi)) +
                             scale_y_continuous(limits = c(pl2.diff.lo, pl2.diff.hi)),
                           p.nc.diff + scale_x_continuous(limits = c(pl3.diff.lo, pl3.diff.hi)) +
                             scale_y_continuous(limits = c(pl3.diff.lo, pl3.diff.hi)),
                           p.none.diff + scale_x_continuous(limits = c(pl4.diff.lo, pl4.diff.hi)) +
                             scale_y_continuous(limits = c(pl4.diff.lo, pl4.diff.hi)),
                           nrow = 2, ncol = 2)


##### 1 family per frailty ######

load("OE_Sim_1Fam.RData")
## O/E ratios
p.cnc.1fam <- ggplot(res.oe.cnc.1fam, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.01, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.03, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "Frailty on all") +
  scale_color_continuous(name = "Frailty") +
  theme(plot.title = element_text(hjust = 0.5))

p.c.1fam <- ggplot(res.oe.c.1fam, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.01, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.03, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "Frailty on carriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.nc.1fam <- ggplot(res.oe.nc.1fam, aes(OE.C, OE.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.01, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.03, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "Frailty on noncarriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.none.1fam <- ggplot(res.oe.none.1fam, aes(OE.C, OE.NC)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = OE.NC.lo, ymax = OE.NC.hi), width = 0.01, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OE.C.lo, xmax = OE.C.hi), height = 0.03, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "O/E, carriers",
       y = "O/E, noncarriers",
       title = "No frailty") +
  theme(plot.title = element_text(hjust = 0.5))


## O-E differences
p.cnc.diff.1fam <- ggplot(res.oe.cnc.1fam, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.002, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.001, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O-E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "Frailty on all") +
  scale_color_continuous(name = "Frailty") +
  theme(plot.title = element_text(hjust = 0.5))

p.c.diff.1fam <- ggplot(res.oe.c.1fam, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.002, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.0008, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O/E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "Frailty on carriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.nc.diff.1fam <- ggplot(res.oe.nc.1fam, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(aes(color = W), size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.003, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.001, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O/E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "Frailty on noncarriers") +
  theme(plot.title = element_text(hjust = 0.5))

p.none.diff.1fam <- ggplot(res.oe.none.1fam, aes(OEdiff.C, OEdiff.NC)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = OEdiff.NC.lo, ymax = OEdiff.NC.hi), width = 0.002, alpha = 0.5) +
  geom_errorbarh(aes(xmin = OEdiff.C.lo, xmax = OEdiff.C.hi), height = 0.002, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "(O/E)/n, carriers",
       y = "(O/E)/n, noncarriers",
       title = "No frailty") +
  theme(plot.title = element_text(hjust = 0.5))


grid_arrange_shared_legend(p.cnc.1fam, p.c.1fam, p.nc.1fam,
                           p.none.1fam, nrow = 2, ncol = 2)
grid_arrange_shared_legend(p.cnc.diff.1fam, p.c.diff.1fam,
                           p.nc.diff.1fam, p.none.diff.1fam, nrow = 2, ncol = 2)

## square plots, same scale
pl1.lo.1fam <- min(c(res.oe.cnc.1fam$OE.C.lo, res.oe.cnc.1fam$OE.NC.lo))
pl1.hi.1fam <- max(c(res.oe.cnc.1fam$OE.C.hi, res.oe.cnc.1fam$OE.NC.hi))
pl2.lo.1fam <- min(c(res.oe.c.1fam$OE.C.lo, res.oe.c.1fam$OE.NC.lo))
pl2.hi.1fam <- max(c(res.oe.c.1fam$OE.C.hi, res.oe.c.1fam$OE.NC.hi))
pl3.lo.1fam <- min(c(res.oe.nc.1fam$OE.C.lo, res.oe.nc.1fam$OE.NC.lo))
pl3.hi.1fam <- max(c(res.oe.nc.1fam$OE.C.hi, res.oe.nc.1fam$OE.NC.hi))
pl4.lo.1fam <- min(c(res.oe.none.1fam$OE.C.lo, res.oe.none.1fam$OE.NC.lo))
pl4.hi.1fam <- max(c(res.oe.none.1fam$OE.C.hi, res.oe.none.1fam$OE.NC.hi))

grid_arrange_shared_legend(p.cnc.1fam + scale_x_continuous(limits = c(pl1.lo.1fam, pl1.hi.1fam)) +
                             scale_y_continuous(limits = c(pl1.lo.1fam, pl1.hi.1fam)),
                           p.c.1fam + scale_x_continuous(limits = c(pl2.lo.1fam, pl2.hi.1fam)) +
                             scale_y_continuous(limits = c(pl2.lo.1fam, pl2.hi.1fam)),
                           p.nc.1fam + scale_x_continuous(limits = c(pl3.lo.1fam, pl3.hi.1fam)) +
                             scale_y_continuous(limits = c(pl3.lo.1fam, pl3.hi.1fam)),
                           p.none.1fam + scale_x_continuous(limits = c(pl4.lo.1fam, pl4.hi.1fam)) +
                             scale_y_continuous(limits = c(pl4.lo.1fam, pl4.hi.1fam)),
                           nrow = 2, ncol = 2)


pl1.diff.lo.1fam <- min(c(res.oe.cnc.1fam$OEdiff.C.lo, res.oe.cnc.1fam$OEdiff.NC.lo))
pl1.diff.hi.1fam <- max(c(res.oe.cnc.1fam$OEdiff.C.hi, res.oe.cnc.1fam$OEdiff.NC.hi))
pl2.diff.lo.1fam <- min(c(res.oe.c.1fam$OEdiff.C.lo, res.oe.c.1fam$OEdiff.NC.lo))
pl2.diff.hi.1fam <- max(c(res.oe.c.1fam$OEdiff.C.hi, res.oe.c.1fam$OEdiff.NC.hi))
pl3.diff.lo.1fam <- min(c(res.oe.nc.1fam$OEdiff.C.lo, res.oe.nc.1fam$OEdiff.NC.lo))
pl3.diff.hi.1fam <- max(c(res.oe.nc.1fam$OEdiff.C.hi, res.oe.nc.1fam$OEdiff.NC.hi))
pl4.diff.lo.1fam <- min(c(res.oe.none.1fam$OEdiff.C.lo, res.oe.none.1fam$OEdiff.NC.lo))
pl4.diff.hi.1fam <- max(c(res.oe.none.1fam$OEdiff.C.hi, res.oe.none.1fam$OEdiff.NC.hi))

grid_arrange_shared_legend(p.cnc.diff.1fam + scale_x_continuous(limits = c(pl1.diff.lo.1fam, pl1.diff.hi.1fam)) +
                             scale_y_continuous(limits = c(pl1.diff.lo.1fam, pl1.diff.hi.1fam)),
                           p.c.diff.1fam + scale_x_continuous(limits = c(pl2.diff.lo.1fam, pl2.diff.hi.1fam)) +
                             scale_y_continuous(limits = c(pl2.diff.lo.1fam, pl2.diff.hi.1fam)),
                           p.nc.diff.1fam + scale_x_continuous(limits = c(pl3.diff.lo.1fam, pl3.diff.hi.1fam)) +
                             scale_y_continuous(limits = c(pl3.diff.lo.1fam, pl3.diff.hi.1fam)),
                           p.none.diff.1fam + scale_x_continuous(limits = c(pl4.diff.lo.1fam, pl4.diff.hi.1fam)) +
                             scale_y_continuous(limits = c(pl4.diff.lo.1fam, pl4.diff.hi.1fam)),
                           nrow = 2, ncol = 2)


