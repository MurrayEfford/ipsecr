###############################################################################
## package 'ipsecr'
## plot3D.IP.R
## 2022-08-25
###############################################################################

plot3D.IP <- function (object, box = 1, oldplot = NULL, plotcentre = TRUE, 
    plotfinal = FALSE, zkludge = -0.2) {
    
    if (!inherits(object, 'ipsecr') || is.null(object$simulations)) {
        stop("object must be of class 'ipsecr'")
    }
    if (is.null(object$simulations)) {
        stop("object does not include saved simulations; use details=list(keep.sim = TRUE)")
    }
    if (!requireNamespace('plot3D')) {
        stop ("plot3D.IP requires package plot3D")
    }
    if (object$details$factorial != "full") {
        stop ("plot3D.IP is not ready for fractional factorial designs")
    }
    pred <- predict(object)
    if (nrow(pred) != 3) {
        stop("plot3D.IP assumes 3 parameters, one session")
    }
    
    #------------------------
    # set up
    parnames <- rownames(pred)
    final <- coef(object)[,'beta']
    proxynames <- rownames(object$variance.bootstrap)
    box1 <- as.data.frame(apply(object$parameters[[box]], 2, range))
    vert1 <- expand.grid(box1)
    cent1 <- apply(box1, 2, mean)
    
    sim1 <- object$simulations[[box]]
    vertex <- rep(c(1:8, rep(9, object$details$centre)), length.out=nrow(sim1))
    sequ <- c(1,2,4,3,1,5,6,8,7,5,NA,2,6,NA,3,7,NA,4,8)
    obs <- object$variance.bootstrap$target
    
    if (is.null(oldplot)) {
        pr <- box1
        pr[,1] <- pr[,1]+c(-0.15,0.1)
        pr[,2] <- pr[,2]+c(0,0.3)
        pr[,3] <- pr[,3]+c(-0.15,0.15)
        sr <- apply(sim1,2,range)
        sr <- sweep(sr, MARGIN=1, FUN='*', STATS = c(0.98,1.02))
    } 
    else {
        pr <- oldplot$pr
        sr <- oldplot$sr
    }
    
    #---------------------------------------
    # experimental design in parameter space
    pmatparm <- plot3D::box3D(
        r = 3, theta = 55, phi = 20, 
        x0 = box1[1,1], y0 = box1[1,2], z0 = box1[1,3],
        x1 = box1[2,1], y1 = box1[2,2], z1 = box1[2,3],
        col = plot3D::gg.col(8, alpha = 0.05), lty = 1, box = T,
        lwd = c(1, 4), main = "", border = 'black',
        xlim = pr[,1], ylim = pr[,2], zlim = pr[,3],
        xlab = parnames[1], ylab = parnames[2], zlab = parnames[3]
    )
    
    # override projection if matching previous plot
    if (!is.null(oldplot)) {
        pmatparm <- oldplot$pmatparm
    }

    # add points
    XY1 <- plot3D::trans3D (x=vert1[,1], y=vert1[,2], z=vert1[,3], pmat = pmatparm)
    points(XY1, bg = plot3D::jet.col(8), cex=1.5, pch=21)  # 'orange'
    if (plotcentre) {
        XY <- plot3D::trans3D(x = cent1[1], y=cent1[2], z=cent1[3], pmat = pmatparm) 
        points(XY, cex=1.6, pch=23, bg = 'blue')
    }
    if (plotfinal) {
        XYfinal   <- plot3D::trans3D(x = final[1], y = final[2], z = final[3], pmat = pmatparm) 
        # ad hoc +zkludge (cannot find bug)
        basefinal <- plot3D::trans3D(x = final[1], y = final[2], z = pr[1,3]+zkludge, pmat = pmatparm)
        lines(c(XYfinal$x, basefinal$x), c(XYfinal$y, basefinal$y))
        points(XYfinal, cex=1.6, pch=22, bg = 'white')
    }

    #-----------------------------------------
    # 3D scatterplot of simulated proxy values
    pmatsim <- plot3D::scatter3D(
        x = sim1[vertex<9,1], y = sim1[vertex<9,2],z = sim1[vertex<9,3],
        r = 3, theta = 55, phi = 20, colvar = vertex[vertex<9],
        xlab = proxynames[1], ylab=proxynames[2], zlab=proxynames[3],
        xlim = sr[,1], ylim = sr[,2], zlim = sr[,3],
        colkey = FALSE)
    
    # override projection if matching previous plot
    if (!is.null(oldplot)) {
        pmatsim <- oldplot$pmatsim
    }
    
    #------------------------------------------
    # add mean at each design point, plus frame
    mn <- apply(sim1,2,tapply, vertex, mean)
    xy <- as.data.frame(plot3D::trans3D(x=mn[,1], y=mn[,2], z=mn[,3], pmat = pmatsim))
    lines(xy[sequ,1], xy[sequ,2], lwd=1.2)
    points(xy[1:8,1], xy[1:8,2], pch=21, bg = plot3D::jet.col(8), cex=1.5)
    
    #-------------------------------
    # observed value of proxy vector
    obsxy <- plot3D::trans3D(x=obs[1], y=obs[2], z=obs[3], pmat = pmatsim)
    baseobs <- plot3D::trans3D(x=obs[1], y=obs[2],z=sr[1,3], pmat = pmatsim)
    lines(c(obsxy$x, baseobs$x), c(obsxy$y, baseobs$y))
    points(obsxy, pch = 22, cex = 1.6, bg = 'yellow')
    
    # silently return projection matrices and axis limits
    invisible(list(pmatparm = pmatparm, pmatsim = pmatsim, pr = pr, sr = sr))
}

# ipsecrdemo
# par(mfrow=c(2,2), oma = c(1,1,3,1))
# oldplot <- plot3D.IP(ipsecrdemo, box = 1)
# plot3D.IP(ipsecrdemo, box = 2, oldplot, plotfinal = TRUE, zkludge=-0.1)
# mtext(outer = TRUE, side = 3, line = 0.5, adj = c(0.2,0.8), cex = 1.1, 
#     c('Parameter space', 'Proxy space'))

