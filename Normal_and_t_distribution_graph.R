curve(dnorm(x,0,1),xlim=c(-3,3),main='Normal Density')

cord.x <- c(-3)
cord.y <- c(0)

cord.x <- c(cord.x,-3) 
cord.y <- c(cord.y,dnorm(-3))

cord.x <- c(cord.x,-2,-2)
cord.y <- c(cord.y,dnorm(-2),0)

polygon(cord.x,cord.y,col='skyblue')

# Second attempt

cord.x <- c(-3,seq(-3,-2,0.01),-2) 
cord.y <- c(0,dnorm(seq(-3,-2,0.01)),0) 
curve(dnorm(x,0,1),xlim=c(-3,3),main='Standard Normal') 
polygon(cord.x,cord.y,col='skyblue')