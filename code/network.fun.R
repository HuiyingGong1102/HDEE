get_interaction <- function(data, col, reduction = FALSE ){
  
  n <- nrow(data)
  clean_data <- data
  gene_list <- list()
  m <- clean_data[,col]
  M <- clean_data[,-col]
  x_matrix <- M
  x_matrix <- as.matrix(x_matrix)
  if (reduction!=FALSE) {
    vec <- sapply(1:length(M[1,]),function(c)cor(m,M[,c]))
    x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
    x_matrix <- as.matrix(x_matrix)
  }
  name <- colnames(clean_data)
  ridge1_cv <- cv.glmnet(x = x_matrix, y = m,alpha = 0)
  best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
  
  fit_res <- cv.glmnet(x = x_matrix, y = m,alpha = 1,
                       penalty.factor = 1/best_ridge_coef,
                       keep = TRUE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
  
  gene_list_one <- list()
  gene_list_one[[1]] <- name[col]
  gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
  gene_list_one[[3]] <- best_alasso_coef1@x[-1]
  gene_list[[col]] <- gene_list_one
  
  return(gene_list_one)
}


#' @title generate a growth curve(ind or dep effect curve)
#' @importFrom orthopolynom legendre.polynomials polynomial.derivatives scaleX
#' @param pars vector of LOP pars
#' @param effect scalar of observed genetic or any other data of a time point
#' @param times vector of time point
#' @param order the order of LOP
#' @param y0 scalar of initial y
#' @return vector of generated curve
get_effect <- function(pars,effect,times,order,y0){
  #Legendre polynomials
  LOP <-  legendre.polynomials(order, normalized=F)
  #derivatives of LOP
  d_LOP_fit <-  sapply(1:length(pars),function(c)
    pars[c]*polynomial.derivatives(LOP)[[c+1]])
  s_time <- scaleX(times,u=-1,v=1)
  h <- diff(s_time) #per step h
  f <- function(x){dy=do.call(sum,polynomial.values(polynomials=d_LOP_fit,x=x));dy}
  dy_fit <- c(y0)
  for (j in 1:(length(times)-1)) {
    k1 <- f(s_time[j])*effect[j]
    k2 <- f(s_time[j]+h[j]/2)*effect[j]
    k3 <- f(s_time[j]+h[j]/2)*effect[j]
    k4 <- f(s_time[j]+h[j])*effect[j]
    y <- dy_fit[j]+h[j]/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
    dy_fit <- c(dy_fit,y)
  }
  return(dy_fit)
}



#' @title calculate least-square for observed and fitted data
#' @param pars matrix of LOP parameters for ind and dep growth curve
#' @param ind the independent growth curve id
#' @param dep the dependent growth curve id
#' @param times vector of time point
#' @param data dataframe of observed data
#' @param order scalar of LOP order
#' @param effect matrix of observed data
#' @return scalar of least-square error
ode_optimize <- function(pars,ind,dep,times,data,order,effect){
  ind_pars <- matrix(pars,ncol=order)[1,]
  dep_pars <- matrix(pars,ncol=order)[-1,]
  inital_value <- data[,ind][1]
  ind_effect <- get_effect(ind_pars,data[,ind],times,order,inital_value)
  if ( is.null(nrow(dep_pars)) ) {
    dep_effect <- get_effect(dep_pars,data[,dep],times,order,0)
    y <- ind_effect+dep_effect
  }else{
    dep_effect <- sapply(1:length(dep), function(c)
      get_effect(dep_pars[c,],data[,dep[c]],times,order,0))
    y <- ind_effect+rowSums(dep_effect)
  }
  ssr <- sum((data[,ind]-y)^2)
  #add penalty
  alpha=0.01
  if(length(which(ind_effect<0))>0){
    # return(sum((data[,ind]-y)^2+sum((ind_effect[ind_effect<0])^2)))
    # return(sum((data[,ind]-y)^2+alpha*(sum(ind_pars^2)+sum(dep_pars^2))))
    return(sum((data[,ind]-y)^2+alpha*(sum(abs(ind_pars))+sum(abs(dep_pars)))))
  }else{
    return(ssr)
  }
}


#' @title estimate parameters for LOP
#' @importFrom stats optim
#' @param effect matrix of observed data
#' @param relationship list contain the result of lasso-based variable election
#' @param times vector of time point
#' @param order scalar of LOP order
#' @return matrix of estimated LOP parameters
get_value <- function(effect,relationship,times,order){
  #input
  ind <- relationship[[1]]
  dep <- relationship[[2]]
  ind_no <- as.numeric(which(colnames(effect)==ind))
  dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effect)==dep[c])))
  init_pars <- rep(0.001,(length(ind_no)+length(dep_no))*order)
  result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,
                  times=times,data=effect,order=order,
                  method = "Nelder-Mead", control=list(maxit=50000,trace=T))
  par_after <- matrix(result$par,length(ind)+length(dep),order)
  return(par_after)
}

#' @title calculation all LOP parameters in parallel mode
#' @import parallel
#' @import pbapply
#' @param data matrix of observed data
#' @param times vector of time point
#' @param order scalar of LOP order
#' @param reduction use n/log(n) dimension reduction
#' @param parallel use parallel computation or not
#' @return list contain variable selection results and LOP parameters for every row
#' @export
get_ode_par <- function(data, times, order, reduction, parallel = TRUE) {
  cat('Start variable selection',sep="\n")
  relationship <- pblapply(1:nrow(data),function(c)
    get_interaction(data, c, reduction = reduction))
  cat('Finish variable selection',sep="\n")
  
  cat('Start ODE solving',sep="\n")
  if (parallel == TRUE) {
    core.number <- detectCores()
    cl <- makeCluster(getOption("cl.cores", core.number))
    clusterEvalQ(cl, {require(orthopolynom)})
    clusterExport(cl, c("data","times","order","reduction","relationship",
                        "get_value","ode_optimize","get_effect","LOP_rk4"),envir=environment())
    lop_par <- pblapply(1:nrow(data),function(c)
      get_value(t(data),relationship[[c]],times,order),cl=cl)
    stopCluster(cl)
  }
  else{
    lop_par <- pblapply(1:nrow(data),function(c)
      get_value(t(data),relationship[[c]],times,order))
  }
  cat('Finish ODE solving',sep="\n")
  return_obj <- list(lop_par = lop_par,
                     relationship = relationship,
                     dataset = data,
                     times = times,
                     order = order)
  return(return_obj)
}

#' @title helper function to convert ODE result
#' @param relationship list contain the result of lasso-based variable election
#' @param par vector conatin LOP parameters
#' @param effect matrix of observed data
#' @param times vector of time point
#' @param order scalar of LOP order
#' @return list contain bunch of useful output result for a row
get_ode_output <- function(relationship, par, effect, times, order){
  #effect <- t(effect)
  output <- list()
  output[[1]] <- relationship[[1]]
  output[[2]] <- relationship[[2]]
  output[[3]] <- par[1,]
  output[[4]] <- par[2:nrow(par),]
  ind_no <- as.numeric(which(colnames(effect)==output[[1]]))
  dep_no <- as.numeric(sapply(1:length(output[[2]]),
                              function(c) which(colnames(effect)==output[[2]][c])))
  inital_value <- effect[,ind_no][1]
  ind_effect <- get_effect(as.numeric(output[[3]]),effect[,ind_no],times,order,inital_value)
  if (length(dep_no)==1) {
    dep_effect <- get_effect(as.numeric(output[[4]]),effect[,dep_no],times,order,0)
  }else{
    dep_effect <- sapply(1:length(dep_no), function(c)
      get_effect(as.numeric(output[[4]][c,]),effect[,dep_no[c]],times,order,0))
    colnames(dep_effect) <- dep_no
  }
  all_effect <- cbind(ind_effect,dep_effect)
  effect_mean <- apply(all_effect,2,mean)
  output[[5]] <- effect_mean
  output[[6]] <- all_effect
  return(output)
}

#' @title convert all ODE result into a list(use when multiply network needed)
#' @param input list contain the ODE solving result from get_ode_par function
#' @return list of output for all data
#' @export
get_all_net <- function(input){
  effect <- input$dataset
  n = ncol(effect)
  par = input$lop_par
  order = input$order
  times = input$times
  relationship = input$relationship
  
  net <- lapply(1:n,function(c)
    get_ode_output(relationship = relationship[[c]],
                   par = par[[c]],
                   effect = effect,
                   times = times,
                   order = order))
  return(net)
}

#' @title helper function of plot decomposition_curve
#' @param input1 list object from get_ode_par function
#' @param input2 list object from get_all_net function
#' @param i scalar of row number
#' @return datafrmae for plot
get_curve_data <- function(input1, input2, i){
  times <- input1$times
  cluster_mean <- input1$dataset
  d1 <- data.frame(input2[[i]][[6]],check.names = F)
  colnames(d1)[-1] <- input2[[i]][[2]]
  d1$sum <- rowSums(d1)
  d1$origin <- as.numeric(cluster_mean[,i])
  d1$x <- times
  return(d1)
}


#' @title plot decomposition_curve
#' @param input1 list object from get_ode_par function
#' @param input2 list object from get_all_net function
#' @param i scalar of i row number to plot
#' @return plot
#' @export
get_decomposition_plot <- function(input1, input2, i){
  col <- c('#ff7171','#0dceda','#9ede73')
  plot_data <- get_curve_data(input1, input2, i)
  ind_name <- input2[[i]][[1]]
  
  d <- plot_data
  p1 <- ggplot(d,aes(x=x))+geom_line(d,mapping=aes(x=x,y=origin),color=col[2],size=1.5)+
    geom_line(d,mapping=aes(x=x,y=ind_effect),color=col[1],size=1.4)
  for (j in 1:(ncol(d)-4)) {
    p1 <- p1+geom_line(aes_string(y=d[,1+j]),color=col[3],size=1.2,alpha=1)+
      annotate('text',x=max(as.numeric(d$x))+0.03,y=d[nrow(d),(j+1)],
               label=colnames(d)[j+1],size=3)
    
  }
  p1 <- p1+ theme(panel.grid = element_blank(),
                  panel.background = element_rect(color = 'black', fill = 'transparent'),
                  legend.title = element_blank()) +
    xlab("Expression Index") + ylab('Effect')+
    theme(plot.margin = unit(c(0,0,0,0),"lines"))+
    geom_hline(yintercept = 0, size = 0.6) +
    annotate('text',x=mean(as.numeric(d$x)),y=round(max(d[,-ncol(d)])),label=ind_name,size=5)
  
  return(p1)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
