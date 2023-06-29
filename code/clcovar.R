AR1.get_mat <- function(par0,light, traits=1, options=list())
{
  par<-par0;
  if (class(par0)=="list")
    par <- unlist(par0);
  
  t_len <- length(light);
  
  Ar.1 <- array(0, dim=c(t_len*traits,t_len*traits));
  for (i0 in 1:traits)
    for (i1 in 1:traits)
    {
      if (i0==i1)
        for (k0 in 1:t_len)
          for (k1 in 1:t_len)
          {
            Ar.1[(i0-1)*t_len+k0,(i1-1)*t_len+k1] <- par[i0*2]^2 * par[i0*2-1]^abs( k0 - k1 );
          }
    }
  
  return(Ar.1);
}

AR1.get_inv_mat <- function(par,light, traits=1, options=list())
{
  Ar.1 <- AR1.get_mat(par,light, traits, options)
  return(solve(Ar.1));
}

AR1.get_mat_det <- function(par,light, traits=1, options=list())
{
  Ar.1 <- AR1.get_mat(par,light, traits, options)
  return(det(Ar.1));
}