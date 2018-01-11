#This function calculates the rich-club coefficient  and the complementary rich-club coefficient 
#of an unweighted network.
#The rich-club coefficient measures the density of connections among the nodes with degree higher than a certain
#threshold k (i.e. the density of the so called 'rich-club connections').  
#Ref:
#Zhou S., Mondragon R.J. (2004) The rich-club phenomenon in the internet topology. IEEE Comm Lett, 8:180-182.

#The complementary rich-club coefficient measures the density of connections among the nodes with degree higher
#than a certain threshold k and the nodes with degree lower/equal than k (i.e. the density of the so called
#'feeder connections'). The complementary rich-club coefficient is useful in order to quantify the extent 
#to which the core is connected to the rest of the network (i.e. it provides further evidence in the case of 
#rich-club ordering). The complementary rich-club coefficient can be normalized just like the rich-club coeff.
#Ref:
#M. Cinelli, G. Ferraro, and A. Iovanella, “Rich-club ordering and the dyadic effect: Two interrelated phenomena”,
#Physica A: Statistical Mechanics and its Applications, vol. 490, pp. 808–818,2018.

library(igraph)

rich.club.coeff.new <- function (g, k = 1) 
{
  if (!is.igraph(g)) {
    stop(sprintf("%s is not a graph object", deparse(substitute(g))))
  }
  if ("degree" %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  }
  else {
    degs <- degree(g)
  }
  E <- ecount(g) 
  Nv <- vcount(g)
  Nk <- sum(degs > k) # n1
  N_mineq_k <- sum(degs <= k) # n0
  
  if (Nk == 0) {                         #
    
    poor.club.nodes <- order(degs)[1 : (Nv - Nk)] #
    poor.club.graph <- induced.subgraph(g, poor.club.nodes) #
    Ek_poor <- ecount(poor.club.graph) # m00
    Ek_segn <- E - Ek_poor #
    phi_segn <- (Ek_segn) / (Nk * N_mineq_k)
    
    return(list(phi = NaN,# graph = make_empty_graph(),
                Nk = 0, 
                Ek = 0, phi_segn = phi_segn , N_mineq_k = N_mineq_k, Ek_segn = Ek_segn,
                Ek_poor = Ek_poor, UBm11 = NaN, phi_new = NaN,
                UBm10 = NaN, phi_segn_new = NaN))
  }
  else {
    rich.club.nodes <- order(degs)[(Nv - Nk + 1):Nv] #Consider the Nk nodes

    rich.club.graph <- induced.subgraph(g, rich.club.nodes)
    Ek <- ecount(rich.club.graph)
    
    poor.club.nodes <- order(degs)[1 : (Nv - Nk)] #
    poor.club.graph <- induced.subgraph(g, poor.club.nodes) #
    Ek_poor <- ecount(poor.club.graph) #
    
    Ek_segn <- E - Ek - Ek_poor #
    
    phi <- graph.density(rich.club.graph)
    phi_segn <- (Ek_segn) / (Nk * N_mineq_k)
    
    ds <- sort(degs, decreasing = T)
    UBm11 <- m11.max.new(vcount(rich.club.graph),ecount(g),ds)
    phi_new <- ecount(rich.club.graph) / UBm11
    
    
    UBm10 <- m10.max.new(vcount(rich.club.graph),vcount(g),ecount(g), ds)
    phi_segn_new <- (Ek_segn) / UBm10
    
    return(list(phi = phi,# graph = rich.club.graph, 
                Nk = Nk, 
                Ek = Ek, phi_segn = phi_segn , N_mineq_k = N_mineq_k, Ek_segn = Ek_segn,
                Ek_poor = Ek_poor, UBm11 = UBm11, phi_new = phi_new,
                UBm10 = UBm10, phi_segn_new = phi_segn_new ))
  }
}
