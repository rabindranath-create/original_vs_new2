#install.packages("spatial")
#install.packages("spatstat")
#install.packages("igraph")

library(spatial)
library(spatstat)
library(igraph)


#########################################################
#####               Basic Functions                 #####
#########################################################
# 1. discretize graph
Graph_Discretized<-function(x,y){
  n=(x+1)*(y+1)
  temp<-make_empty_graph(n, directed = FALSE)
  #create graph with n verticles, no edges
  
  G<-set_vertex_attr(temp, "name", value = as.character(1:n))
  
  # adds or modifies a vertex attribute called "name" for the graph temp.
  
  
  A<-as_adjacency_matrix(G, sparse = FALSE)
  #adjacency matrix 
  
  #dim(A)=n by n
  #assume that x,y greater than 1.
  for(j in 1:(y-1)){
    for(i in 1:(x-1)){
      
      A[1+i+j*(x+1),i+j*(x+1)]<-1 #left
      A[1+i+j*(x+1),2+i+j*(x+1)]<-1 #right
      A[1+i+j*(x+1),1+i+(j-1)*(x+1)]<-1 #down
      #i for horizontal, j for vertical
      A[1+i+j*(x+1),1+i+(j+1)*(x+1)]<-1 #up
      
      A[1+i+j*(x+1),i+(j+1)*(x+1)]<-sqrt(2) #left up corner
      A[1+i+j*(x+1),2+i+(j+1)*(x+1)]<-sqrt(2) #right up corner
      A[1+i+j*(x+1),i+(j-1)*(x+1)]<-sqrt(2) #left down corner
      A[1+i+j*(x+1),2+i+(j-1)*(x+1)]<-sqrt(2) #right down corner
      
    }#end of inner loop
  }#end of outer loop
  
  #j=0
  for(i in 1:(x-1)){
    A[1+i,i]<-1; A[1+i,2+i]<-1; A[1+i,1+i+(x+1)]<-1;
    A[1+i,i+(x+1)]<-sqrt(2); A[1+i,2+i+(x+1)]<-sqrt(2); 
  }
  #j=y
  for(i in 1:(x-1)){
    A[1+i+y*(x+1),i+y*(x+1)]<-1; A[1+i+y*(x+1),2+i+y*(x+1)]<-1; A[1+i+y*(x+1),1+i+(y-1)*(x+1)]<-1;
    A[1+i+y*(x+1),i+(y-1)*(x+1)]<-sqrt(2); A[1+i+y*(x+1),2+i+(y-1)*(x+1)]<-sqrt(2); 
  }
  #i=0
  for(j in 1:(y-1)){
    A[1+j*(x+1),1+(j-1)*(x+1)]<-1; A[1+j*(x+1),2+j*(x+1)]<-1; A[1+j*(x+1),1+(j+1)*(x+1)]<-1;
    A[1+j*(x+1),2+(j-1)*(x+1)]<-sqrt(2); A[1+j*(x+1),2+(j+1)*(x+1)]<-sqrt(2); 
  }
  #i=x
  for(j in 1:(y-1)){
    A[1+x+j*(x+1),1+x+(j-1)*(x+1)]<-1; A[1+x+j*(x+1),1+x+(j+1)*(x+1)]<-1; A[1+x+j*(x+1),x+j*(x+1)]<-1;
    A[1+x+j*(x+1),x+(j-1)*(x+1)]<-sqrt(2); A[1+x+j*(x+1),x+(j+1)*(x+1)]<-sqrt(2); 
  }
  #4 Corners
  A[1,2]<-1; A[1,2+x]<-1; A[1,3+x]<-sqrt(2);
  A[1+y*(x+1),1+(y-1)*(x+1)]<-1; A[1+y*(x+1),2+y*(x+1)]<-1; A[1+y*(x+1),2+(y-1)*(x+1)]<-sqrt(2);
  A[1+x,x]<-1; A[1+x,2*(1+x)]<-1; A[1+x,2*(1+x)-1]<-sqrt(2);
  A[(1+x)*(1+y),(1+x)*y]<-1; A[(1+x)*(1+y),(1+x)*(1+y)-1]<-1; A[(1+x)*(1+y),(1+x)*y-1]<-sqrt(2);
  
  #G<-graph.adjacency(A,mode=c("undirected"),weighted=TRUE)
  G<-graph_from_adjacency_matrix( A,mode=c("undirected"),weighted=TRUE)
  return(G)
}



##########
#Examples
##########

#G<- Graph_Discretized(2,2)
#dev.new()
#plot(G, edge.label = round(E(G)$weight, 2))


#summary(G)
#class(G)
#id <- tkplot(G)
#coords <- tkplot.getcoords(id)
#plot(G, layout = coords)

#nodes <- data.frame(id = V(G)$name)
#edges <- as_data_frame(G, what = "edges")

#visNetwork(nodes, edges) %>% visEdges(arrows = "to")

##########
#Examples
##########





# 3. get index for certain coordinates
Index_Coordinates<-function(m,x,y){
  # m is the coordinates of points of interests
  temp<-1+m[1]+m[2]*(x+1)
  return(temp)
}


# 4. generate coordinates for lattice
Lattice_Vertices <- function(x,y){
  LatticeCoordinates<-matrix(nrow=(x+1)*(y+1),ncol=2)
  temp<-rep(0:x,y+1)	
  # 0 1 ... x repeat y+1 times 
  LatticeCoordinates[,1]<- temp
  temp<-rep(0,0)
  for(i in 0:y){
    temp<-c(temp,rep(i,x+1))
  }
  LatticeCoordinates[,2]<-temp
  temp<-rep(0,0)
  return(LatticeCoordinates)
}
#index -> coordinate 


# 5. distance between two points
Dist_Euclidean <- function(point1,point2){
  result <- as.numeric(sqrt((point1[1]-point2[1])^2+(point1[2]-point2[2])^2))
  return(result)
}

# 6. get all lattice coordinates within a circle
Lattice_inCircle<-function(c,r){
  # input is c=(cx,cy) x,y coordinate of a circle. r-radius
  cx<-c[1]; cy<-c[2]
  cx1<-floor(cx); cy1<-floor(cy)
  r1<-ceiling(r)+1 # no need to add 2?
  temp<-c()
  for(j in -r1:r1){
    for(i in -r1:r1){
      #Center(cx,cy) and latticepoint(cx1+i,cy1+j)
      if((cx-(cx1+i))^2+(cy-(cy1+j))^2 <= r^2){
        #cx1+i, cy1+j within circle?
        temp<-c(temp,c(cx1+i,cy1+j))	
      }
    }#inner loop
  }#outer loop
  if(length(temp)>=2)
    lattices<-matrix(temp,ncol=2,byrow=TRUE)
  else
    return(0)
  return(lattices)
}
#Lattice_inCircle(c(2.5, 2.5), 1)










# 2. intersected edges
Intersect_Obs <- function(c,r,x,y){
  # input is c=(cx,cy) x,y coordinate of a circle. r-radius
  cx<-c[1]; cy<-c[2]
  cx1<-floor(cx); cy1<-floor(cy)
  r1<-ceiling(r)+2
  coor_info <- Lattice_Vertices(x,y)
  
  # points outside the boundary of circle
  X_temp <- matrix(c(rep(-r1:r1,times=length(seq(-r1,r1,by=1))),
                     rep(-r1:r1,each=length(seq(-r1,r1,by=1)))),ncol=2)
  Intersect_temp1 <- function(vector_ij){
    x_i<-vector_ij[1];x_j<-vector_ij[2]
    if(r^2< (cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 & (cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 <=(r+sqrt(2))^2){
      return(c(cx1+x_i,cy1+x_j))	
    }
  }
  temp<-as.numeric(unlist(apply(X_temp,1,Intersect_temp1)))  
  case2<-matrix(temp,ncol=2,byrow=TRUE)
  
  # points inside and on the boundary of circle
  Intersect_temp2 <- function(vector_ij){
    x_i<-vector_ij[1];x_j<-vector_ij[2]
    if((r-sqrt(2))^2<=(cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 &  (cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 <= r^2){
      return(c(cx1+x_i,cy1+x_j))	
    }
  }
  temp<-as.numeric(unlist(apply(X_temp,1,Intersect_temp2)))
  case3<-matrix(temp,ncol=2,byrow=TRUE)
  
  #I am going to form edgelist matrix intersecting the circle
  el<-matrix(0,ncol=2)
  n2=nrow(case2); n3=nrow(case3)
  
  #temp1<-matrix(0,ncol=4)
  
  if(n2!=0 & n3!=0){
    X_temp2 <- matrix(c(rep(1:n3,times=length(seq(1,n2,by=1))),
                        rep(1:n2,each=length(seq(1,n3,by=1)))),ncol=2)
    Intersect_temp3 <- function(vector_ij){
      x_i<-vector_ij[1];x_j<-vector_ij[2]
      if(Dist_Euclidean(case2[x_j,],case3[x_i,])==1 || Dist_Euclidean(case2[x_j,],case3[x_i,])==sqrt(2)){
        e1<-which(coor_info[,1]==case2[x_j,1]&coor_info[,2]==case2[x_j,2])
        e2<-which(coor_info[,1]==case3[x_i,1]&coor_info[,2]==case3[x_i,2])
        return(sort(c(e1,e2),decreasing=F))
      }
    }
    el_vector <- as.numeric(unlist(apply(X_temp2,1,Intersect_temp3)))
    el <- matrix(el_vector,ncol=2,byrow=T)
  }#if statement
  
  #for(i in 2:89){points(x=c(temp1[i,1],temp1[i,3]),y=c(temp1[i,2],temp1[i,4]),type="l",col="red")}
  return(el)
}

#matrix with colum 2 one index -> one index 







# Main algorithm - DT
# 1. update graph weights according to DT algorithm
Update_graph_intersect_DT<-function(g,x,y,circle_info,r){
  #read circle center x,y coordinate c-cost,p-prabability,True or False Obstacles
  #circles=read.csv("example1.csv",header=FALSE)
  n <- nrow(circle_info)
  elg <- as_data_frame(g,what="edges")
  colnames(elg) <- c("From","To","Cost")
  int_info <- matrix(0,ncol=nrow(circle_info),nrow=nrow(elg))
  for(i in 1:n){
    el<-Intersect_Obs(t(circle_info[i,1:2]),r,x,y)
    n1=nrow(el)
    dt <- Dist_Euclidean(as.numeric(circle_info[i,1:2]),c(50,1))
    for(k in 1:n1){el[k,]<-sort(el[k,],decreasing=FALSE)} #sort the elements
    
    for(j in 1:n1){
      index=which((elg[,1]==el[j,1] & elg[,2]==el[j,2]))
      elg[index,3] <- elg[index,3]+0.5*(circle_info[i,3]+(dt/(1-circle_info[i,4]))^(-log(1-circle_info[i,4])))
      int_info[index, i] <- 1
    }#inner loop
  }#outer loop
  
  updateg=graph_from_data_frame(elg,directed=FALSE)
  output <- list(G_info=updateg, Int_info=int_info)
  return(output)
}



DT_Alg <- function(obs_info){
  # generate obstacle info
  x <- 100; y <- 100; r <- 4.5
  # begin the loop to travel from s to t
  s <- 10151
  t <- 152
  # create graph - based on W
  vertice_list <- Lattice_Vertices(x,y)
  G_original <- Graph_Discretized(x,y)
  output_Ginfo <- Update_graph_intersect_DT(G_original, x, y, obs_info, r)
  G_ed <- output_Ginfo$G_info
  Int_info <- output_Ginfo$Int_info
  df_edge_ed <- as_data_frame(G_ed, what="edges")
  # some record vectors
  length_total <- 0 # record Euclidean length
  cost_total <- 0 # record disambiguation cost 
  reach_t <- F
  path_record <- s
  D_record <- c()
  while(reach_t!=T){
    # implement shortest path algorithm
    output <- shortest_paths(G_ed,
                             which(vertex.attributes(G_ed)$name==as.character(s)),
                             which(vertex.attributes(G_ed)$name==as.character(t)),
                             weights=df_edge_ed$Cost, output="both",algorithm="dijkstra")
    V_list <- c(as.numeric(attributes(output$vpath[[1]])$names))
    # follow the path until the disambiguation state
    for (i in 2:length(V_list)){
      edge_ind_temp <- which(df_edge_ed$from==min(V_list[(i-1):i])&df_edge_ed$to==max(V_list[(i-1):i]))
      edge_length <- Dist_Euclidean(as.numeric(vertice_list[V_list[i-1],1:2]),as.numeric(vertice_list[V_list[i],1:2]))
      if (sum(Int_info[edge_ind_temp,])!=0){
        #cross obstacle? 
        D_state <- V_list[i-1]
        D_record <- c(D_record, D_state)
        break
      } else{
        #no cross, just go 
        length_total <- length_total+edge_length
        path_record <- c(path_record,V_list[i])
        D_state <- NULL
      }
    }
    if(is.null(D_state)){
      # if didn't run into any obstacle and reach target
      reach_t=T
      output_final <- list(Optimal_found=T,Length_total=length_total,Cost_total=cost_total,
                           Optimal_path=path_record, Disambiguate_state=D_record)
    } else{
      # run into obstacle
      # update start to current disambiguation state
      # subtract one disambiguation 
      reach_t=F
      s <- D_state
      # adjust graph
      # determine which obstacle
      obs_ind_temp <- which(Int_info[edge_ind_temp,]==1)
      if (length(obs_ind_temp)==1){
        # add cost of disambiguation
        cost_total <- cost_total+obs_info[obs_ind_temp,3]
        if(obs_info$status[obs_ind_temp]==1){
          # adjust based on true obstalce
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
        } else{
          dt <- Dist_Euclidean(as.numeric(obs_info[obs_ind_temp,1:2]),c(50,1))
          # adjust based on false obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp]==1),3]-
                                                                    0.5*(obs_info[obs_ind_temp,3]+(dt/(1-obs_info[obs_ind_temp,4]))^(-log(1-obs_info[obs_ind_temp,4]))))
          Int_info[which(Int_info[,obs_ind_temp]==1),obs_ind_temp] <- 0 
        }
      } else{
        dist_temp <- rep(0,length(obs_ind_temp))
        for(i in 1:length(obs_ind_temp)){
          dist_temp[i] <- Dist_Euclidean(as.numeric(vertice_list[D_state,1:2]),obs_info[obs_ind_temp[i],1:2])
        }
        obs_ind_temp2 <- obs_ind_temp[which.min(dist_temp)]
        # add cost of disambiguation
        cost_total <- cost_total + obs_info[obs_ind_temp2,3]
        if (obs_info$status[obs_ind_temp2]==1){
          # true obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
        } else{
          dt <- Dist_Euclidean(as.numeric(obs_info[obs_ind_temp2,1:2]),c(50,1))
          # false obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3]-
                                                                     0.5*(obs_info[obs_ind_temp2,3]+(dt/(1-obs_info[obs_ind_temp2,4]))^(-log(1-obs_info[obs_ind_temp2,4]))))
          Int_info[which(Int_info[,obs_ind_temp2]==1),obs_ind_temp2] <- 0
        }
      }
      G_ed <- graph_from_data_frame(df_edge_ed,directed = FALSE)
    }
  }
  return(output_final)
}



#########
#fixed DT 
########

Update_graph_intersect_DT_fixed<-function(g,x,y,circle_info,r, dt){
  #read circle center x,y coordinate c-cost,p-prabability,True or False Obstacles
  #circles=read.csv("example1.csv",header=FALSE)
  n <- nrow(circle_info)
  elg <- as_data_frame(g,what="edges")
  colnames(elg) <- c("From","To","Cost")
  int_info <- matrix(0,ncol=nrow(circle_info),nrow=nrow(elg))
  for(i in 1:n){
    el<-Intersect_Obs(t(circle_info[i,1:2]),r,x,y)
    n1=nrow(el)
    for(k in 1:n1){el[k,]<-sort(el[k,],decreasing=FALSE)} #sort the elements
    
    for(j in 1:n1){
      index=which((elg[,1]==el[j,1] & elg[,2]==el[j,2]))
      elg[index,3] <- elg[index,3]+0.5*(circle_info[i,3]+(dt/(1-circle_info[i,4]))^(-log(1-circle_info[i,4])))
      int_info[index, i] <- 1
    }#inner loop
  }#outer loop
  
  updateg=graph_from_data_frame(elg,directed=FALSE)
  output <- list(G_info=updateg, Int_info=int_info)
  return(output)
}



DT_fixed_Alg <- function(obs_info, dt){
  # generate obstacle info
  x <- 100; y <- 100; r <- 4.5
  # begin the loop to travel from s to t
  s <- 10151
  t <- 152
  # create graph - based on W
  vertice_list <- Lattice_Vertices(x,y)
  G_original <- Graph_Discretized(x,y)
  output_Ginfo <- Update_graph_intersect_DT_fixed(G_original, x, y, obs_info, r, dt)
  G_ed <- output_Ginfo$G_info
  Int_info <- output_Ginfo$Int_info
  df_edge_ed <- as_data_frame(G_ed, what="edges")
  # some record vectors
  length_total <- 0 # record Euclidean length
  cost_total <- 0 # record disambiguation cost 
  reach_t <- F
  path_record <- s
  D_record <- c()
  while(reach_t!=T){
    # implement shortest path algorithm
    output <- shortest_paths(G_ed,
                             which(vertex.attributes(G_ed)$name==as.character(s)),
                             which(vertex.attributes(G_ed)$name==as.character(t)),
                             weights=df_edge_ed$Cost, output="both",algorithm="dijkstra")
    V_list <- c(as.numeric(attributes(output$vpath[[1]])$names))
    # follow the path until the disambiguation state
    for (i in 2:length(V_list)){
      edge_ind_temp <- which(df_edge_ed$from==min(V_list[(i-1):i])&df_edge_ed$to==max(V_list[(i-1):i]))
      edge_length <- Dist_Euclidean(as.numeric(vertice_list[V_list[i-1],1:2]),as.numeric(vertice_list[V_list[i],1:2]))
      if (sum(Int_info[edge_ind_temp,])!=0){
        D_state <- V_list[i-1]
        D_record <- c(D_record, D_state)
        break
      } else{
        length_total <- length_total+edge_length
        path_record <- c(path_record,V_list[i])
        D_state <- NULL
      }
    }
    if(is.null(D_state)){
      # if didn't run into any obstacle and reach target
      reach_t=T
      output_final <- list(Optimal_found=T,Length_total=length_total,Cost_total=cost_total,
                           Optimal_path=path_record, Disambiguate_state=D_record)
    } else{
      # run into obstacle
      # update start to current disambiguation state
      # subtract one disambiguation 
      reach_t=F
      s <- D_state
      # adjust graph
      # determine which obstacle
      obs_ind_temp <- which(Int_info[edge_ind_temp,]==1)
      if (length(obs_ind_temp)==1){
        # add cost of disambiguation
        cost_total <- cost_total+obs_info[obs_ind_temp,3]
        if(obs_info$status[obs_ind_temp]==1){
          # adjust based on true obstalce
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
        } else{
          
          # adjust based on false obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp]==1),3]-
                                                                    0.5*(obs_info[obs_ind_temp,3]+(dt/(1-obs_info[obs_ind_temp,4]))^(-log(1-obs_info[obs_ind_temp,4]))))
          Int_info[which(Int_info[,obs_ind_temp]==1),obs_ind_temp] <- 0 
        }
      } else{
        dist_temp <- rep(0,length(obs_ind_temp))
        for(i in 1:length(obs_ind_temp)){
          dist_temp[i] <- Dist_Euclidean(as.numeric(vertice_list[D_state,1:2]),obs_info[obs_ind_temp[i],1:2])
        }
        obs_ind_temp2 <- obs_ind_temp[which.min(dist_temp)]
        # add cost of disambiguation
        cost_total <- cost_total + obs_info[obs_ind_temp2,3]
        if (obs_info$status[obs_ind_temp2]==1){
          # true obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
        } else{
          
          # false obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3]-
                                                                     0.5*(obs_info[obs_ind_temp2,3]+(dt/(1-obs_info[obs_ind_temp2,4]))^(-log(1-obs_info[obs_ind_temp2,4]))))
          Int_info[which(Int_info[,obs_ind_temp2]==1),obs_ind_temp2] <- 0
        }
      }
      G_ed <- graph_from_data_frame(df_edge_ed,directed = FALSE)
    }
  }
  return(output_final)
}

####################################
####################################
####################################
####################################
####################################

######################################
edit_distance <- function(A, B) {
  nz <- length(A)
  mz <- length(B)
  
  dp <- matrix(0, nrow = nz + 1, ncol = mz + 1)
  dp[, 1] <- 0:nz
  dp[1, ] <- 0:mz
  
  for (iz in 1:nz) {
    for (jz in 1:mz) {
      if (A[iz] == B[jz]) {
        dp[iz + 1, jz + 1] <- dp[iz, jz]
      } else {
        dp[iz + 1, jz + 1] <- 1 + min(
          dp[iz, jz + 1],
          dp[iz + 1, jz],
          dp[iz, jz]
        )
      }
    }
  }
  
  return(dp[nz + 1, mz + 1])
}

sim <- function(A, B) {
  dzz <- edit_distance(A, B)
  return(1 - dzz / max(length(A), length(B)))
}


######################################
######################################
######################################
######################################


########################

#compare algorithm 

########################

x <- 100; y <- 100; r <- 4.5
# begin the loop to travel from s to t
s <- 10151
t <- 152
vertice_list <- Lattice_Vertices(x,y)
G_original <- Graph_Discretized(x,y)


#####################################
#####################################
#test the algorithm 

#obs_info <- read.csv("pattern/rho_0.5_1.csv")
#dt <- 25


DT_Alg_Save <- function(obs_info){
  output_Ginfo <- Update_graph_intersect_DT(G_original, x, y, obs_info, r)
  
  G_ed <- output_Ginfo$G_info
  Int_info <- output_Ginfo$Int_info
  
  ##from to cost
  df_edge_ed <- as_data_frame(G_ed, what="edges")
  
  explored_node <- s

  #used for total traversal costs calculation 
  length_total <- 0 # record Euclidean length
  cost_total <- 0 # record disambiguation cost 
  
  #for traversal overall
  seq_of_edges <- c()
  
  #for margin
  margin_edges <- c()
   
  #Kendall(ranking edges each time, compute Knedall after finishing) 
  #indexed by nodes 
  ranking_list_edges <- list()
  #rankings_algo1[["q"]] <- is_boundary[order(edge_scores)]
  
  
  
  ############
  #while loop 
  ###############
  
  while(explored_node != t){
  
  ########
  #DT finding optimal action 
  ########
  #DT 
    is_boundary <- which(df_edge_ed[[1]] == explored_node | df_edge_ed[[2]] == explored_node)
    outside_vertex <- ifelse(df_edge_ed[[1]][is_boundary] == explored_node, df_edge_ed[[2]][is_boundary], df_edge_ed[[1]][is_boundary])
    
    all_distances_to_t <- distances(
      G_ed, V(G_ed),    # all nodes as source
      which(vertex.attributes(G_ed)$name==as.character(t)),
      weights = df_edge_ed$Cost,
      algorithm = "dijkstra"
    )
    
    edge_scores <- df_edge_ed$Cost[is_boundary] + 
      as.numeric(all_distances_to_t[ as.numeric(outside_vertex) ])
    
    #edges used for traversal overall
    action_edge <- is_boundary[which.min(edge_scores)]
    seq_of_edges <- c(seq_of_edges, action_edge)
    
    #for margin
    sorted_scores <- sort(edge_scores)
    margin_edges <- c(margin_edges, abs(sorted_scores[2]- sorted_scores[1]))
    
    #for lists for Kendall 
    #indexed by nodes (in character)
    ranking_list_edges[[as.character(explored_node)]] <- is_boundary[order(edge_scores)]
  
  #working
  
    
  
  ####  
  ####now the whether disambiguation 
  ####  
  #for DT
  #######
    obs_ind_temp <- which(Int_info[action_edge,]==1)
    
    if(length(obs_ind_temp) !=0){
      if (length(obs_ind_temp)==1){
        # add cost
        cost_total <- cost_total+obs_info[obs_ind_temp,3]
        
        if(obs_info$status[obs_ind_temp]==1){
          # adjust based on true obstalce
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
          #ends here, current node kept and no length added 
          
        }else{
          #false 
          dt_nonfixed <- Dist_Euclidean(as.numeric(obs_info[obs_ind_temp,1:2]),c(50,1))
          # adjust based on false obstacle
          
          #edges index which which(Int_info[,obs_ind_temp]==1) is true, corss the obsacle
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp]==1),3]-
                                                                    0.5*(obs_info[obs_ind_temp,3]+(dt_nonfixed/(1-obs_info[obs_ind_temp,4]))^(-log(1-obs_info[obs_ind_temp,4]))))
          Int_info[which(Int_info[,obs_ind_temp]==1),obs_ind_temp] <- 0 
          
          #adding the length
          will_be_nex_node <- as.numeric(outside_vertex[match(action_edge, is_boundary)])
          
          edge_length <- Dist_Euclidean(as.numeric(vertice_list[explored_node, 1:2]),as.numeric(vertice_list[will_be_nex_node,1:2]))
          length_total <- length_total + edge_length
          
        }
        
      }else{
      ##########
      #multiple obstacle crossed 
      ##########  
        #keep the Euclidean distance, which is cloested 
        dist_temp <- rep(0,length(obs_ind_temp))
        for(i in 1:length(obs_ind_temp)){
          dist_temp[i] <- Dist_Euclidean(as.numeric(vertice_list[as.numeric(df_edge_ed[action_edge, 1]), 1:2]),obs_info[obs_ind_temp[i],1:2])
        }
        obs_ind_temp2 <- obs_ind_temp[which.min(dist_temp)]
        # add cost of disambiguation
        cost_total <- cost_total+obs_info[obs_ind_temp2,3]
        if (obs_info$status[obs_ind_temp2]==1){
          # true obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
          #ends here, current node kept and no length added 
          
        } else{
          dt_nonfixed <- Dist_Euclidean(as.numeric(obs_info[obs_ind_temp2,1:2]),c(50,1))
          # false obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3]-
                                                                     0.5*(obs_info[obs_ind_temp2,3]+(dt_nonfixed/(1-obs_info[obs_ind_temp2,4]))^(-log(1-obs_info[obs_ind_temp2,4]))))
          Int_info[which(Int_info[,obs_ind_temp2]==1),obs_ind_temp2] <- 0
          
          #adding the length
          will_be_nex_node <- as.numeric(outside_vertex[match(action_edge, is_boundary)])
          
          edge_length <- Dist_Euclidean(as.numeric(vertice_list[explored_node, 1:2]),as.numeric(vertice_list[will_be_nex_node,1:2]))
          length_total <- length_total + edge_length
          
        }
        
      }
      
    }else{
      #does not cross any obstacles 
      #add length 
     
      will_be_nex_node <- as.numeric(outside_vertex[match(action_edge, is_boundary)])
      
      #adding the length 
      edge_length <- Dist_Euclidean(as.numeric(vertice_list[explored_node, 1:2]),as.numeric(vertice_list[will_be_nex_node,1:2]))
      length_total <- length_total + edge_length
      #next 
      
      explored_node <- will_be_nex_node
    }
    
    G_ed <- graph_from_data_frame(df_edge_ed,directed = FALSE)
  
  
  }
    
  
  ###############################
  #summarize the data for return 
  ###############################
  
  total_trav_cost <- length_total + cost_total
  
    
    
  avg_margin_edges <- mean(margin_edges)
    
  output_final <- list(trav_cost = total_trav_cost, 
                         rank_edges = ranking_list_edges, 
                         full_edge = seq_of_edges, 
                         avg_margin = avg_margin_edges)
    
  return(output_final)   
    
  
}



###improved one ends

####################################
####################################
#end test
####################################
####################################










##################################
#fDT SAve
##################################

fDT_Alg_Save <- function(obs_info,dt){
  output_Ginfo <- Update_graph_intersect_DT_fixed(G_original, x, y, obs_info, r, dt)
  
  G_ed <- output_Ginfo$G_info
  Int_info <- output_Ginfo$Int_info
  
  ##from to cost
  df_edge_ed <- as_data_frame(G_ed, what="edges")
  
  explored_node <- s
  
  #used for total traversal costs calculation 
  length_total <- 0 # record Euclidean length
  cost_total <- 0 # record disambiguation cost 
  
  #for traversal overall
  seq_of_edges <- c()
  
  #for margin
  margin_edges <- c()
  
  #Kendall(ranking edges each time, compute Knedall after finishing) 
  #indexed by nodes 
  ranking_list_edges <- list()
  #rankings_algo1[["q"]] <- is_boundary[order(edge_scores)]
  
  
  
  ############
  #while loop 
  ###############
  
  while(explored_node != t){
    
    ########
    #DT finding optimal action 
    ########
    #DT 
    is_boundary <- which(df_edge_ed[[1]] == explored_node | df_edge_ed[[2]] == explored_node)
    outside_vertex <- ifelse(df_edge_ed[[1]][is_boundary] == explored_node, df_edge_ed[[2]][is_boundary], df_edge_ed[[1]][is_boundary])
    
    all_distances_to_t <- distances(
      G_ed, V(G_ed),    # all nodes as source
      which(vertex.attributes(G_ed)$name==as.character(t)),
      weights = df_edge_ed$Cost,
      algorithm = "dijkstra"
    )
    
    edge_scores <- df_edge_ed$Cost[is_boundary] + 
      as.numeric(all_distances_to_t[ as.numeric(outside_vertex) ])
    
    #edges used for traversal overall
    action_edge <- is_boundary[which.min(edge_scores)]
    seq_of_edges <- c(seq_of_edges, action_edge)
    
    #for margin
    sorted_scores <- sort(edge_scores)
    margin_edges <- c(margin_edges, abs(sorted_scores[2]- sorted_scores[1]))
    
    #for lists for Kendall 
    #indexed by nodes (in character)
    ranking_list_edges[[as.character(explored_node)]] <- is_boundary[order(edge_scores)]
    
    #working
    
    
    
    ####  
    ####now the whether disambiguation 
    ####  
    #for DT
    #######
    obs_ind_temp <- which(Int_info[action_edge,]==1)
    
    if(length(obs_ind_temp) !=0){
      if (length(obs_ind_temp)==1){
        # add cost
        cost_total <- cost_total+obs_info[obs_ind_temp,3]
        
        if(obs_info$status[obs_ind_temp]==1){
          # adjust based on true obstalce
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
          #ends here, current node kept and no length added 
          
        }else{
          #false 
          # adjust based on false obstacle
          
          #edges index which which(Int_info[,obs_ind_temp]==1) is true, corss the obsacle
          df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp]==1),3]-
                                                                    0.5*(obs_info[obs_ind_temp,3]+(dt/(1-obs_info[obs_ind_temp,4]))^(-log(1-obs_info[obs_ind_temp,4]))))
          Int_info[which(Int_info[,obs_ind_temp]==1),obs_ind_temp] <- 0 
          
          #adding the length
          will_be_nex_node <- as.numeric(outside_vertex[match(action_edge, is_boundary)])
          
          edge_length <- Dist_Euclidean(as.numeric(vertice_list[explored_node, 1:2]),as.numeric(vertice_list[will_be_nex_node,1:2]))
          length_total <- length_total + edge_length
          
        }
        
      }else{
        ##########
        #multiple obstacle crossed 
        ##########  
        #keep the Euclidean distance, which is cloested 
        dist_temp <- rep(0,length(obs_ind_temp))
        for(i in 1:length(obs_ind_temp)){
          dist_temp[i] <- Dist_Euclidean(as.numeric(vertice_list[as.numeric(df_edge_ed[action_edge, 1]), 1:2]),obs_info[obs_ind_temp[i],1:2])
        }
        obs_ind_temp2 <- obs_ind_temp[which.min(dist_temp)]
        # add cost of disambiguation
        cost_total <- cost_total+obs_info[obs_ind_temp2,3]
        if (obs_info$status[obs_ind_temp2]==1){
          # true obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
          #ends here, current node kept and no length added 
          
        } else{
          # false obstacle
          df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- pmax(0, df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3]-
                                                                     0.5*(obs_info[obs_ind_temp2,3]+(dt/(1-obs_info[obs_ind_temp2,4]))^(-log(1-obs_info[obs_ind_temp2,4]))))
          Int_info[which(Int_info[,obs_ind_temp2]==1),obs_ind_temp2] <- 0
          
          #adding the length
          will_be_nex_node <- as.numeric(outside_vertex[match(action_edge, is_boundary)])
          
          edge_length <- Dist_Euclidean(as.numeric(vertice_list[explored_node, 1:2]),as.numeric(vertice_list[will_be_nex_node,1:2]))
          length_total <- length_total + edge_length
          
        }
        
      }
      
    }else{
      #does not cross any obstacles 
      #add length 
      
      will_be_nex_node <- as.numeric(outside_vertex[match(action_edge, is_boundary)])
      
      #adding the length 
      edge_length <- Dist_Euclidean(as.numeric(vertice_list[explored_node, 1:2]),as.numeric(vertice_list[will_be_nex_node,1:2]))
      length_total <- length_total + edge_length
      #next 
      
      explored_node <- will_be_nex_node
    }
    
    G_ed <- graph_from_data_frame(df_edge_ed,directed = FALSE)
    
    
  }
  
  
  ###############################
  #summarize the data for return 
  ###############################
  
  total_trav_cost <- length_total + cost_total
  
  
  
  avg_margin_edges <- mean(margin_edges)
  
  output_final <- list(trav_cost = total_trav_cost, 
                       rank_edges = ranking_list_edges, 
                       full_edge = seq_of_edges, 
                       avg_margin = avg_margin_edges)
  
  return(output_final)   
  
  
}



