DFSCOLOR <- numeric(0)
DFSBACKEDGE <- numeric(0)
ORDERVERTICES <- numeric(0)

DFSVisit <- function(A,i){
  DFSCOLOR[i] <<- 1
  for (j in 1:dim(A)[1]){
    if(A[i,j]!=0){
      if (DFSCOLOR[j]==0){
        DFSVisit(A,j)
      }
      else{
        if(DFSCOLOR[j]==1){
          ## It's a back edge: list to remove
          DFSBACKEDGE[i,j] <<- 1
        }
      }
    }
  }
  DFSCOLOR[i] <<- 2
  ORDERVERTICES <<- c(i, ORDERVERTICES)
}

DFS <- function(A){
  ## A = Adjacency matrix sorted such that the producers are at the start of the nodes.
  S <- dim(A)[1]
  DFSCOLOR <<- rep(0,S)
  DFSBACKEDGE <<- matrix(0,S,S)
  ORDERVERTICES <<- numeric(0)
  for (i in 1:S){
    if (DFSCOLOR[i]==0){
      DFSVisit(A,i)
    }
  }
  return (A-DFSBACKEDGE)
}
