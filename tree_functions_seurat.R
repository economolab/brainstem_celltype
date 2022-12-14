#' @include generics.R
#'
NULL

cluster.ape <- paste(
  "Cluster tree functionality requires 'ape'",
  "please install with 'install.packages('ape')'"
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Depth first traversal path of a given tree
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
# @param path Path through the tree (for recursion)
# @param include.children Include children in the output path
# @param only.children Only include children in the output path
# @return Returns a vector representing the depth first traversal path
#
DFT <- function(
  tree,
  node,
  path = NULL,
  include.children = FALSE,
  only.children = FALSE
) {
  if (only.children) {
    include.children = TRUE
  }
  children <- which(x = tree$edge[, 1] == node)
  child1 <- tree$edge[children[1], 2]
  child2 <- tree$edge[children[2], 2]
  if (child1 %in% tree$edge[, 1]) {
    if (!only.children) {
      path <- c(path, child1)
    }
    path <- DFT(
      tree = tree,
      node = child1,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child1)
    }
  }
  if (child2 %in% tree$edge[, 1]) {
    if (!only.children) {
      path <- c(path, child2)
    }
    path <- DFT(
      tree = tree,
      node = child2,
      path = path,
      include.children = include.children,
      only.children = only.children
    )
  } else {
    if (include.children) {
      path <- c(path, child2)
    }
  }
  return(path)
}

# Function to return all internal (non-terminal) nodes in a given tree
#
# @param tree Tree object (from ape package)
#
# @return Returns a vector of all internal nodes for the given tree
#
GetAllInternalNodes <- function(tree) {
  return(c(tree$edge[1, 1], DFT(tree = tree, node = tree$edge[1, 1])))
}

# Function to get all the descendants on a tree of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants of the given node
#
GetDescendants <- function(tree, node, curr = NULL) {
  if (is.null(x = curr)) {
    curr <- vector()
  }
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  curr <- c(curr, daughters)
  w <- which(x = daughters >= length(x = tree$tip))
  if (length(x = w) > 0) {
    for (i in 1:length(x = w)) {
      curr <- GetDescendants(tree = tree, node = daughters[w[i]], curr = curr)
    }
  }
  return(curr)
}

# Function to get all the descendants on a tree left of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants left of the given node
#
GetLeftDescendants <- function(tree, node) {
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  if (daughters[1] <= (tree$Nnode + 1)) {
    return(daughters[1])
  }
  daughter.use <- GetDescendants(tree, daughters[1])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

# Function to get all the descendants on a tree right of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants right of the given node
#
GetRightDescendants <- function(tree, node) {
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  if (daughters[2] <= (tree$Nnode + 1)) {
    return(daughters[2])
  }
  daughter.use <- GetDescendants(tree = tree, node = daughters[2])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

# Merge childen of a node
#
# Merge the childen of a node into a single identity class
#
# @param object Seurat object
# @param node.use Merge children of this node
# @param rebuild.tree Rebuild cluster tree after the merge?
# @param ... Extra parameters to BuildClusterTree, used only if rebuild.tree = TRUE
#
# @seealso \code{BuildClusterTree}
#
#
# @examples
# data("pbmc_small")
# PlotClusterTree(object = pbmc_small)
# pbmc_small <- MergeNode(object = pbmc_small, node.use = 7, rebuild.tree = TRUE)
# PlotClusterTree(object = pbmc_small)
#
MergeNode <- function(object, node.use, rebuild.tree = FALSE, ...) {
  CheckDots(..., fxns = 'BuldClusterTree')
  object.tree <- object@cluster.tree[[1]]
  node.children <- DFT(
    tree = object.tree,
    node = node.use,
    include.children = TRUE
  )
  node.children <- intersect(x = node.children, y = levels(x = object@ident))
  children.cells <- WhichCells(object = object, ident = node.children)
  if (length(x = children.cells > 0)) {
    object <- SetIdent(
      object = object,
      cells.use = children.cells,
      ident.use = min(node.children)
    )
  }
  if (rebuild.tree) {
    object <- BuildClusterTree(object = object, ...)
  }
  return(object)
}

# Function to check whether a given node in a tree has a child (leaf node)
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns a Boolean of whether the given node is connected to a terminal leaf node

NodeHasChild <- function(tree, node) {
  children <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  return(any(children %in% tree$edge[, 2] && !children %in% tree$edge[, 1]))
}

# Function to check whether a given node in a tree has only children(leaf nodes)
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns a Boolean of whether the given node is connected to only terminal leaf nodes

NodeHasOnlyChildren <- function(tree, node) {
  children <- tree$edge[which(x = tree$edge[, 1] == node), ][, 2]
  return(!any(children %in% tree$edge[, 1]))
}




################################################################################
# MY OWN WRITTEN FUNCTIONS 
################################################################################

# a function to get all terminal leaf nodes

GetLeavesOnly<-function(object,tree, node){
  if(NodeHasOnlyChildren(tree, node)==TRUE){
    descendants<-GetDescendants(tree, node)  # subtract 1 from it because it adds 1 to create all the nodes
  }
  else{
    descendants<-GetDescendants(tree, node)
    for (i in descendants){
      if (i>length(levels(object))){
        descendants<-descendants [! descendants %in% i]              # Remove list element with !=
      }
    }
    
    
  }
  return(descendants-1)
}






