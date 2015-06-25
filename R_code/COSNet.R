COSNet <-
function(W, labeling, cost=0){	
	unlabeled <- which(labeling == 0);	
	pos.labels <- sum(labeling>0)
	### Step 1		            	
   #computes points for the separation phase   	
   points <- generate_points(W, unlabeled, labeling);
	### Step 2
	# compute optimal parameters	
	opt_parameters <- optimizep(points$pos_vect[-unlabeled], 
	                  points$neg_vect[-unlabeled], labeling[-unlabeled]);
	alpha <- opt_parameters$alpha;
	c <- opt_parameters$c;
	Fscore <- opt_parameters$Fscore;	
	rm(opt_parameters);	
	### Step 3
	#compute predictions and scores for unlabeled nodes
	res <- runSubnet(W, labeling, alpha, c,  cost);
	out <- list(alpha = alpha, c = c, Fscore = Fscore, pred = res$state, 
	            scores = res$scores, iter = res$iter);
	return(out);			
}

.onLoad <- function(libname=.libPaths(), pkgname="COSNet")
       library.dynam("COSNet", pkgname, libname);
	   
.onAttach <- function(libname=.libPaths(), pkgname="COSNet")
 cat("COSNet: Cost-Sensitive algorithm for binary classification in graphs.\n");
