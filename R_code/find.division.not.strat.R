find.division.not.strat <-
function(vett, n_fold)
{	
	n <- length(vett);
	if(n > n_fold){
		# get the actual seconds to set seed
		set.seed(as.integer(substring(Sys.time(), 18 ,19)));	
		# makes a permutation a random permutations of indexes before dividing its in various folds	
		d <- sample(vett, n);
		k <- floor(n/n_fold);	
		fold_list <- lapply(1:n_fold, function(x, d, k)
		              {return(d[(k*(x-1)+1):(k*x)])}, d=d, k=k);
		fold_list[[n_fold]] <- d[(k*(n_fold-1)+1):n];	
		return(fold_list);
	}
	else
		stop(paste(sep="", "The number of folds must be less or equal to the ",
		                     "length of the vector to be partitioned"));
}
