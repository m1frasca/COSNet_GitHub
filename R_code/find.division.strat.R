find.division.strat <-
function(labels, vett, n_fold)
{
	n <- length(vett);
	if(n > n_fold){
		# get the actual seconds to set seed
		set.seed(as.integer(substring(Sys.time(), 18 ,19)));		
		pos = which(labels > 0);
		npos <- length(pos)	
		vett.pos <- vett[pos];
		if(npos >= n_fold){	#stratified partition			
			d <- vett[-pos];
			d <- sample(d, (n-npos));
			k <- round((n-npos)/n_fold);
			kp <- round((npos)/n_fold);	
			fold_indices <- lapply(1:n_fold, function(x, d, k, kp, vett.pos)
			                   { return(c(d[(k * (x - 1) + 1):(k * x)], 
			                            pos[(kp * (x - 1) + 1):(kp * x)])) }, 
			                  d = d, k = k, kp = kp, vett.pos = vett.pos);
			fold_indices[[n_fold]] <- c( d[(k*(n_fold-1)+1):(n-npos)], 
			                              vett.pos[(kp*(n_fold-1)+1):npos]);		
		}
		else # not stratified 
			fold_indices <- find.division.not.strat(1:n, n_fold)

		return (fold_indices);	
	}
	else
		stop(paste(sep="", "The number of folds must be less or equal to the ",
		                     "length of the vector to be partitioned"));	

}
