optimizep <-
function(pos_vect, neg_vect, training_labels) {	
	alpha <- 0; c <- 0;
	Fscore <- -1.0; pos_halfplane <- 0;
	N <- length(training_labels);
	# function to optimize the model parameters										
	out2 <- .C("error_minimization",as.double(pos_vect),as.double(neg_vect),
	         as.integer(training_labels), as.integer(N), as.double(alpha), 
	          as.double(c), as.double(Fscore), as.integer(pos_halfplane), 
	           PACKAGE="COSNet");	
	res <- list(alpha = out2[[5]], c = out2[[6]], Fscore = out2[[7]], 
	       pos_half = out2[[8]]);
#	res[[1]] <- out2[[5]];	res[[2]] <- out2[[6]];
#	res[[3]] <- out2[[7]];	res[[4]] <- out2[[8]];
	if(res$pos_half < 0)
		return (res)
	else	
		out2 <- optimize_pos_above(pos_vect, neg_vect, training_labels, res)	
	return (out2);						
}
