generate_points <-
function(W, unlabeled, labeling) {
	n <- nrow(W);
	N <- length(labeling[-unlabeled]);
	N_pos <- length(which(labeling[-unlabeled] == 1));
	pos_rate <- N_pos/N;					
	p <- length(unlabeled);
	# generating random labels for unlabeled nodes according to the distribution B(p, N_pos/N) 			
	y_rand <- generate_labels(p, pos_rate);		
	labeling[unlabeled] <- y_rand;			
	labeling[labeling <=0] <- 0;
	pos_vect <- W %*% labeling;	
	neg_vect <- W %*% (1-labeling)
	res <- list(pos_vect = pos_vect, neg_vect = neg_vect);
	return (res);
}
