generate_labels <-
function(n, pos_rate)
{
	set.seed(as.integer(substring(Sys.time(), 15 ,16)));
	# generates n random number between 0 and 1				
	y_rand <- runif(n, 0, 1)
	y_rand_tmp <- y_rand;
	y_rand[y_rand_tmp <= pos_rate] <- 1;
	y_rand[y_rand_tmp > pos_rate] <- -1;	
	return(y_rand);					
}
