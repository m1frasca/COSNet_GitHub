reg_data <-
function(W, theta, eta, M, m, pos_num) {	
	N <- length(theta);
	a <- 1/(M - m);
	b <- -m/(M - m);
	W <- W - (2 * eta * (a^2));
	theta <- theta + eta*a*(2*b*(N - 1) - 2*pos_num + 1);	
	diag(W) <- 0;	
	res <- list(W = W, theta = theta);
	return (res);
}
