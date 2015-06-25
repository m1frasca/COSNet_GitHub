runSubnet <-
function(W, labeling, alpha_value, c_value, cost) { 
	unlabeled <- which(labeling == 0);			
	p <- length(unlabeled);
	n <- length(labeling);
	# number of positives in training
	N_pos <- length(which(labeling[-unlabeled] == 1));	
	# defining neuron values by considering the angle "alpha_value"	
	pos_state <- sin(alpha_value);
	neg_state <- -cos(alpha_value);					
	labeling[labeling == 1] <- pos_state;
	labeling[labeling == -1] <- neg_state;
	# thresholds for unlabeled nodes deriving from optimization phase
	thresholds <- rep(c_value, p);	
	# computing the contribution for unlabeled nodes from those labeled 
	thresholds.top <- W[unlabeled, ] %*% labeling;
	w_u <- W[unlabeled, unlabeled];				
	# setting final thresholds for unlabeled nodes			
	theta <- unlist(thresholds.top);		        	
	theta <- thresholds - theta;	
	diag(w_u) <- 0;
	# computing data due to regularization only when cost > 0. The regularization changes according to the position of the positive half-plane
	if( (cost > 0) && (alpha_value > pi/4)) {  
		eta <- (cost) * abs(tan( (alpha_value - pi/4)*2 ));
		data <- reg_data(w_u, theta, eta, pos_state, neg_state,
		       ceiling((N_pos/(n-p))*p))
		w_u <- data$W; 
	    theta <- data$theta;
		rm(data);	
    }	
	state <- vector("numeric", length = p);	
	genes <- rownames(W);
	scores <- vector("numeric", length = p);	
	# network dynamics		
	iter=0;	
	repeat
	{					
		iter <- iter + 1;
		out1 <- .C("global_update_under", as.double(w_u), as.double(theta),
		        as.double(pos_state), as.double(neg_state), as.double(state), 
		         as.double(scores), as.integer(p), PACKAGE="COSNet")
		new_state <- out1[[5]];		
		n_changes <- sum(new_state != state);
		if(n_changes == 0){
			scores <- out1[[6]]				
			break
		}	
		else		
			state <- new_state							
	}
	names(state) <- genes[unlabeled];
	names(scores) <- genes[unlabeled];	
	th.res <- (pos_state + neg_state)/2
	state[which(state > th.res)] <- 1
	state[which(state <= th.res)] <- -1
	res <- list(state = state, scores = scores, iter=iter);	
	return (res); 
}
