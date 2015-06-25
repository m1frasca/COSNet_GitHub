cosnet.cross.validation <-
function(labels, W, nfolds, cost=0.0001)
{		
	n <- nrow(labels);
	m <- ncol(labels);	
	genes <- rownames(W);
	predictions <- matrix(0, nrow = n, ncol = m);        
	scores <- matrix(0, nrow = n, ncol = m);        	
	rownames(scores) <- rownames(predictions) <- genes;
	labels[labels <= 0] <- -1;	
	tmp.sum <- apply(W, 1, sum);
	isolated <- which(tmp.sum == 0);	
	labels.noiso <- labels;
	if(length(isolated) > 0){
		W <- W[-isolated, -isolated]
		labels.noiso <- as.matrix(labels.noiso[-isolated, ])
    	n <- nrow(W)
    }
	for(i in 1:m){		
		labeling <- as.vector(labels.noiso[, i]);	
		fold_indexes <- find.division.strat(labeling, 1:n, nfolds);	
		nodes <- rownames(labels.noiso);
		for(j in 1:length(fold_indexes)){		
			test <- as.vector(fold_indexes[[j]]);				
			labelingTMP <- labeling;
			labelingTMP[test] <- 0;			
			out <- COSNet(W, labelingTMP, cost)	
			predictions[nodes[sort(test)], i] <- out$pred;
			scores[nodes[sort(test)], i] <- out$scores;		
		}
	}
	classes <- colnames(labels);
	if(length(classes) > 0)
		colnames(scores) <- colnames(predictions) <- classes	
	out <- list(labels = labels, predictions = predictions, scores = scores);
	return (out);
}
