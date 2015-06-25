optimize_pos_above <-
function(pos_vect, neg_vect, training_labels, res)
{		
    n <- length(pos_vect);
    zero.abscissa <- which(pos_vect == 0); # nodes with no positive neighbors
    n.point <- length(zero.abscissa)
		#*************************************************************************
		# Step 1******************************************************************
    ordinate.zero.abscissa <- neg_vect[zero.abscissa]; # ordinates of points with no positive neighbors
		# indices (in zero.ascissa) of positive nodes with no positive neighbors 
    positive.ordinate.zero.abscissa<-which(training_labels[zero.abscissa] > 0) 
		# permutation of "ordinate.zero.abscissa" corresponding to the increasing sorted sequance
    ordinate.zero.abscissa.order <- order(ordinate.zero.abscissa);	
		# sorted "ordinate.zero.abscissa"
    ordinate.zero.abscissa.sort <- sort(ordinate.zero.abscissa);	
    n.pos <- length(positive.ordinate.zero.abscissa);
    best.F <- 0; best.index <- 0; TP  <- FP <- 0; FN <- n.point;
		# selecting the indices of positive nodes with absissa 0 in the sorted vector of the corresponding ordinates. 
    ord.pos.ind <- sort(which(ordinate.zero.abscissa.order %in% 
                    positive.ordinate.zero.abscissa))	
      # if there is at least a positive node with ascissa 0
    if(n.pos > 0){
			# scans all the lines passing for a positive point and with negative slopes, and for each of it computes the corresponding F-score
			# by considering the positive half-plane below the line			
        for(i in 1:n.pos){
            act.pos <- ord.pos.ind[i] # position in the sorted oridinates vector
				# the ordinate of the actual point
            act.ord <- ordinate.zero.abscissa.sort[act.pos]
				# it means that the line with negative slope passing for this point classifies in wrong all the negatives under the 
				# line (the positive hal-plane is below the line). 
            TP <- i; FN <- n.pos - TP; FP <- act.pos - TP
            if((TP + FP)!=0 & TP + FN != 0)
                act.F <- (2*(TP/(TP + FP))*(TP/(TP + FN)))/(TP/(TP + FP) + 
                         TP/(TP + FN))
            else   act.F <- 0					
            if(act.F > best.F)	{
                best.F <- act.F;	best.index <- act.pos
            } 
		}			
			# *******************************************************************************************
			# Step 2*************************************************************************************
			# after found the positive point P (corresponding to best.index) with abscissa 0 corresponding to the best Fscore,
			# it computes the Fscore corresponding to the lines crossing P and all the ramaining points with ascissa > 0. 
        x_1 <- 0; y_1 <- ordinate.zero.abscissa.sort[best.index] 
        slopes <- vector();
			# number of points with abscissa > 0
        n.non.zero <- n - n.point
		# computing slopes of all the lines crossing P and all the ramaining points with ascissa > 0
        for (j in 1:n.non.zero){
            x_2 <- (pos_vect[-zero.abscissa])[j];
            y_2 <- (neg_vect[-zero.abscissa])[j]
            slopes[j] <- (y_2 - y_1)/((x_2 - x_1))
        }
        sub.labels <- training_labels[-zero.abscissa]
        slopes.order <- order(slopes);  slopes.sort <- sort(slopes)
			# positive points among those with abscissa > 0
        positive.up <- which(training_labels[-zero.abscissa] > 0)
        TP <- which(ord.pos.ind == best.index)
        FN <- n.pos - TP; FP <- best.index - TP;			
        best.F<-(2*(TP/(TP + FP))*(TP/(TP + FN)))/(TP/(TP + FP) + TP/(TP + FN));
        FN <- FN + length(positive.up)			
        best.slope <- tan(1.58) # almost vertical line as default
        number.slopes.filtered <- 0
        # computing the best Fscores among these lines scannig in an incresing order of slopes 
        for(j in 1:n.non.zero){
            if(slopes.sort[j] > -1) { # considers just the lines with slope in ]-Inf, -1]
                next # angle not acceptable
            }
            else	{
                number.slopes.filtered <- number.slopes.filtered + 1
                if(sub.labels[slopes.order[j]] > 0){
                    TP <- TP + 1;	FN <- FN - 1					
                }	
                else	FP <- FP + 1					
                if((TP + FP)!=0 & TP + FN != 0)
                    act.F <- (2*(TP/(TP + FP))*(TP/(TP + FN)))/(TP/(TP + FP)
                              + TP/(TP + FN))
                else act.F <- 0
                if(act.F > best.F) {
                    best.F <- act.F;	best.slope <- slopes.sort[j];
                }			
            }
        }
			# *******************************************************************************************
			# Step 3 *************************************************************************************
			#after found the optimal angle, computing intercept
        intercepts <- vector();
        for (j in 1:n){
            x_2 <- pos_vect[j]; y_2 <- neg_vect[j]
            intercepts[j] <- y_2 - best.slope*x_2
        }
        intercepts.order<-order(intercepts); intercepts.sort <- sort(intercepts)			
        FP <- 0; TP <- 0; FN <- sum(training_labels>0)						
        best.F <- 0; best.intercept <- 0;
        for(j in 1:n){				
            if(training_labels[intercepts.order[j]] > 0)	{				
                TP <- TP + 1; FN <- FN - 1					
            }	
            else	FP <- FP + 1							
            if(TP > 0)
                act.F <- (2*(TP/(TP + FP))*(TP/(TP + FN)))/(TP/(TP + FP) 
                          + TP/(TP + FN))
            else 	act.F <- 0
            if(act.F > best.F) {
                best.F <- act.F;	best.intercept <- intercepts.sort[j]
             }		
        }
        res$alpha <- atan(best.slope)+pi;
        res$c <- -best.intercept*cos(res$alpha);
        res$Fscore <- best.F
        return (res)
    }
    else	{
			# No positive points with abscissa = 0
        return (res)		
    }
		
}
