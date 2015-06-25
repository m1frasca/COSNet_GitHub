

#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include <math.h>
#include <float.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

//typedef unsigned int R_NativePrimitiveArgType;

/* 
  This function computes, for each index "i", the weighted sum of its positive and negative neighbors. The pairwise weights are contained in the
   matix "W", and the neighborhood of a node i is given by the indices j  for which the i,j-th component of W is greater than 0.   

  INPUT
  W:	square matrix, whose components are in the [0,1] interval. The i,j-th component is the weight between node i and node j.
		 The components of the diagonal are 0.
  labels: vector of node labels : 1 for positive, -1 for negative items
  pos_vect: output vector for positive neighborhood
  neg_vect: output vector for negative neighborhood
  n: number of items (dimension of W and labels) 
*/

/****************************************************************************** 
	Function to swap the value of x and y
	INPUT
	x : first integer to swap
	y : second integer to swap
********************************************************************************/
void swap(int *x,int *y)
{
   int temp; temp = *x; *x = *y; *y = temp;
}


/*****************************************************************
Function to realize a random in place permutation of element in vector "a". The input variable "size" represent the size of "a"

INPUT
  a: integer vector. In out variable
  size: integer representing the number of elements in a
****************************************************************************/ 
void permute(int * a, int size)
{
	register int i;
	int change;
	srand(time(0));
	for(i = 0; i < size; i++)
	{
		change = rand()%size;
		swap(&a[i], &a[change]);		
	}	
}


/********************************************************************************
 Function to realize the update of node "i" in the dynamics of Hopfield network. The update rule is the signum function of its 
   "internal energy" (score).

INPUT
  W: symmetric connection weight matrix 
  theta_k: activation threshold of neuron "k_"
  pos_state: positive activation value  
  neg_state: negative  activation value  
  k_:index of the current neuron
  state: input/output vector of current network  state
  scores: input/output vector containing the internal energy of neurons
  N: number of neurons

*/
void update_i_under(double *W, double theta_k, double pos_state, double neg_state, int k_, double *state, double *scores, int N)
{	
	double int_act = 0.0;
	register int  j;	
	for(j = 0; j < N; j++)				
		int_act = int_act + (W[ ( (k_)*(N) ) + j] * state[j]);	
	scores[k_] = int_act - theta_k;
	if(scores[k_] > 0)
		state[k_] = pos_state;
	else	
		state[k_] = neg_state;		
}


/******************************************************************************** 
	Function to asynchronously update all the neurons in the Hopfield network in a random order when the postive. 
	
INPUT
  W: symmetric connection weight matrix
  theta: vector of activation thresholds
  pos_state: positive activation value  
  neg_state: negative  activation value    
  state: input/output vector of current network  state
  scores: input/output vector containing the internal energy of neurons
  size: vector of the neurons scores 

*/
void global_update_under(double *W, double *theta, double *pos_state, double * neg_state, double * state, double * scores, int * size)
{
	register int i;
	const int size_ = *size;
	int update_sequence[size_];
	srand(time(0));	
	// randomly permuting updating sequence
	for(i = 0; i < size_; update_sequence[i] = i, i++);			
	permute(update_sequence, size_);	
	for(i = 0; i < size_; i++)	
		update_i_under(W, theta[update_sequence[i]], *pos_state, *neg_state, update_sequence[i], state, scores,*size);	

}


/******************************************************************
  Function to compute the number of positive elements in a vector "a". 

INPUT
	a : integer vector
	size : length of a

OUTPUT
	N_pos: Number of elements in a greater than 0
*/
int count_positives(int * a, int size)
{
	register int j;
	int N_pos = 0;
	for(j = 0 ; j < size ; j++)
	{
		if(a[j] > 0)
			N_pos++;		
	}
	return N_pos;
}


/******************************************************************  
  Function to compute the angles (formed with x axes) of each line crossing the origin and each point (pos_vect[i], neg_vect[i])

INPUT    
	pos_vect: vector in which position "i" contains the weighted sum of positive neighbors of item "i"
   neg_vect: vector in which position "i" contains the weighted sum of negative neighbors of item "i"
	order_thetas: output vector which will contain the indices from 0 to "size"-1
   thetas: output vector which will contain the computed angles
	size: number of points
*/
void compute_angles(double * pos_vect, double * neg_vect, int * order_thetas, double * thetas, int size)
{
	register int j;
	double m, x1, y1;
	for(j = 0; j < size; j++)
	{		
		order_thetas[j] = j;
		x1 = pos_vect[j];	y1 = neg_vect[j];
		// excluding lanes parallel to y axes
		if((x1) != 0)
		{	
			m = (double)y1/x1;			
			thetas[j] = atan(m);
		}	
		else			
			thetas[j] = (double)M_PI/2 - DBL_MIN;
	}
}

// computes the Fscore of "tp" true positive, "fn" false negative, "fp" false positive
double compute_F(int tp, int fn, int fp)
{
	double prec, recall, tmp_F;
	if(tp + fp != 0)
		prec = (double)tp/(tp + fp);
	else
		prec = 0;
	if(tp + fn != 0)
		recall = (double)tp/(tp + fn);	
	else
		recall = 0;
	if(prec + recall != 0)
		tmp_F = (2*prec*recall)/(prec+recall);
	else
		tmp_F = 0;
	return tmp_F;	
}

/* Function to compute the intercepts of each line whose angle (with x axis) is "theta_best" crossing each point "i" whose coordinates are
  pos_vect[i] and neg_vect[i].

INPUT    
	pos_vect: vector in which position "i" contains the weighted sum of positive neighbors of item "i"
	neg_vect: vector in which position "i" contains the weighted sum of negative neighbors of item "i"
	order_c_values: output vector which will contain the indices from 0 to "size"-1
   c_values: output vector which will contain the computed intercepts	
	theta_best: optimum angle previously computed
	size: number of points
*/	
void compute_c(double * pos_vect, double * neg_vect, int * order_c_values, double * c_values, double theta_best, int size)
{
	register int j;
	double x1, y1;
	for(j = 0; j < size; j++)
	{			
		order_c_values[j] = j;
		x1 = pos_vect[j];	y1 = neg_vect[j];
		c_values[j] = y1 - tan(theta_best)*x1;	
	}
}

/****************************************************************************** 
	Function to swap the value of x and y
	INPUT
	x : first double to swap
	y : second double to swap
********************************************************************************/
void swap_d(double *x,double *y)
{
   double temp;
   temp = *x;
   *x = *y;
   *y = temp;
}

/****************************************************************************** 
	Function to swap the value of x and y
	INPUT
	x : first integer to swap
	y : second integer to swap
*******************************************************************************
void swap(int *x,int *y)
{
   int temp;
   temp = *x;
   *x = *y;
   *y = temp;
}
*/
/***************************************************************************************************** 
  Partition Function in the quickSort Algorithm. It returns the indices j such that all the elements before j in the vector
	 a are less than or equal to those after j.	 
*/
int partition( double a[], int indices[], int l, int r) {
   register int i, j;
   double pivot;
   pivot = a[l];
   i = l; j = r+1;		
   while(1)
   {
   	do ++i; while( a[i] <= pivot && i <= r );
   	do --j; while( a[j] > pivot );
   	if( i >= j ) break;
   	swap_d(&a[i], &a[j]);
   	swap(&indices[i], &indices[j]);     	
   } 
   swap_d(&a[l], &a[j]);
   swap(&indices[l], &indices[j]);
   return j;
}

/**********************************************************************************
   QuickSort algorithm. It sorts in an increasing order the elements in the vector a.
*/ 
void quicksort( double a[], int indices[], int l, int r)
{
   int j;
   if( l < r ) 
   {   	
       j = partition( a, indices, l, r);
       quicksort( a, indices, l, j-1);
       quicksort( a, indices, j+1, r);
   }	
}


/*********************************************************************************************
 Function to check if tmp_mean > opt_mean. In such case it updates accordingly the values of opt_tp, opt_fn, opt_fp and theta_best. 
	INPUT
	tmp_hmean:	 actual computed value.
	opt_hmean:	 actual best value
	theta_best:	 actual best computed angle,
	opt_tp:		 actual value of true positive corresponding to opt_hmean
	opt_fp: 		 actual value of false positive corresponding to opt_hmean
	opt_fn: 		 actual value of false negative corresponding to opt_hmean
	tp: 			actual value of true postive
	fp:			actual value of false positive
	fn:			actual value of false negative	 
*/
void check_update(double tmp_hmean, double *opt_hmean, double *theta_best, double angle,
									 int *opt_tp, int *opt_fp, int *opt_fn, int tp, int fp, int fn){
	if(tmp_hmean > *opt_hmean){				
		*opt_hmean = tmp_hmean;
		*theta_best = angle;
		*opt_tp = tp;
		*opt_fn = fn;
		*opt_fp = fp;				
	}
}

/*********************************************************************************************
 Function to check whether the highest Fscore between the optimum Fscore when the positive half-plane is above the separation line
  (opt_hmean_over) and below the line (opt_hmean_under). The highest value is assigned to max_F and the corresponding angle to theta_best. 
  The variable "pos_halfplane" will contain 1 if the best value correspond to positive half-plane above the separation line, -1 otherwise
INPUT
  opt_hmean_over:		optimal value of Fscore when the positive half-plane is above the separation line
  opt_hmean_under:	optimal value of Fscore when the positive half-plane is below the separation line
  theta_best_over:	angle corresponding to the optimal Fscore when the positive half-plane is above the separation line
  theta_best_under:	angle corresponding to the optimal Fscore when the positive half-plane is above the separation line
  pos_halfplane:		output variable which will contain 1 if the best value of Fscore correspond to positive half-plane above the separation
  							 line, -1 otherwise
  max_F:					output variable which will contain the maximum Fscore
  theta_best:			output variable which will contain the angle corresponding to the maximum Fscore
*/
void check_halfplane(double opt_hmean_over, double opt_hmean_under, double theta_best_over, double theta_best_under, int *pos_halfplane, double *max_F, double *theta_best){
	if(opt_hmean_over > opt_hmean_under)
	{
		*pos_halfplane = 1;					
		*max_F = opt_hmean_over;
		*theta_best = theta_best_over;		
	}
	else
	{
		*pos_halfplane = -1;				
		*max_F = opt_hmean_under;		
		*theta_best = theta_best_under;	
	}
}

/********************************************************************************
    Function to compute the intercept of the line which maximize the Fscore among the lines whose slope is tan(theta_best) and the intercepts
	 those contained in the vector c_values. The lines are scanned according to their intercepts in an increasing order; and for each line 
	 the initial value of true positive(tp), false positive(fp) and false negative (fn) is updated according to label("labels" vector) of the
	  point associated with the current line. 
INPUT
  size: 				number of lines to be analyzed
  tp:					initial value of true positive
  fp:					initial value of false positive
  tn:					initial value of true negative
  fn:					initial value of false negative
  labels:			vector containing the labels of the points associated with the lines
  c_values:			vector containing the intercepts of the lines to be scanned
  order_c_values:	vector containing the intercepts sorted in increasing orderof the lines to be scanned
  max_F:				output variable containing the optimal Fscore value
  c_best:			output variable containing the intercept corresponding to the optimal Fscore value, double 
  theta_best:		the optimal angle previously computed
*/ 
void compute_best_c(int size, int tp, int fp, int tn, int fn, int *labels, double *c_values,
	                          int *order_c_values, double *max_F, double *c_best, double theta_best){
	register int i, h;
	int cnt = 0, pos_labels = 0, neg_labels = 0;
	double tmp_hmean_under;
	*max_F = compute_F(tp, fn, fp);	
	for(i = 0; i < size; i++)
	{		
		cnt = 0; pos_labels = 0; neg_labels = 0;
		// counting the number of collinear points
		while(c_values[i] == c_values[i + cnt + 1] && (i + cnt + 1) < size)cnt++;
		for(h = 0; h <= cnt; h++)
		{			
			if(labels[ order_c_values[i + h] ] > 0)
				pos_labels++;
			else
				neg_labels++;				
 		} 	
 		// updating fscore	 		
		tp += pos_labels; fn -= pos_labels; fp += neg_labels; tn -= neg_labels;	 			
 		// compute the F-score relative to the current line when the positive half-plane is that under the line
		tmp_hmean_under = compute_F(tp, fn, fp);		
		// check whether current hmean is greater than actual maximum Fscore
		if(tmp_hmean_under > *max_F){					
			*max_F = tmp_hmean_under;
			*c_best = c_values[i];				
		}		
		i = i + cnt;
	}
}


/********************************************************************************
 	This function compute the line among those crossing each labeled point (pos_vect[i], neg_vect[i]) which correspond to the highest Fscore.
 	 The points labels are contained in the vector "labels. The algorithm is made up by two steps:
 			Step 1) compute the line (hence its angle alpha) which maximizes the Fscore among those crossing the origin and 
 			 		the points (pos_vect[i], neg_vect[i])  
 			Step 2) compute the line (hence its intercept) which maximizes the Fscore among those having slope tan(alpha) and crossing
 			 		the points (pos_vect[i], neg_vect[i])  
	During the optinization, for each line two possibility are considered: the positive half-plane is above or below the line. The variable
	 "positive_halfplane" will contain the choice which correspond to the highest Fscore.
	 
INPUT    
	pos_vect: 				vector in which position "i" contains the weighted sum of positive neighbors of item "i"
	neg_vect: 				vector in which position "i" contains the weighted sum of negative neighbors of item "i"
	labels:	 				vector containg the item labels 
	n: 						number of points		 
	theta:					output variable which will contain the angle formed with the x axis by the line corresponding to the highest Fscore
	c: 						output variable which will contain the optimal threshold c = -q*cos(theta), where q is the intercept corresponding 
									to the optimal line
   opt_hmean: 				output variable containing the Fscore corresponding to the optimum line 
   positive_halfplane: 	output variable containing the position of the positive half-plane: 1 over, -1 under the line
*/
void error_minimization(double *pos_vect, double *neg_vect, int *labels, int *n, 
           double *theta, double *c, double *opt_hmean, int *positive_halfplane)
{	
	int N_pos = 0, pos_halfplane = 0, N_neg, tp_o, fn_o, fp_o, tn_o, tp_u, fn_u, fp_u, tn_u, 
		cnt, pos_labels, neg_labels, opt_fp_o = 0, opt_fp_u = 0, opt_fn_o=0, opt_fn_u=0, opt_tp_o=0, opt_tp_u=0;
	register int i, h;
	const int n_ = (*n);
	int order_thetas[n_], order_c_values[n_];	
	double max_F = 0, theta_best = 0, c_best = 0, opt_hmean_over = 0, opt_hmean_under = 0;	
	double theta_best_over, theta_best_under=0,tmp_hmean_under, tmp_hmean_over;		 
	double thetas[n_], c_values[n_];		
	// finding the number of positive labels
	N_pos = count_positives(labels, n_);	
	N_neg = n_ - N_pos;
	// initial errors when positive halfplane 'over' the line
	tp_o = N_pos; fp_o = N_neg; tn_o = 0; fn_o = 0;	
	tmp_hmean_over = compute_F(tp_o, fn_o, fp_o);	
	// initial errors when positive halfplane 'under' the line
	tp_u = 0; fp_u = 0; tn_u = N_neg; fn_u = N_pos; tmp_hmean_under = 0.0;
	// computing the angles of each line passing through the origin and a point of the training set	
	compute_angles(pos_vect, neg_vect, order_thetas, thetas, n_);				
	// sorting angles and their indices
	quicksort(thetas, order_thetas, 0, (n_)-1);	
	// scanning ordered angles to find the optimum line	
	for(i = 0; i < n_; i++)
	{		
		if(thetas[i] >= 1.57)break;
		cnt = 0; pos_labels = 0; neg_labels = 0;
		// counting the number of collinear points
		while(thetas[i] == thetas[i + cnt + 1]) cnt++;
		if(i != (n_-1)){
			for(h = 0; h <= cnt; h++)
			{
				if(labels[ order_thetas[i + h] ] > 0)
					pos_labels++;
				else
					neg_labels++;	 
 			} 		
 		}
 		// updating actual errors
 		tp_o -= pos_labels; fn_o += pos_labels;	fp_o -= neg_labels; tn_o += neg_labels;
		tp_u += pos_labels; fn_u -= pos_labels;	fp_u += neg_labels; tn_u -= neg_labels; 		
 		tmp_hmean_over = compute_F(tp_o, fn_o, fp_o); 
		tmp_hmean_under = compute_F(tp_u, fn_u, fp_u);			
		// check whether current F-scores is greater than actual maximum Fscores
		check_update(tmp_hmean_under, &opt_hmean_under, &theta_best_under, thetas[i], 
									&opt_tp_u, &opt_fp_u, &opt_fn_u, tp_u, fp_u, fn_u);
		check_update(tmp_hmean_over, &opt_hmean_over, &theta_best_over, thetas[i], 
									&opt_tp_o, &opt_fp_o, &opt_fn_o, tp_o, fp_o, fn_o);		
 		// increment in order to avoid to consider again collinear points
 		i = i + cnt;
	}	
	// choosing the optimum half-plane, fscore and angle
	check_halfplane(opt_hmean_over, opt_hmean_under, theta_best_over, theta_best_under, &pos_halfplane, &max_F, &theta_best);
								
// ------- Step 2: computing best intercept---------------		
	compute_c(pos_vect, neg_vect, order_c_values, c_values, theta_best, n_);	
	// sorting intercepts and their indices	
	quicksort(c_values, order_c_values, 0, n_-1);
	tp_u = 0; fp_u = 0; tn_u = N_neg; fn_u = N_pos;		
	compute_best_c(*n, tp_u, fp_u, tn_u, fn_u, labels, c_values, order_c_values, &max_F, &c_best, theta_best);			
		
	theta_best = theta_best + DBL_MIN;		
	*opt_hmean = max_F;	
	*theta = theta_best;
	*c = -c_best * cos(*theta);
	*positive_halfplane = (int)pos_halfplane;
}	



static R_NativePrimitiveArgType errorm_t[8] = {REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType globalup_t[7] = {REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};

/* Registering functions to be used into R */
static const R_CMethodDef cMethods[] = {
  {"error_minimization", (DL_FUNC) &error_minimization, 8, errorm_t},
  {"global_update_under", (DL_FUNC) &global_update_under, 7, globalup_t},
//   {"error_minimization", (DL_FUNC) &error_minimization, 8},
//   {"global_update_under", (DL_FUNC) &global_update_under, 7},
   {NULL, NULL, 0}
};


void
R_init_COSNet(DllInfo *info)
{
   R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
