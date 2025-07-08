#set AMPLZg; 		# set with the VIRTUAL product categories in Group g (defined as being equal to the number of product categories - some may endup empty)
#set AMPLPRODg; 		# set with the product categories in Group g

#param AMPLP_prod {prod in AMPLPRODg};								# Avg unit prfit margin for product category p
#param AMPLASSOCIATION_prod {prod1 in AMPLPRODg, prod2 in AMPLPRODg};	# Association rate between product categories p and r
#param AMPLN {p in AMPLZg};											# The number of products in each group (if the number of categories is even - all N's =2; if odd - last N=1, other N's=2)

var P {p in AMPLZg};										# Avg unit prfit margin for product category p
var X {p in AMPLZg, prod in AMPLPRODg} binary;					# 1 iff product category 'prod' is assigned to virtual category 'p'
var W {p in AMPLZg, prod1 in AMPLPRODg, prod2 in AMPLPRODg} binary;	# 1 iff product categories 'prod1' & 'prod2' BOTH are assigned to virtual category 'p'
var i_cat {p in AMPLZg};									# Inter-category mpulse purchase rate for product category p


maximize Obj : (  sum {p in AMPLZg} ( P[p]*i_cat[p] )  );


# These may have to written as the first few constraints, as the virtual categories are used there onwards:
subject to st31 {p in AMPLZg}: (sum {prod in AMPLPRODg} X[p,prod]) = AMPLN[p];				# Two product categories may be assigned to a virtual category
subject to st32 {prod in AMPLPRODg}: (sum {p in AMPLZg} X[p,prod]) = 1;					# Every product MUST be assigned to a virtual category	


subject to st33 {p in AMPLZg}: i_cat[p] = 0.5 * ( sum {prod1 in AMPLPRODg, prod2 in AMPLPRODg: prod1<>prod2}  ( (AMPLASSOCIATION_prod[prod1,prod2])*(W[p,prod1,prod2]) ) ) ;	# The purchase impulse of virtual category 'p' is now based on the association between the product categories in it

subject to st34 {p in AMPLZg, prod1 in AMPLPRODg, prod2 in AMPLPRODg: prod1<>prod2}: (W[p,prod1,prod2]) <= 0.5*( X[p,prod1] + X[p,prod2] )   ;	# 1st of 2 constraints ensuring that W is 1, iff both corresponding X's are equal to 1 [see 470 notes or page F9 for details]
subject to st35 {p in AMPLZg, prod1 in AMPLPRODg, prod2 in AMPLPRODg: prod1<>prod2}: (W[p,prod1,prod2]) >= ( X[p,prod1] + X[p,prod2] - 1 )   ;	# 2nd of 2 constraints ensuring that W is 1, iff both corresponding X's are equal to 1 [see 470 notes or page F9 for details]


subject to st36 {p in AMPLZg}: P[p]=(sum {prod in AMPLPRODg} X[p,prod]*AMPLP_prod[prod]);		# Adding the profit margins of the product(s) assigned to each virtual category

