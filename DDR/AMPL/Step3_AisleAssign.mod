#set AMPLAISLES;			# set with the indices of store AMPLAISLES
#set AMPLGROUPS;			# set with the indices of AMPLGROUPS of Product Categories
#set AMPLFEASIBLECASES within {AMPLGROUPS,AMPLAISLES}; 	# set containing the (g,a) cases which have been deemed 'feasible for assignment' based on our criteria [done on MATLAB]
#set AMPLINFEASIBLECASES within {AMPLGROUPS,AMPLAISLES};	# set containing the (g,a) cases which have been deemed 'INfeasible for assignment' based on our criteria [done on MATLAB]

#param AMPLPSI {g in AMPLGROUPS, a in AMPLAISLES};		# Optimal profits achievable (BASED ON SSAP(g,a) by assigning category group 'g' to aisle 'a')
#param AMPLi {g in AMPLGROUPS, a in AMPLAISLES};			# The purchase impulse of assigning group 'g' to aisle 'a' (this is based on the group located in the same aisle, the previous time peiod)

var T {g in AMPLGROUPS, a in AMPLAISLES} binary;	# 1 iff product category group 'g' is assigned to aisle 'a'



maximize Obj : (  sum {a in AMPLAISLES, g in AMPLGROUPS: (g,a) in AMPLFEASIBLECASES} ( AMPLPSI[g,a]*T[g,a]*AMPLi[g,a] )  );



subject to st31 {g in AMPLGROUPS}: (sum {a in AMPLAISLES: (g,a) in AMPLFEASIBLECASES} T[g,a]) = 1 ;

subject to st32 {a in AMPLAISLES}: (sum {g in AMPLGROUPS: (g,a) in AMPLFEASIBLECASES} T[g,a]) = 1 ;

# Setting 0 for infeasible T's (this is just for our convenience, no need to show)
subject to st33 {a in AMPLAISLES, g in AMPLGROUPS: (g,a) in AMPLINFEASIBLECASES}: T[g,a] = 0 ;
