#set AMPLZg; 		# set with the product categories in Group g
#set AMPLEa; 		# set with the shelf segments in Aisle a
#set AMPLEb; 		# set with the indices of the shelves in Aisle a
#set AMPLEb1;		# set with the shelf segments in 'SHELF 400' of Aisle a
#set AMPLEb2;		# set with the shelf segments in 'SHELF 500' of Aisle a
# More Ebx's would be defined based on the number of elements in Eb (i.e. the number of shelves)

#param AMPLALPHAa;			# Smallest shelf-segment index within aisle a
#param AMPLBETAa;			# Largest shelf-segment index within aisle a 
#param AMPLFIRSTSHELFa;		# Smallest shelf index within aisle a (extracted from the set Eb)
#param AMPLSECONDSHELFa; 		# Largest shelf index within aisle a (extracted from the set Eb)
#param AMPLP {p in AMPLZg};		# Avg unit prfit margin for product category p
#param AMPLi {p in AMPLZg};		# Impulse purchase rate for product category p
#param AMPLl {p in AMPLZg};		# Min shelf-space requirement for product category p (as a length)
#param AMPLu {p in AMPLZg};		# Max shelf-space requirement for product category p (as a length)
#param AMPLPHI {p in AMPLZg};	# Smallest product facing length within product category p
#param AMPLc {e in AMPLEa};		# Shelf space available along shelf segment e (as a length)
#param AMPLk {e in AMPLEa};		# Customer traffic density along shelf segment e (based on location)

var Vpe {p in AMPLZg, e in AMPLEa};						# Visibility of product category p over shelf-segment e
var Vp {p in AMPLZg};								# Total visibility of product category p
var Spe {p in AMPLZg, e in AMPLEa} >= 0;				# Space allocated to product category p over shelf-segment e
var Hpe {p in AMPLZg, e in AMPLALPHAa..(AMPLBETAa-1)} >= 0;	# Auxiliary variable for organizing over multiple segments
var Ype {p in AMPLZg, e in AMPLEa} binary;				# 1 iff product category p is assigned to shelf-segment e
var Qpb {p in AMPLZg, b in AMPLEb} binary;				# 1 iff product category p is assigned to shelf b (a product cannot be assigned over multiple shelves)



maximize Obj : (  sum {p in AMPLZg} ( AMPLP[p]*AMPLi[p]*Vp[p] )  );



subject to st20 {p in AMPLZg, e in AMPLEa}: Vpe[p,e] = (AMPLk[e]*Spe[p,e])/AMPLc[e] ;

subject to st21 {p in AMPLZg}: Vp[p] = (sum {e in AMPLEa} Vpe[p,e]) ;


subject to st22a {p in AMPLZg}: (sum {e in AMPLEa} Spe[p,e]) <= AMPLu[p] ;

subject to st22b {p in AMPLZg}: (sum {e in AMPLEa} Spe[p,e]) >= AMPLl[p] ;


subject to st23 {e in AMPLEa}: (sum {p in AMPLZg} Spe[p,e]) <= AMPLc[e] ;


subject to st24a {p in AMPLZg, e in AMPLEa}: Spe[p,e] <= ( min(AMPLc[e],AMPLu[p]) )*Ype[p,e] ;

subject to st24b {p in AMPLZg, e in AMPLEa}: Spe[p,e] >= AMPLPHI[p]*Ype[p,e] ;


# Ensuring that all products are assigned and all the segments are filled up.
# This is needed when using REAL DATA, cz there may be some cases where no profitability/supports are not available for an entire department even.
subject to st24c : (sum {p in AMPLZg, e in AMPLEa} Spe[p,e]) = (sum {e in AMPLEa} AMPLc[e]) ;
subject to st24d {p in AMPLZg}: (sum {e in AMPLEa} Ype[p,e]) >= 1 ;


subject to st25 {p in AMPLZg, e in AMPLEa, f in AMPLEa, j in AMPLEa: e<f<j}: Spe[p,f] >= AMPLc[f]*( Ype[p,e] + Ype[p,j] - 1 ) ;
# 'i' of the original formulation is replaced by 'f' here, to avoid confusing this index with the parameter 'i'

subject to st26 {p in AMPLZg, e in AMPLALPHAa..(AMPLBETAa-1)}: Hpe[p,e] >= ( Ype[p,e+1] + Ype[p,e] - 1 ) ;

subject to st27 {e in AMPLALPHAa..(AMPLBETAa-1)}: (sum {p in AMPLZg} Hpe[p,e]) <= 1 ;

# Might have to have MATLAB direct to multiple AMPL files depending on how many shelves in the considered aisle
# Else  it will be challenging to define the different versions of each constaint: 28a for 1 shelf, 28a/b/c for 3 shelves

subject to st28a {p in AMPLZg, b in AMPLFIRSTSHELFa..AMPLFIRSTSHELFa}: (1 - (sum {e in AMPLEb1} Ype[p,e]) ) <= 0 + 10000*(1-Qpb[p,b]) ;
subject to st28b {p in AMPLZg, b in AMPLSECONDSHELFa..AMPLSECONDSHELFa}: (1 - (sum {e in AMPLEb2} Ype[p,e]) ) <= 0 + 10000*(1-Qpb[p,b]) ;  

subject to st29a {p in AMPLZg, b in AMPLFIRSTSHELFa..AMPLFIRSTSHELFa}: (sum {e in AMPLEb1} Ype[p,e]) <= 0 + 10000*Qpb[p,b] ;
subject to st29b {p in AMPLZg, b in AMPLSECONDSHELFa..AMPLSECONDSHELFa}: (sum {e in AMPLEb2} Ype[p,e]) <= 0 + 10000*Qpb[p,b] ;  

subject to st30 {p in AMPLZg, r in AMPLFIRSTSHELFa..AMPLFIRSTSHELFa, s in AMPLSECONDSHELFa..AMPLSECONDSHELFa}: ( Qpb[p,r] + Qpb[p,s] ) = 1 ;  