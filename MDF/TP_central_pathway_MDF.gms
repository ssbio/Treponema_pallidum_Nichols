*************************************************************
********MAX/MIN DRIVING FORCE ANALYSIS FOR TREPONEMA*********
*************GLYCOLYSYS BIOSYNTHESIS PATHWAY*****************
********************NIAZ BAHAR CHOWDHURY*********************
*************************************************************

$INLINECOM /*  */
$ONEMPTY

OPTIONS

	limrow = 1000
       	optCR = 0
        optCA = 0
        iterlim = 100000
        decimals = 7
        reslim = 100000
        work = 5000000;

*********Defining Sets**************************************
SETS

	i					set of metabolites

$include "metabolites.txt"	

	j					set of reactions

$include "reactions.txt"
;
*************************************************************

***********Defining Parameters*******************************
PARAMETERS

	S(i,j)					stoichiometric matrix

$include "sij.txt"

	delta_G_o(j)				maximum flux of v(j)

$include "deltaGo.txt"

	Cmax(i)					highest concentration

$include "cmax.txt"

	Cmin(i)					Lowest concentration

$include "cmin.txt"

	R					Gas constant
/0.008314/

	T					Temperature
/310.15/
;
**************************************************************

*********Defining Equations***********************************
EQUATIONS

	objective				objective function
	constraint_1(j)				constranit 1 of the formulation
	constraint_2(j)				constranit 2 of the formulation
	upper_bound(i)				upper bound
	lower_bound(i)				lower bound
	ATPtoADP				ATP to ADP ratio fixed according to Noor et al. 2014
	NADHtoNAD				NADH to NAD ratio fixed according to Noor et al. 2014
;
**************************************************************

*********Defining Variables***********************************
POSITIVE VARIABLES

	x(i)					concentration of a metabolite

FREE VARIABLES

	v(j)					reaction flux
	B					Maximum driving force
	deltaG(j)				driving force in the given condition
	Z					objective value
;
****************************************************************

***************Defining Model***********************************
objective..			Z =e= B;

constraint_1(j)..		- deltaG(j) =g= B;

constraint_2(j)..		deltaG(j) =e= delta_G_o(j) + R * T * sum(i, S(i,j) * x(i));

upper_bound(i)..		x(i) =l= log(Cmax(i));

lower_bound(i)..		x(i) =g= log(Cmin(i));

ATPtoADP..			x('C00002') =e= 10 * x('C00008');

NADHtoNAD..			x('C00004') =e= 0.1 * x('C00003');

Model MDF /all/;
******************************************************************

**********Solving Model*********************

solve MDF using lp maximizing Z;

********************************************
****************Output File*****************
FILE RESULTS1 /DELTAG.txt/;

PUT RESULTS1;

PUT "reaction      deltaG"/;

LOOP(j,
	
	PUT j.tl:0:100,"    ", deltaG.l(j):20:5/;
		
);

PUTCLOSE;

FILE RESULTS2 /CONCENTRATION.txt/;

PUT RESULTS2;

PUT "metabolites      concentration"/;

LOOP(i,
	
	PUT i.tl:0:100,"    ", x.l(i):20:8/;
		
);

PUTCLOSE;
**********************************************