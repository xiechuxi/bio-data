(* ::Package:: *)

(* Mathematica Raw Program *)

generateRandomParams[regulatorsMask_, seed_]:=
Module[{L=Length[regulatorsMask],i,j},
SeedRandom[seed];
{Table[If[regulatorsMask[[k]]==0,0.0,RandomReal[{0.1,1}]],{k,1,L}],
(*This generate a vector*)
Table[If[regulatorsMask[[j]]==0,
ConstantArray[0.0,L],
Table[If[i==j,0.0,RandomReal[{0.1,1}]],{i,1,L}]],
{j,1,L}]
}
	(*This generate a matrix*)
]
(* genericEqns generates system of equations for n mutually repressing proteins. They are
   generic in that the parameters remain as indexed symbols b[j] and c[j,i], except that
   c[i,i] is replaced by zero, implementing the ban on auto-repression.*)
genericEqns[n_]:=
Module[{i,j},
Table[b[i]*m[i]*Fold[Times,Table[m[j]*If[j==i,0,c[j,i]]+1,{j,1,n}]]-1==0,{i,1,n}]
]
(* mutantEqns takes a genePresenceMask indicating which genes were deleted in the cells 
   in which these mRNA measurements were taken. It generates a set of generic equations 
   and then replaces the mRNA level variables for the deleted genes with zeros, ensuring 
   that any actual measurments made on those deleted genes are ignored and any solver 
   doesn't have to deal with them.
*)
mutantEqns[genePresenceMask_]:=
Module[{i,j,n,eqns},
n=Length[genePresenceMask];
eqns=Table[b[i]*m[i]*Fold[Times,
Table[If[genePresenceMask[[j]]==0,0,m[j]]
*If[j==i,0,c[j,i]]+1,{j,1,n}]]-1==0,{i,1,n}];
Delete[eqns,Position[genePresenceMask,0]]
]
(* mutantExpVars returns a list of variables representing the expression levels of
   genes that remain in the genotype described by genePresenceMask. *)
mutantExpVars[genePresenceMask_]:=
Module[{i,originalexp},
originalexp=Table[m[i],{i,1,Length[genePresenceMask]}];
Delete[originalexp,Position[genePresenceMask,0]]
]
(* expressionEqns replaces the symbolic parameters of the generic equations with actual values
   from a parameter matrix, returning a final set of equations with fixed parameters but variables
   for mRNA levels. Solving these for the mRNA variables yields a set of expression levels that
   are consistent with the parameters.*)
expressionEqns[eqns_, {bVector_, cMatrix_}]:=
Module[{i,eb,j,k},
eb=eqns/.b[i_]:>bVector[[i]];
eb/.c[j_,k_]:>cMatrix[[j,k]]
]
(* expressionProfile finds a set of steady state mRNA expression levels that are consistent with
   a set of expression equations and a genotype as indicated by genePresenceMask. To do this, you
   will use the built in function FindRoot to find a solution to the set of simultaneous equations
   you have constructed with the function expressionEqns. IMPORTANT: FindRoot allows you to specify
   a starting, minimum, and maximum value for each unknown (in this case the mRNA levels). You need
   to do this, making use of the bounds specified in the Introduction section of the assignment notebook. 
   
   The result of FindRoot is a list of replacement rules for expression variables. These are then used 
   to construct a vector of length n containing the expression levels for genes that are present and 
   zeros for genes that are deleted.*)
expressionProfile[expressionEqns_, genePresenceMask_]:=
Module[{modifiedexpressionEqns,i,vartable},
modifiedexpressionEqns=Delete[expressionEqns,Position[genePresenceMask,0]];
vartable=Table[{m[i],1,1,100},{i,1,Length[expressionEqns]}];
FindRoot[modifiedexpressionEqns,Delete[vartable,Position[genePresenceMask,0]]]
]
(* expressionMatrix takes parameters and a list of genePresenceMasks, one for each strain that we
   want an expression profile for. For each mask, it makes a set of equations and calls
   expressionProfile with those equations and the mask. This returns a list of rules that is
   applied to the expression level variables to produce actual values, which are returned.
   The return value is a matrix in which each row represents the expression levels for one
   strain (as described by one genePresenceMask). *)
expressionMatrix[params_, genePresenceMasks_]:=
	(*This ensures that b, c, m are undefined here and in all functions called from here.*)
	Block[{b, c,i,eqns,j,orgmatrix,m},
	eqns=expressionEqns[genericEqns[Length[params[[1]]]],params];
	(*generate expressionequations*)
	orgmatrix=Table[
	Table[m[i],{i,1,Length[params[[1]]]}]
	/.expressionProfile[eqns,genePresenceMasks[[j]]],
	{j,1,Length[genePresenceMasks]}];
	orgmatrix*genePresenceMasks
	]
(*Part 4: Parameter inference.*)
	
(* inferParameters takes a list of expression profiles and a 
   matched list of genePresenceMasks, one for each strain that we have an expression 
   profile for. For each profile-mask pair, it makes a set of equations with the c parameters
   as unknown variables. It is important that the the right-hand-side of each equation be zero.
   A vector of the left-hand-sides of these equations is then dotted with itself to make a
   sum of squared errors, where the error is the deviation from zero. NMminize is then used to 
   find values of the c[i,j] parameters that are consistent with the equations, or as close to
   consistent as possible. We must use NMiminize instead of NSolve because the number of 
   equations may be larger than the number of variables and NSolve won't accept over-determined 
   system of equations *)

(* If the level of regulator j is 0 in all experiments then c[j,i] will never appear in the equations 
   and therefore can't be included in the list of variables to solve for.*)

inferParameters[regulatorsMask_, geneExpressionProfiles_, genePresenceMasks_]:=
	(*Make sure b, c, m are undefined here and in all functions called from here.*)
	Block[{b, c, m,eqns1,i,j,k,solveqns,parameter,parameter1,replacerule,parameter2},
	eqns1=Table[mutantEqnsmod[genePresenceMasks[[i]]],{i,1,Length[genePresenceMasks]}];
	eqns1=Table[parameterEqns[eqns1[[j]],geneExpressionProfiles[[j]]],{j,1,Length[geneExpressionProfiles]}];
	eqns1=Transpose[eqns1];
	solveqns=Table[parameterSumSquaredError[eqns1[[k]]],{k,1,Length[eqns1]}];
	parameter={Table[b[i],{i,1,Length[regulatorsMask]}],
	Table[Table[If[i==j,0,c[j,i]],{i,1,Length[regulatorsMask]}],{j,1,Length[regulatorsMask]}]};
	
	(*modify the parameter so that spcial cases parameters are excluded*)
	parameter1={parameter[[1]]*
	Table[Max[Transpose[genePresenceMasks][[k]]],{k,1,Length[Transpose[genePresenceMasks]]}],
	parameter[[2]]*regulatorsMask*
	Table[Max[Transpose[genePresenceMasks][[k]]],{k,1,Length[Transpose[genePresenceMasks]]}]};
	parameter1=Flatten[parameter1];
	parameter1=DeleteCases[parameter1,0];
	replacerule=NMinimize[solveqns,parameter1];
	parameter2=parameter/.replacerule[[2]];
	parameter2=parameter2/.c[j_,i_]:>0;
	parameter2=parameter2/.b[i_]:>0;
	parameter2
	 ]

(* parameterEqns replaces the symbolic mRNA expression levels of the generic equations with actual values
   from a gene expression matrix, returning a final set of equations with fixed mRNA levels but variables
   for parameters. Solving these for the parameter variables yields a set of parameters that
   are consistent with the expression levels.*)

parameterEqns[eqns_, geneExpressionProfile_]:=
Module[{i},
eqns/.m[i_]:>geneExpressionProfile[[i]]
]


parameterSumSquaredError[eqns_]:=
Module[{rightside,leftside,i},
rightside=Table[eqns[[i]][[1]],{i,1,Length[eqns]}];
leftside=Table[eqns[[i]][[2]],{i,1,Length[eqns]}];
Dot[rightside-leftside,rightside-leftside]
]

(*SECTION Some auxiliary functions to support testing parameter inference.*)
mutantEqnsmod[genePresenceMask_]:=
Module[{i,j,n,eqns},
n=Length[genePresenceMask];
eqns=Table[b[i]*If[genePresenceMask[[i]]==0,0,m[i]]*Fold[Times,
Table[If[genePresenceMask[[j]]==0,0,m[j]]
*If[j==i,0,c[j,i]]+1,{j,1,n}]]-1==0,{i,1,n}]
]
(* Generate a random set of parameters, generate expression profiles for those 
   parameters in wild-type cells and cells with all genes deleted one at a time.
   Then use these expression levels to infer the parameters. Compare the inferred 
   parameters to the ones used to generate the expression levels. If the 
   difference is less than or equal to 10^-4 replace it by zero. Return the resulting matrix of differences. 
   
   It should contain nothing but zeros.*)		 

testParameterInference[n_, seed_]:=
	With[{regulatorsMask=ConstantArray[1, {n}],
		  genePresenceMasks=allSinglesAndWTPresenceMatrix[n]},
		With[{params=generateRandomParams[regulatorsMask, seed]},
			With[{expressionMatrix=expressionMatrix[params, genePresenceMasks]},
				With[{diff=params-inferParameters[regulatorsMask,
					 	                  		  expressionMatrix,
					 	     			          genePresenceMasks]},
				Round[diff, 0.0001]
				]]]]

allSinglesAndWTPresenceMatrix[n_]:=
	Join[ConstantArray[1,{1,n}],
		 ConstantArray[1, {n,n}] - IdentityMatrix[n]]
		 
