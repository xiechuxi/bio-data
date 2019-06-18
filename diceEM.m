(* ::Package:: *)

 (* This defines the accessors for the parameter set of the dice model. It is stored as a list
    {type1Prob, type2Prob, faceProbs1, faceProbs2}. The accessors make the code more readable. *)
  type1Prob[parameters_] := parameters[[1]];
  type2Prob[parameters_] := parameters[[2]];
  faceProbs1[parameters_] := parameters[[3]];
  faceProbs2[parameters_] := parameters[[4]];
 
 (* diceEM just sets up for the actual EM algorithm by counting the frequencies of the faces
    in each trial and providing initial random or arbitrary values to diceEMIterator.*)
 (* A trial is the result of drawing a die at random from the bag and rolling it n times. *)
 diceEM[trials_, maxIterations_, accuracy_, randomSeed_:314, debug_:False]:=
	Module[{numFaces, binCountsList, initialFaceProbs1, initialFaceProbs2,rad1,rad2},
		(*SeedRandom[] initializes the random number generator so you can get deterministic
		  results for testing purposes.*)
		SeedRandom[randomSeed];
		(* For strictly defined EM, which uses maximum likelihood estimations, any die faces
		   that do not occur anywhere in the input can be treated as though they don't exist.*)
		numFaces = Max[trials[[1]]];
		(*initialize binCountsList here, by mapping the builtin function BinCounts over the trials.
		  Each entry in the resulting binCountsList should be a list of the face frequencies from one trial.
		  E.g., if there are four faces on the dice, one entry might be {3, 0, 1, 1}, indicating 
		  that on this trial face 1 was rolled 3 times, face 2 was not rolled at all, etc.  *)
		binCountsList = Map[BinCounts[#,{1,Max[#]+1}]&,trials];
	    (* Initialize initialFaceProbs1 and initialFaceProbs2 here. Use RandomReal and make it so
	       so that all entries are within a factor of two of one another. The built in functions
	       Normalize and Total may also be useful here. *)
	    (* Put your code here.*)
	    rad1=Table[RandomReal[],numFaces];
	    rad2=Table[RandomReal[],numFaces];
	    initialFaceProbs1 = rad1/Total[rad1];
	    initialFaceProbs2 = rad2/Total[rad2];
   diceEMIterator[binCountsList, 
		           numFaces, 
		           {0.45, 0.55, initialFaceProbs1, initialFaceProbs2}, 
		           maxIterations, 
		           accuracy]]
(* diceEMIterator implements the outer loop of the EM algorithm.
   It calls updateProbs on each iteration. *)		           

diceEMIterator[binCountsList_, numFaces_, initParams_, maxIterations_, accuracy_, debug_:False]:=
	Module[{oldParamEstimates, newParamEstimates,d},
		(* Initialize the local variables. *)
		oldParamEstimates = initParams;
		(* Loop here until either maxIterations has been reached or the accuracy goal has been
		   met. The accuracy goal is met when the sum, over all estimated parameters, of the absolute 
		   values of the changes from one iteration to the next is less than accuracy. *)
		(*I suggest using Do to iterate here.
		   On each iteration, call updateProbs, passing in the old values, to set the new values. 
		   Then test whether the termination conditions have been met (Break[] breaks a loop).
		   Finally, if termination conditions have not been met, set old values to be the same 
		   as the new values.
		   *)
	    (* Put your code here.*)
	    Do[newParamEstimates=updateProbs[binCountsList, oldParamEstimates];
	    If[Total[Total[Abs[newParamEstimates-oldParamEstimates]]]<=accuracy,
	    Break[],oldParamEstimates=newParamEstimates],
	    {d,1,maxIterations}];
		(*At the end, return the estimated parameters with the less likely die first.*)
		If[type1Prob[newParamEstimates] <= type2Prob[newParamEstimates],
		   newParamEstimates,
		   {type2Prob[newParamEstimates], 
		    type1Prob[newParamEstimates], 
		    faceProbs2[newParamEstimates], 
		    faceProbs1[newParamEstimates]}]
	 ]
   
(* updateProbs does the actual EM calculations. Implementing this is your main task. *)
updateProbs[binCountsList_, oldParamEstimates_, debug_:False] :=
	Module[{posteriors,
		    (* type1Count and type2Count are the expected number of times a type1 or type2
		       die was drawn.*) 
		    type1Count, type2Count,
		    (* faceCounts1 is the expected number of times each face was rolled on a die 
		       of type 1.Likewise for faceCounts2.*) 
		    faceCounts1, faceCounts2,i,j,a,b,c,fc1,fc2},
		(*Create list of posterior probabilities of a Type1 die having been rolled on each draw 
		   by calling your dicePosteriors, which you should paste in to this file. *)
		(* Put your code here.*)
		posteriors=Table[dicePosterior[binCountsList[[a]],type1Prob[oldParamEstimates],
		type2Prob[oldParamEstimates],faceProbs1[oldParamEstimates],faceProbs2[oldParamEstimates]],{a,1,Length[binCountsList]}];
		
		(* Now use the posteriors to calculate EXPECTED number of times each die type was drawn. *) 
		(* Put your code here.*)
		type1Count= Total[Total[Table[binCountsList[[i]]*posteriors[[i]],{i,1,Length[binCountsList]}]]];
		type2Count= Total[Total[Table[binCountsList[[j]]*(1-posteriors[[j]]),{j,1,Length[binCountsList]}]]];
		(* Now use the posteriors to calculate EXPECTED number of times each face was rolled
		   on each die typep. *) 
		(* Put your code here.*)
		fc1=Table[binCountsList[[i]]*posteriors[[i]],{i,1,Length[binCountsList]}];
		fc2=Table[binCountsList[[j]]*(1-posteriors[[j]]),{j,1,Length[binCountsList]}];
		faceCounts1= Total[fc1];
		faceCounts2= Total[fc2]; 
		
				(* Finally, use these counts to compute maximum likelihood estimates for the parameters and 
		   return these estimates in a list: {newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2} *)
		(* Put your code here.*)
		{type1Count/(type1Count+type2Count),type2Count/(type1Count+type2Count),faceCounts1/Total[faceCounts1],faceCounts2/Total[faceCounts2]}
	]

(* Make sure to include your diceSample and dicePosterior functions here. If yours don't work get working
   copies from the TAs. *)
diceSample[numType1_, numType2_, type1_, type2_, draws_, rollsPerDraw_] :=
Module[{numType1d=EmpiricalDistribution[type1->Range[Length[type1]]](*distribution of type1 dice1*),
numType2d=EmpiricalDistribution[type2->Range[Length[type2]]](*distribution of type1 dice2*),
h},
Table[If [1==RandomVariate[BernoulliDistribution[numType1/(numType1+numType2)]],
RandomVariate[numType1d,rollsPerDraw],
RandomVariate[numType2d,rollsPerDraw]],{h,1,draws}]]

dicePosterior[binCounts_, type1Prior_, type2Prior_, faceProbs1_, faceProbs2_] := 
Module[{(*replace all zeros^zeros with one*)
P1=Fold[Times,Diagonal[Table[faceProbs1^binCounts[[m]],{m,1,Length[binCounts]}]/.Indeterminate-> 1]], 
P2=Fold[Times,Diagonal[Table[faceProbs2^binCounts[[k]],{k,1,Length[binCounts]}]/.Indeterminate-> 1]]},
type1Prior*P1/(type1Prior*P1+type2Prior*P2)];
(*myRound is a hack to get around a problem with rounding numbers in \
Mathematica. Please don't bother to try to understand it.*)

myRound[x_, n_] :=
  N[IntegerPart[Round[x,10^-n]*10^n] / 10^n];
