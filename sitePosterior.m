(* ::Package:: *)

(*A draw is the observed bases of a drawn sequence, which is subsequently replaced.
  sitePosterior calculates the posterior probability of a bound site versus a non-bound site, based on the
  observed sequence in the draw, the proportion of bound vs. non-bound sequences in the bag, and the 
  probabilities of observing each base in a sequence for bound and non-bound sites.
  The single number returned is the posterior probability of a bound site.*)
sitePosterior[sequence_, sitePrior_, backgroundPrior_, siteProbs_, backgroundProbs_] := 
Module[{Psequencesite=Diagonal[Table[siteProbs[[i]][[sequence[[j]]]],{i,1,Length[siteProbs]},{j,1,Length[sequence]}]],
Pbackground=Table[backgroundProbs[[sequence[[k]]]],{k,1,Length[sequence]}]},
If[Fold[Times,1,Psequencesite]==0&Fold[Times,1,Pbackground]==0,Pbackground=Pbackground,Pbackground=ReplaceAll[0->1][Pbackground]];
Fold[Times,1,Psequencesite]*sitePrior/
(Fold[Times,1,Psequencesite]*sitePrior+
Fold[Times,1,Pbackground]*backgroundPrior)
		   ]
