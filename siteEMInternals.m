(* ::Package:: *)

(* Put your code for sequencePosteriors, updateMotifPrior, and updatePFMCounts,
   along with any other functions you write to implement them, in this file.
   Don't forget that you can use your implementation of sitePosterior as part of your
   implementation of sequencePosteriors or, if you prefer, you can request a correct
   version from the TAs. *)
sitePosterior[sequence_, sitePrior_, backgroundPrior_, siteProbs_, backgroundProbs_] := 
Module[{Psequencesite=Diagonal[Table[siteProbs[[i]][[sequence[[j]]]],{i,1,Length[siteProbs]},{j,1,Length[sequence]}]],
Pbackground=Table[backgroundProbs[[sequence[[k]]]],{k,1,Length[sequence]}]},
If[Fold[Times,1,Psequencesite]==0&Fold[Times,1,Pbackground]==0,Pbackground=Pbackground,Pbackground=ReplaceAll[0->1][Pbackground]];
Fold[Times,1,Psequencesite]*sitePrior/
(Fold[Times,1,Psequencesite]*sitePrior+
Fold[Times,1,Pbackground]*backgroundPrior)
		   ];

sequencePosteriors[inputsequence_,oldMotifPrior_, oldPFM_, backgroundFreqs_]:= 
Module[{i,W=Length[oldPFM],L=Length[inputsequence],backgroundPrior},
backgroundPrior=(1-oldMotifPrior);
Table[
sitePosterior[Take[inputsequence,{i,i+W-1}],oldMotifPrior,backgroundPrior,oldPFM,backgroundFreqs],{i,1,L-W+1} ]];

updateMotifPrior[normalizedPosteriors_]:= 
Module[{i,
N=Length[normalizedPosteriors]},
Total[Table[normalizedPosteriors[[i]]/Length[normalizedPosteriors[[i]]],{i,1,N}],{1,N}]/N];

updatePFMCounts[motifLength_, input_, normalizedPosteriors_, motifPseudocounts_, erasers_]:= 
Module[{j,k,m,i,erasers1},
If[TrueQ[erasers],erasers1=erasers,erasers1=Map[ConstantArray[1.0, #] &, Map[Length, input]]];
Table[
(Total[Table[If[input[[j]][[k]]==i,1,0]*erasers1[[j]][[k]],{j,1,Length[input]},{k,m,m+Length[input[[j]]]-motifLength}]*
normalizedPosteriors,{1,Length[normalizedPosteriors]}]+motifPseudocounts[[i]])/(Total[normalizedPosteriors,{1,Length[normalizedPosteriors]}]
+Total[motifPseudocounts]),{m,1,motifLength},{i,1,4}]
]

