(* ::Package:: *)

 <<tools`
 
createPHMM[msaFile_]:=
	Module[{sequences, matches, matchStates, numMatch, pHMM, i},
		(*read in the msaFile*)
		sequences = parsingMSA[msaFile];
		(*determine which positions in the sequences are match states*)
		matches = Select[Map[Select[#, StringQ] &, Transpose[sequences]], Length[#] < Length[sequences]/2 &];
		matchStates = Sort[DeleteDuplicates[Flatten[Map[Position[Map[Select[#, StringQ] &, Transpose[sequences]], #] &, matches]]]];
		numMatch = Length[matchStates];
		(*get skeleton PHMM, assuming alphabet is DNA (4 letters)*)
		pHMM=skeletonPHMM[numMatch, 4];
		(*loop through all sequences, updating the PHMM counts*)
		For[i=1, i<=Length[sequences], i++,
			pHMM=updatePHMM[sequences[[i]], matchStates, pHMM]
		];
		(*normalize the PHMM counts to return PHMM probabilities*)
		normalizePHMM[pHMM]
	]

(*
Reminder about the PHMM structure:
Insert0 state emission probabilities (list length = length of alphabet)
Transition probabilities for B -> M1, B -> I0, B -> D1; I0 -> M1, I0 -> I0
A 3 D emission array, 
	where the first level corresponds to the index of the states (so if we have 10 match states, we expect length of 10 for the emission array), 
	the second level corresponds to the state type (Insert, or Match, so we expect length of 2 for the emission subarrays), 
	and the last level corresponds to the alphabet emission probabilities (so for the DNA alphabet, we expect length of 4 for the emission sub - subarray)
A 2 D transition array, 
	where the first level corresponds to the index of the the states (so if we have 10 match states, we expect length of 10 for the transition array), 
	and the second level corresponds to the transition probabilities to the next states 
		(so we expect 7 values for each subarray in the transition array corresponding to 
		Mk -> Mk + 1, Mk -> Ik, Mk -> Dk + 1; Ik -> Mk + 1, Ik -> Ik; Dk -> Mk + 1, Dk -> Dk + 1)

update PHMM pseudocounts based on given sequence
*)
updatePHMM[sequence_, matchStates_, pHMM_]:=
	Module[{pHMMIndex, currentState, sequenceIndex,
		initialInsertEmissions, initialTransitions, emissions, transitions,e1,t1,countI0,j,k},
		
		{initialInsertEmissions, initialTransitions, emissions, transitions}=pHMM;
		
		countI0=0;
		Do[If[MemberQ[{1,2,3,4},sequence[[i]]]&&MemberQ[{1,2,3,4},sequence[[i+1]]],countI0++,countI0=countI0],{i,1,matchStates[[1]]-2}];
		
		initialInsertEmissions=initialInsertEmissions+BinCounts[Take[sequence,matchStates[[1]]-1],{1,5}];
		
		initialTransitions=initialTransitions+{If[matchStates[[1]]==1,1,0],If[matchStates[[1]]!=1 && sequence[[1]]!="-",1,0],
		If[sequence[[matchStates[[1]]]]=="-",1,0],If[matchStates[[1]]!=1 && sequence[[matchStates[[1]]-1]]!="-",1,0],countI0};
		
		(*trim the sequence with "-" not at matchstates position*)
		sequence=Select[sequence,#!="-"||MemberQ[matchStates,Position[sequence,#]]&];
		
		(*modify matchstates to take out non delete states "-"*)
	    Do[If[sequence[[Length[sequence]+1-k]]=="-"&&!MemberQ[matchStates,Length[sequence]+1-k],
	    matchStates=Join[TakeWhile[matchStates,#<Length[sequence]+1-k&],Select[matchStates,#>Length[sequence]+1-k&]-1],
	    matchStates=matchStates],{k,1,Length[sequence]}];
		
		emissions=emissions+Table[{BinCounts[sequence[[matchStates[[e1]]]],{1,5}],
		If[matchStates[[e1+1]]-matchStates[[e1]]==1,{0,0,0,0},
		BinCounts[Take[sequence,{matchStates[[e1]]+1,matchStates[[e1+1]]-1}],{1,5}]]},{e1,1,Length[matchStates]}];
		
	    
		transitions=transitions+Table[{If[sequence[[matchStates[[t1]]]]!="-"&&sequence[[matchStates[[t1+1]]]]!="-"&&matchStates[[t1]]+1==matchStates[[t1+1]],1,0],
		If[sequence[[matchStates[[t1]]]]!="-"&&matchStates[[t1]]+1!=matchStates[[t1+1]],1,0],
		If[sequence[[matchStates[[t1]]]]!="-"&&sequence[[matchStates[[t1]]+1]]==sequence[[matchStates[[t1+1]]]]=="-",1,0],
		If[sequence[[matchStates[[t1+1]]]]!="-"&&matchStates[[t1+1]]-1!=matchStates[[t1]],1,0],
		If[matchStates[[t1]]+1!=matchStates[[t1+1]],matchStates[[t1+1]]-matchStates[[t1]]-2,0],
		If[sequence[[matchStates[[t1]]]]=="-"&&matchStates[[t1]]+1==matchStates[[t1+1]],1,0],
		If[sequence[[matchStates[[t1]]]]=="-"&&sequence[[matchStates[[t1+1]]]]=="-",1,0]}
		,{t1,1,Length[matchStates]}];
		
		(*
		Walk through given sequence.
		Note that transitions up to the first match state will update the initial transitions list of the PHMM
		All following transitions will update the 2D transitions array
		*)
	
		(*Update the transitions into the end state (Mk+1)*)
		(*return updated pseudocounts*)
		{initialInsertEmissions, initialTransitions, emissions, transitions}
	]
	
	
(*Normalize the parts of the PHMM to get probability values instead of counts
Keep in mind that probabilities should not add up to one for every sublist, but rather for:
	transitions out of the same state
	emissions from the same state
*)
normalizePHMM[{initialInsertEmissions_, initialTransitions_, emissions_, transitions_}]:=
	Module[{normalizedInitialInsertEmissions, normalizedInitialTransitions, normalizedEmissions, normalizedTransitions,e1,e2,t1},
	
	normalizedInitialInsertEmissions=initialInsertEmissions/Total[initialInsertEmissions];
	
	normalizedInitialTransitions=Join[Take[initialTransitions,{1,3}]/Total[Take[initialTransitions,{1,3}]],
	Take[initialTransitions,{4,5}]/Total[Take[initialTransitions,{4,5}]]];
	
	normalizedEmissions=Table[emissions[[e1,e2]]/Total[emissions[[e1,e2]]],{e1,1,Length[emissions]},{e2,1,2}];
	
	normalizedTransitions=Table[Join[Take[transitions[[t1]],{1,3}]/Total[Take[transitions[[t1]],{1,3}]]
	,Take[transitions[[t1]],{4,5}]/Total[Take[transitions[[t1]],{4,5}]]
	,Take[transitions[[t1]],{6,7}]/Total[Take[transitions[[t1]],{6,7}]]],{t1,1,Length[transitions]}];
	
	
		(*Return the normalized PHMM*)
		{normalizedInitialInsertEmissions, normalizedInitialTransitions, normalizedEmissions, normalizedTransitions}
	]
	

(*extra credit parser to stockholm file*)
parser[msaFile_]:=
	Module[{input, seqs},
	input=Import["sequence.txt"];
	seqs=input/.{"-"->"."}]
