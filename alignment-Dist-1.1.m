(* ::Package:: *)

(* 
In this lab, you will implement the Smith-Waterman algorithm for local alignment of
two DNA sequences. 

INPUT
 - The actual inputs to smithWaterman are file paths to two fasta files containing
   DNA sequences and one file containing a scoring matrix. smithWaterman reads these
   files in and sets seq1, seq2, and scoringMatrix to the results.
 - seq1 and seq2 are lists of integers representing DNA sequences, e.g., {2,3,4,1,2}
 - scoringMatrix is a 5x5 symmteric matrix giving the -log2 probabilities of each an 
   alignment column containing each pair consisting of two DNA bases or one base and 
   an insertion symbols, represented by number 5. The score for a pair of two insertion
   symbols should be -Infinity.
OUTPUT
- A matrix whose columns are printed representations of alignment columns. The matrix
  should be displated by passing it to prettyPrintAlignment. 
  *)

smithWaterman[fastaFilePath1_, fastaFilePath2_, scoringMatrixPath_] := 
	Module[{seq1, seq2, scoringMatrix},
		seq1=readFasta[fastaFilePath1];
		seq2=readFasta[fastaFilePath2];
		scoringMatrix=readScoringMatrix[scoringMatrixPath];
		alignmentDisplayMatrix[swTraceback[swBuildMatrix[seq1, seq2, scoringMatrix],
										   seq1, seq2, scoringMatrix]
	]]

(*swBuildMatrix returns a list whose first element is the coordinates of a top-scoring cell in the 
  Smith-Waterman matrix. The second element is the Smith Waterman matrix itself. No traceback
  pointers are kept.*)
swBuildMatrix[seq1_, seq2_, scoringMatrix_] :=
	Module[{matrix, index1, index2, topScoringCell={1, 1},M,j,i,m,k},
		(*The SW matrix is larger than sequence by 1 in each dimension to accomodate the
		  first row and column.*)
		matrix=ConstantArray[0, {Length[seq1]+1, Length[seq2]+1}];
		(*Put your code here.*)
		Do[M[i,1]:=0,{i,1,Length[seq1]+1}];
		Do[M[1,j]:=0,{j,1,Length[seq2]+1}];
		
		M[index1_,index2_]:=M[index1,index2]=Max[0,M[index1-1,index2-1]+scoringMatrix[[seq1[[index1-1]],seq2[[index2-1]]]],
		M[index1-1,index2]+scoringMatrix[[seq1[[index1-1]],5]],
		M[index1,index2-1]+scoringMatrix[[5,seq2[[index2-1]]]]];
		
		M[Length[seq1]+1,Length[seq2]+1];
		matrix=Table[M[m,k],{m,1,Length[seq1]+1},{k,1,Length[seq2]+1}];
		topScoringCell=Take[Position[matrix,Max[matrix]],-1][[1]];
		{topScoringCell, matrix}
	] 

(* swTraceback returns a list of pairs representing alignment columns. The elements of each pair are the 
   integers 1-5, with 5 representing a gap in the alignment.
   
   Implementation: Starting topScoringCell, trace back by checking whether the score in the current cell
   is what you would have gotten from each of the previous cells. The current cell should be updated to
   a previous cell that could yield the score of the current cell. 
   
   It is possible that more than one of the three previous adjacent cells would yield the score in the 
   current cell. If an alignment without gap and an alignment with a gap are available, always prefer the 
   one without a gap. If a gap in either sequence is possible, prefer to put the gap in seq2 over seq1.
   
   You may find AppendTo useful for adding alignment columns to the end of the list of columns. AppendTo
   is much faster than PrependTo but you'll have to use Reverse to reverse the order of the columns
   before returninging it. If you really want it to be fast you can look up and use Reap and Sow but this
   is not required.
 *)
swTraceback[{topScoringCell_, swMatrix_}, seq1_, seq2_, scoringMatrix_]:=
	Module[{alignmentColumns={}, 
			index1=topScoringCell[[1]], 
			index2=topScoringCell[[2]], 
			currentCellScore=Extract[swMatrix, topScoringCell],maxnext,l},
		While[currentCellScore > 0,
			(*Put your code here.*)
			maxnext=Max[swMatrix[[index1-1,index2]],swMatrix[[index1,index2-1]]];
			 If[maxnext>currentCellScore,
			 If[swMatrix[[index1-1,index2]]==maxnext,(alignmentColumns=AppendTo[alignmentColumns,{seq1[[index1-1]],5}])&&
			 (index1=index1-1),(alignmentColumns=AppendTo[alignmentColumns,{5,seq2[[index2-1]]}])&&
			 (index2=index2-1)]
			 ,(alignmentColumns=AppendTo[alignmentColumns,{seq1[[index1-1]],seq2[[index2-1]]}])&&(index1=index1-1)&&(index2=index2-1)];
			 currentCellScore=swMatrix[[index1,index2]]];
			 Reverse[alignmentColumns]
		]

(*prettyPrintAlignment prints an alignment of two DNA sequences in the tradition, BLAST-like format.*)
prettyPrintAlignment[displayMatrix_]:=
	Print[Grid[displayMatrix, Spacings -> {0, 0}]];
 
(* alignmentDisplayMatrix converts a list of alignment columns into the three three rows needed to 
	print it in the traditional, BLAST-like format. This matrix is passed into prettyPrintAlignment.
*)
alignmentDisplayMatrix[columns_]:=
	{
	columns[[All, 1]] /. {1->"A", 2->"C", 3->"G", 4->"T", 5->"_"},
	Map[If[#[[1]]==#[[2]], "|", " "] &,columns],
	columns[[All, 2]] /. {1->"A", 2->"C", 3->"G", 4->"T", 5->"_"}
	}


(* readScoreMatrix takes a pathname of a tab separated file that should contain a 5x5 matrix,
   with one line for each row of the matrix.*)
readScoringMatrix[filePath_] := 
	Module[{matrix}, 		
	matrix = Import[filePath];
	(*Probably should put a bunch of checks for the right dimensions and contents here.*)
	matrix
];

	
(* readFasta reads a fasta file and outputs the nucleotide sequence converted to numbers
INPUT
- fastaFile is a string representing the path to fasta file

OUTPUT
- input is a list of bases in the file indicated by fastaFile.  
  bases are translated form ACGT to 1234.
  e.g., {1,3,2,4,2}
*)
readFasta[fastaFile_]:=
	Flatten[Map[Characters, Import[fastaFile]] 
		   /. {"A"->1, "C"->2, "G"->3, "T"->4}
		   ]
