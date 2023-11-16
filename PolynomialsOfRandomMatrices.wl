(* ::Package:: *)

(* :Title: PolynomialsOfRandomMatrices.wl *)

(* :Context: PolynomialsOfRandomMatrices` *)

(* :Authors: Torben Krueger and Yuriy Nemish *)

(* :Summary:
   Analysis of the Dyson equation and simumations for (self-adjoint) polynomials 
   of Wigner matrices and block-correlated random matrices
 *)

(* :Copyright: \[Copyright] 2023 Torben Krueger and Yuriy Nemish *)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 13.1 *)

(* :Sources:
   Linearization and minimization algorithms are discussed in 
   L.Erdos, T.Kruger, Yu.Nemish. Local laws for polynomials of Wigner matrices
   (arXiv:1804.11340)
*)

(* :Warnings:
   - NC variables are assumed to be self-adjoint
*)

(* :Limitations:
   
*)

(* :Discussion:
   
*)

(* :Requirements:
   NC`
   NCAlgebra`
*)

(* :Examples:
   Typical examples in PolynomialsOfRandomMatrices_Examples.nb
*)


(* NC`and NCAlgebra`has to be loaded before loading PolynomialsOfRandomMatrices`
	For example, do the following:
		- AppendTo[$Path,"path1/"]; where path1 is the path to the folder containing NC package
		- <<NC`
		- <<NCAlgebra`
		- <<"path/to/PolynomialsOfRandomMatrices.wl" *)



BeginPackage["PolynomialsOfRandomMatrices`"];


(* 1. Generating random matrices, random variables, polynomials *)
(* ************************************************************************************ *)
Ginibre::usage = "Ginibre[n] generates a (normalized, i.e. with variance 1/n) 
	Ginibre matrix of size n x n (with complex entries)";
GinibreReal::usage = "GinibreReal[n] generates a (normalized, i.e. with variance 1/n) 
	Ginibre matrix of size n x n (with real entries)";
IIDUnitDisk::usage = "IIDUnitDisk[n] generates a list if n complex numbers uniformly 
	distributed on the unit disk";
Wigner::usage = "Wigner[n] generates a (normalized, i.e. with variance 1/n) Wigner matrix 
	of size n x n (with complex entries)";
WignerReal::usage = "WignerReal[n] generates a (normalized, i.e. with variance 1/n) 
	Wigner matrix of size n x n with real entries";
GeneratePolynomial::usage = "Function GeneratePolynomial[VariableList, StructureVec, Type,
	PPrecision_:MachinePrecision] generates polynomial with randomly chosen coefficients of 
	type 'Type' with ~t_1 terms of degree 1, t_2 terms of degree 2 etc, (maybe a bit more 
	due to symmetrization) in variables from VariableList (in usual ordering i.e. x, y, xx, 
	xy, yx, yy, xxx...). 
	Type can be \"Integer\", \"ComplexInteger\", \"Real\", \"Complex\". Notice that z is 
	always a special variable. Notice that x, y, ... are self-adjont variables";
(* *********************************************************************************** *)
(* 2. MDE *)
(* *********************************************************************************** *)
SOperator2::usage = "SOperator2[Linearization,z,R] returns \[CapitalGamma][R], where \[CapitalGamma] is the 
	superoperator constructed from the Linearization matrix with spectral parameter z";
IterativelySolveMDE::usage = "IterativelySolveMDE[SOperator,ConstantMatrix,DirectionMatrix,z,M0,K,\[Epsilon]] 
	solves the MDE -M^{-1} = zJ-A+S[M] at a fixed spectral parameter z; 
	SOperator is the Superoperator S given as a function with matrices as input and output,
	ConstantMatrix is the constant matrix A, DirectionMatrix is the direction matrix J,
	z is the spectral parameter, M0 is the initial condition for the iteration, 
	K is the maximum number of steps in the iteration, \[Epsilon] is the accuracy"; 
IterativelySolveMDE2::usage = "IterativelySolveMDE[Linearization,z,Z,K,\[Epsilon]] solves the 
	MDE -M^{-1} = zJ-A+S[M] at a fixed spectral parameter z=Z; 
	exactly as IterativelySolveMDE, but SOperator, ConstantMatrix and DirectionMatrix 
	are determined from Linearization, and M0 = \[ImaginaryI]I; z is the spectral parameter (variable),
	Z is the value of the spectral parameter (number), K is the maximum number of steps in
	 the iteration, \[Epsilon] is the accuracy"; 
IterativelySolveMDEonInterval::usage = "IterativelySolveMDEonInterval[SOperator,ConstantMatrix,DirectionMatrix,E1,E2,\[Eta],L,K,\[Epsilon]] 
	solves the MDE along a parameter interval E+i\[Eta] with E \[Element] [E_1,E_2]; 
	Arguments are the same as for IterativelySolveMDE; 
	the extra arguments are E1 (left edge of the interval), E2 (right edge of the interval), 
	\[Eta] (Fixed value for Im[z]), L (Number of points in the interval on which the solution 
	is evaluated)";
IterativelySolveMDEonInterval2::usage = "IterativelySolveMDEonInterval2[Linearization,z,E1,E2,\[Eta],L,K,\[Epsilon]] 
	solves the MDE obtained from the Linearization matrix along a parameter interval E+i\[Eta] 
	with E \[Element] [E_1,E_2]; Arguments are the same as for IterativelySolveMDE2; 
	the extra arguments are E1 (left edge of the interval), E2 (right edge of the interval), 
	\[Eta] (Fixed value for Im[z]), L (Number of points in the interval on which the solution 
	is evaluated); Linearization matrix should already contain variable z";
(* *********************************************************************************** *)
(* 3. Linearization of polynomials *)
(* *********************************************************************************** *)
BigLinearization::usage = "BigLinearization[polynomial] returns the initial (big) 
	linearization of the polynomial in NC self-adjoint variables";
MinimalLinearization::usage = "Construct the minimal linearization in NC self-adjoint vars";
MinimalIndices::usage = "MinimalIndices[Pol,TType,PPrecision:MachinePrecision] returns 
	indices that give the minimal linearization; it's a modification of 
	MinimalLinearization that returns only the set of indices and not the linearization 
	itself";
ReconstructPolynomial::usage = "ReconstructPolynomial[linearization, \"Check\" (optional), 
	polynomial1 (optional), Precision (optional)]reconstructs polynomial2 from the given 
	linearization and (if we choose 'Check') compares it with given polynomial1, i.e., 
	compares the monomial coefficients and returns 'Reconstruction OK' if the max 
	difference is < 10^-5. The function returns ('Check'-> reconstructed polynomial2,
	conclusion, max difference of monomial coefficients between reconstructed polynomial2 
	and polynomial1; 'NoCheck' -> reconstructed polynomial). If the reconstruction is 
	not accurate enough, increase the precision)";
SuperOperatorToMatrix::usage = "SuperOperatorToMatrix[SuperOperator, dim] computes the 
	matrix for the given Superoperator; if SuperOperator is applied to matrix R of size n x n,
	then dim = n";


Begin["`Private`"];


Needs["NC`"]
Needs["NCAlgebra`"]


(* ::Chapter:: *)
(*1. Generating random matrices, random variables, polynomials*)


(* 1. Generating random matrices, random variables, polynomials *)

Ginibre[n_] :=
    (RandomVariate[NormalDistribution[], {n, n}] + I * RandomVariate[
        NormalDistribution[], {n, n}]) / Sqrt[2 * n];

GinibreReal[n_] :=
    RandomVariate[NormalDistribution[], {n, n}] / Sqrt[n];

IIDUnitDisk[n_] :=
    Block[{l, z},
        l = {};
        Until[
            Length[l] >= n
            ,
            z = RandomComplex[{-1 - I, 1 + I}];
            If[Abs[z] < 1,
                AppendTo[l, z]
            ];
        ];
        Return[l];
    ];

Wigner[n_] :=
    Block[{A},
        A = (RandomVariate[NormalDistribution[], {n, n}] + I * RandomVariate[
            NormalDistribution[], {n, n}]) / Sqrt[2 * n];
        Return[(A + ConjugateTranspose[A]) / Sqrt[2]];
    ];

WignerReal[n_] :=
    Block[{A},
        A = RandomVariate[NormalDistribution[], {n, n}] / Sqrt[n];
        Return[(A + Transpose[A]) / Sqrt[2]];
    ];

GeneratePolynomial[VariableList_,StructureVec_,CoefType_,PPrecision_:MachinePrecision]:=
	Block[{Poly,MainList,start,end,ProductList,CCoefficientList,LLength},
		MainList=Table[{VariableList[[i]]},{i,Length[VariableList]}];
		LLength=Sum[Length[VariableList]^k,{k,Length[StructureVec]}];
		start=1;end=Length[MainList]; 
		While[Length[MainList]<=LLength,
			For[i=start,i\[LessSlantEqual]end,i++,
				MainList=Join[MainList,Table[Append[MainList[[i]],VariableList[[j]]],
					{j,Length[VariableList]}]];
			];
		start=end+1;
		end=Length[MainList];
		];
		ProductList=Table[NonCommutativeMultiply@@MainList[[i]],{i,LLength}];

		Switch[CoefType,
		"Integer",
			CCoefficientList
				=Table[RandomChoice[Join[Range[-10,-1],Range[1,10]]],{i,LLength}];,
		"ComplexInteger",(*complex integer coefficients*)
			CCoefficientList
				=Table[RandomChoice[Join[Range[-10,-1],Range[1,10]]]
					+RandomInteger[{-10,10}]*I,{i,LLength}];,
		"Real",
			CCoefficientList=Table[RandomReal[{-10,10}],{j,LLength}];,
		"Complex",
			CCoefficientList
				=Table[(RandomReal[{-10,10}]+I*RandomReal[{-10,10}])/2,{j,LLength}];,
		"RealHighPrecision",
			CCoefficientList
				=Table[RandomReal[{-10,10},WorkingPrecision->PPrecision],{j,LLength}];,
		"ComplexHighPrecision",
			CCoefficientList
				=Table[(RandomReal[{-10,10},WorkingPrecision->PPrecision]
					+I*RandomReal[{-10,10},WorkingPrecision->PPrecision])/2,{j,LLength}];
		];

		start=0;end=0;
		For[i=1,i\[LessSlantEqual]Length[StructureVec],i++,
			start=end+1; end=end+Length[VariableList]^i;
			For[j=start,j<=end,j++,
				CCoefficientList[[j]]
					=CCoefficientList[[j]]*RandomVariate[
						BernoulliDistribution[StructureVec[[i]]/Length[VariableList]^i]];
			];
		];

		Poly=Sum[CCoefficientList[[i]]*ProductList[[i]],{i,LLength}];
		Poly=(Poly+ReplaceAll[aj[Poly], (* aj[x] is the adjoint of x (from NCAlgebra) *)
				Table[aj[VariableList[[i]]]-> VariableList[[i]],{i,Length[VariableList]}]])/2;
		Return[NCE[Simplify[Poly]]]
	];


(* ::Chapter:: *)
(*2. MDE*)


(* 2. MDE *)

SOperator2[Linearization_, z_(*spectral parameter*), R_] :=
    Block[{VariableList, LL},
        VariableList = DeleteCases[DeleteDuplicates[Cases[Linearization,
             _Symbol, Infinity]], z];
        LL = N[Table[D[Linearization, VariableList[[i]]], {i, Length[
            VariableList]}]];
        Return[Simplify[Total[Table[LL[[i]] . R . LL[[i]], {i, Length[LL]}]]]];
            
    ];
    
IterativelySolveMDE[
	SOperator_(*Superoperator S given as a function with matrices as input and output*), 
    ConstantMatrix_(*The constant matrix A*), 
    DirectionMatrix_(*The direction matrix J*), 
    z_(*The spectral parameter*), 
    M0_(*initial condition for the iteration*), 
    K_(*Maximum number of steps in the iteration*), 
    \[Epsilon]_(*Accuracy*)
    ] :=
    Block[{M1 = M0, M2, d, k = 1},
        M2 = -Inverse[z DirectionMatrix - ConstantMatrix + SOperator[
            M1]] / 2 + M1 / 2;
        d = Max[Abs[Flatten[M2 - M1]]];
        While[
            d > \[Epsilon] && k <= K
            ,
            k = k + 2;
            M1 = -Inverse[z DirectionMatrix - ConstantMatrix + SOperator[
                M2]] / 2 + M2 / 2;
            M2 = -Inverse[z DirectionMatrix - ConstantMatrix + SOperator[
                M1]] / 2 + M1 / 2;
            d = Max[Abs[Flatten[M2 - M1]]];
        ];
        Return[M2];
    ];

IterativelySolveMDE2[
	Linearization_(*Linearization matrix with zJ*), 
    z_(*spectral parameter variable*), 
    Z_(*The spectral parameter*), 
    M0_(*Initial condition for the iteration*), 
    K_(*Maximum number of steps in the iteration*), 
    \[Epsilon]_(*Accuracy*)
    ] :=
    Block[{M1, M2, d, k = 1, VariableList, ConstantMatrix, DirectionMatrix},
        VariableList = DeleteCases[DeleteDuplicates[Cases[Linearization,
             _Symbol, Infinity]], z];
        ConstantMatrix = N[Linearization /. Join[Table[VariableList[[
            i]] -> 0, {i, Length[VariableList]}], {z -> 0}]];
        DirectionMatrix = -N[D[Linearization, z]];
        M1 = M0;
        M2 = -Inverse[Z DirectionMatrix - ConstantMatrix + SOperator2[
            Linearization, z, M1]] / 2 + M1 / 2;
        d = Max[Abs[Flatten[M2 - M1]]];
        While[
            d > \[Epsilon] && k <= K
            ,
            k = k + 2;
            M1 = -Inverse[Z DirectionMatrix - ConstantMatrix + SOperator2[
                Linearization, z, M2]] / 2 + M2 / 2;
            M2 = -Inverse[Z DirectionMatrix - ConstantMatrix + SOperator2[
                Linearization, z, M1]] / 2 + M1 / 2;
            d = Max[Abs[Flatten[M2 - M1]]];
            If[Mod[k,100001]==0,Print[k,d]]
        ];
        Return[M2];
    ];

IterativelySolveMDEonInterval[
	SOperator_, 
	ConstantMatrix_, 
	DirectionMatrix_,
    E1_(*left edge of the interval*), 
    E2_(*right edge of the interval*),
    \[Eta]_(*Fixed value for Im[z]*), 
    L_(*Number of points in the interval on which the solution is evaluated*), 
    K_(*Maximum number of steps used in the computation of the values at each point E=Re[z]*), 
    \[Epsilon]_(*Desired
     precision*)
     ] :=
    Block[{EE, M0 = N[I * IdentityMatrix[Length[DirectionMatrix]]], M = {}},
        For[l = 0, l <= L, l++,
            EE = N[E1 + (E2 - E1) l / L];
            M0 = IterativelySolveMDE[SOperator, ConstantMatrix, DirectionMatrix,
                 EE + \[Eta] I, M0, K, N[\[Epsilon]]];
            (*The initial value for the iteration at the next E start
                 with the computed value of the solution for the previous E. 
                 This increases the speed if L is large*)
            M = Append[M, {EE, M0}]
        ];
        Return[M];
    ];

IterativelySolveMDEonInterval2[
	Linearization_, 
	z_, 
	E1_(*left edge of the interval*), 
	E2_(*right edge of the interval*), 
	\[Eta]_(*Fixed value for Im[z]*), 
	L_(*Number of points in the interval on which the solution is evaluated*), 
	K_(*Maximum number of steps used in the computation of the values at each 
		point E=Re[z]*), 
    \[Epsilon]_(*Desired precision*)
    ] :=
    Block[{EE, M0 = N[I * IdentityMatrix[Length[Linearization]]], M ={}},
        For[l = 0, l <= L, l++,
            EE = N[E1 + (E2 - E1) l / L];
            M0 = IterativelySolveMDE2[Linearization, z, EE + \[Eta] I, M0,
                 K, N[\[Epsilon]]];
            (*The initial value for the iteration at the next E start
                 with the computed value of the solution for the previous E. This increases
                 the speed if L is large*)
            M = Append[M, {EE, M0}]
        ];
        Return[M];
    ];



(* ::Chapter:: *)
(*3. Linearization of polynomials*)


(* 3. Linearization of polynomials *)

SymmetrizeCanonicalLinearization[Linearization_]:=
	Block[{dim=Length[Linearization](*Dimension of the linearization*),
	VariableList=DeleteDuplicates[Cases[Linearization,_Symbol,Infinity]]
	(*List of variables that appear in the linearization*)},
		Return[
			Simplify[
			ArrayFlatten[
				{
					{
					{{Linearization[[1,1]]}},
						Linearization[[1;;1,2;;dim]],
						ConjugateTranspose[Linearization[[2;;dim,1;;1]]]
					},
					{
					ConjugateTranspose[Linearization[[1;;1,2;;dim]]],
						SparseArray[{},{dim-1,dim-1}],
						2*ConjugateTranspose[Linearization[[2;;dim,2;;dim]]]
					},
					{
					Linearization[[2;;dim,1;;1]],
						2*Linearization[[2;;dim,2;;dim]],
						SparseArray[{},{dim-1,dim-1}]
					}
				}
				],
			Table[Element[VariableList[[i]],Reals],{i,Length[VariableList]}](*We assume self-adjointness of the variables*)
			]
		];
	];
	
BigLinearization[Polynomial_]:=
	Block[{Pol=NCExpand[Polynomial](*The polynomial is fully expanded*),dim},
		If[AtomQ[Pol],(*If the polynomial is only a variable or a number...*)
			Return[SparseArray[{{1,1}->Pol},{1,1}]//Normal]; 
			(*  then the linearization is simply the polynomial*)
			,
			If[Head[Pol]===NonCommutativeMultiply,(*If the polynomial is a monomial 
			of variables with coefficient 1, then we generate the standard linearization 
			for monomials.*)
				dim=Length[Pol];(*Number of factors*)
				Return[SparseArray[Join[Table[{i,dim+1-i}->Part[Pol,i],{i,dim}],
					Table[{i+1, dim+1-i}->-1,{i,dim-1}]],{dim,dim}]//Normal];
				,
				If[Head[Pol]===Times,(*If the polynomial is a monomial of variables 
				with arbitrary coefficient, then we reduce the linearization to the 
				linearization for monomials with coefficient 1. The coefficient is 
				distributed so that it appears symmetric if it is real. In fact if 
				{{0,c},{b,U}} is the linearization of p then 
				{{0,Sqrt[Abs[a]]*c},{Sqrt[Abs[a]]*b,(Abs[a]/a)*U}} is a linearization 
				of a*p [before a/Abs[a]!!!]*)
					Block[{PolLin},
						PolLin=BigLinearization[Part[Pol,2](*This is the monomial with 
						coefficient 1*)];
						If[Length[PolLin]==1 (*If the linearization is trivial, i.e. 
						the linearization of p is p, then we return a*p, where a is the 
						coefficient*),
							Return[Part[Pol,1](*This is the coefficient*)*PolLin//Normal];
							,
							PolLin[[1,Length[PolLin]]]
								=Sqrt[Abs[Part[Pol,1](*This is the coefficient*)]]
									*PolLin[[1,Length[PolLin]]];
							PolLin[[Length[PolLin],1]]
								=Sqrt[Abs[Part[Pol,1]]]*PolLin[[Length[PolLin],1]];
							(*PolLin[[2;;Length[PolLin],2;;Length[PolLin]]]
								=Part[Pol,1]/Abs[Part[Pol,1]]
									*PolLin[[2;;Length[PolLin],2;;Length[PolLin]]];*)
							PolLin[[2;;Length[PolLin],2;;Length[PolLin]]]
								=Abs[Part[Pol,1]]/Part[Pol,1]
									*PolLin[[2;;Length[PolLin],2;;Length[PolLin]]];
							Return[PolLin//Normal];
						];
					];
					,
					If[Head[Pol]===Plus,(*If the polynonial is a sum of monomials with 
						coefficients, then ...*)
						Block[{LinList,FirstElement,FirstRow,FirstColumn,Minor,
						BlockNumber,Linearization,
						VariableList=DeleteDuplicates[Cases[Pol,_Symbol,Infinity]]},
							LinList=Table[BigLinearization[Part[Pol,i]],{i,Length[Pol]}];
								(*The list of linearizations of the monomials with 
								coefficients*)
							FirstElement=Total[Flatten[Select[LinList,Length[#]==1&]]];
								(*The trivial monomials (variables or numbers) are summed 
								up and put into the (1,1)-entry*)
							LinList=Select[LinList,Length[#]>1&];
								(*We keep the linearizations of the non-trivial 
								monomials *)
							If[LinList=={},
								Return[{{FirstElement}}];
								,
								dim=Length/@LinList-1;
									(*List of the dimensions-1 of the linearizations 
									of the non-trivial monomials*)
								BlockNumber=Length[LinList];
									(*Number of linearizations of the non-trivial 
									monomials*)
								Minor=#[[2;;Length[#],2;;Length[#]]]&/@LinList;
									(*For the linearizations {{0,c},{b,U}} of the 
									monomials this returns the list of U's*)
								Minor=ArrayFlatten[
										Table[
											If[i==j,Minor[[i]],SparseArray[{},{dim[[i]],dim[[j]]}]
											],{i,BlockNumber},{j,BlockNumber}
										]
									];
									(*A block diagonal matrix with U's on the diagonal*)
								FirstColumn=Transpose[
									{Flatten[#[[2;;Length[#],1;;1]]&/@LinList,2]}
									];
									(*A column vector created by putting the columns b of the linearizations 
									{{0,c},{b,U}} of the monomials one after another*)
								FirstRow={Flatten[#[[1;;1,2;;Length[#]]]&/@LinList,2]};
									(*A row vector created by putting the rows c of the 
									linearizations {{0,c},{b,U}} of the monomials one 
									after another*)
								Linearization=ArrayFlatten[
									{
										{SparseArray[{{1,1}->FirstElement},{1,1}],FirstRow},
										{FirstColumn,Minor}
									}
									];
								(*The block matrix {{FirstElement,FirstRow},{FirstColumn,Minor}}*)
								If[
									Assuming[Table[Element[VariableList[[i]],Reals],{i,Length[VariableList]}],
										(*We assume that the variables are self-adjoint and check if the
										linearization is Hermitian*)
										HermitianMatrixQ[Linearization]
										],
									Return[Linearization];,(*If the linearization is already self-adjoint, 
										then we keep it as is*)
									Return[SymmetrizeCanonicalLinearization[Linearization](*If the 
										linearization is not self-adjoint, then we symmetrize it*)];
								];
							];
						];
						,
						Return["Not a proper polynomial"]
							(*If non of the above is applied, then the input is not a proper polynomial*)
					];
				];
			];
		];
	];

ExtendToBasis[VectorList_(*A family of vectors V*)]:=
	Block[{BasisList=VectorList(*We start generating the family spanning the entire 
		space with the given family*),
		dim=Length[VectorList[[1]]](*dimension of the entire space*)
		},
		Table[
			If[MatrixRank[Append[BasisList,SparseArray[{i}->1,dim]//Normal]]
				>Length[BasisList],
				(*If the space spanned by BasisList does not contain e_i, then ...*)
				BasisList=Append[BasisList,SparseArray[{i}->1,dim]//Normal]
				(*We add e_i to the BasisList (At the end a minimal family E has been added)*)
			]
		,{i,dim}
		];(*We go through all canonical basis vectors e_i one by one*)
		Return[BasisList];(*We return the spanning family of V+E*)
	];

GenerateSpace[Vector_,MatrixList_,TType_,PPrecision_:MachinePrecision]:=
	(* given Vector v and MatrixList M_i this program returns the basis of the space 
	spanned by {v, M_iv, M_jM_iv ,...} and the corresponding indices *)
	Block[{AbstractBasis,VectorBasis,NewVector,AbstractAll,VectorAll,
		start,end,startAll,endAll,matrixRankTolerance,$MinPrecision},
			
		Switch[TType,
			"Exact", matrixRankTolerance=0;,
			"MachinePrecision", $MinPrecision=MachinePrecision; matrixRankTolerance=(0.1)^8;,
			"HighPrecision", $MinPrecision=PPrecision; matrixRankTolerance=(0.1)^8;
			];
		
		AbstractAll={{}};(* this will be a complete set of indices in which 
			we will be looking for the basis indices *)
		VectorAll={Vector}; (* this will be a complete set of vectors *)
		AbstractBasis={{}};(*  AbstractBasis=\tilde{I}_U contains the emptyset... *)
		VectorBasis={Vector}; (* which corresponds to the vector in the VectroBasis *)
		start=1;(* start gives the position of the first basis vector of "length" 
			L ..; Initially VectorBasis1 has only one vector of "length" 0 *)
		end=1;(* end gives the position of the last basis vector of "length" 
			L ..; Initially VectorBasis1 has only one vector of "length" 0 *)
		startAll=1;(* startAll gives the position of the first index of "length" 
			L ..; in AbstractAll *) 
		endAll=1;(* endAll gives the position of the last index of "length" L .. 
			in AbstractAll  *)

			
		While[start<=end, (* We stop the iteration if no new "longer" basis vectors 
					were added after doing the following*)
					For[i=startAll,i<=endAll,i++,(* we take all the vectors of the highest 
						order computed so far (that were in fact added on the previous step) ..*)
						For[k=1,k<=Length[MatrixList],k++, (* and we compute vectors of 
							higher degree ..*)
							NewVector=MatrixList[[k]] . VectorAll[[i]]; 		
							If[MatrixRank[Append[VectorBasis,NewVector],Tolerance->matrixRankTolerance]
								>Length[VectorBasis], (* and if this vector is linearly 
									independent of the other basis vectors that are 
									already in VectorBasis.. *)
								VectorBasis=Append[VectorBasis,NewVector]; 
									(* we add the new vector to VectorBasis ..*)
								AbstractBasis=Append[AbstractBasis,Prepend[AbstractAll[[i]],k]]; 
									(* and add the corresponding index to AbstractBasis *)
							];
							AbstractAll=Append[AbstractAll,Prepend[AbstractAll[[i]],k]]; 
								(* and we add all the indices to AbstractAll *)
							VectorAll=Append[VectorAll,NewVector]; (* we add all the 
								vectors we computed to the VectroAll *)
						];
					];
					start=end+1;(* The basis vectors of the biggest "length" should 
						start here.. *)
					end=Length[AbstractBasis]; (* and end here *)
					startAll=endAll+1;(* the newly added indices of the highest degree 
						start here.. *)
					endAll=Length[AbstractAll]; (* and end here *)
				];
				
		VectorBasis=Transpose[VectorBasis]; (* The basis vectors will be 
					considered as columns*)
				
	Return[{AbstractBasis,VectorBasis}];
	
	];

MinimalLinearization[Pol_,TType_,PPrecision_:MachinePrecision]:=
	Block[{InLin,VariableList,ConstantMatrix,ConstantMatrixInverse,LinearizationMatrixList,
			E1,AbstractBasis1,VectorBasis1,PU,PUInverse,PUL,DualVector,AbstractBasis2,
			VectorBasis2,VectorBasis3,A,NewConstantMatrix,NewLinearizationMatrixList,
			NewLinearization,W1,W,NewVector,$MinPrecision},
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Initial linearization *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	
	Switch[TType, (* define InLin : Initial (big) linearization \[Rule] L and set the precision*) 
			"Exact", InLin=BigLinearization[Pol];,
			"MachinePrecision", $MinPrecision=MachinePrecision; InLin=N[BigLinearization[Pol]];,
			"HighPrecision", $MinPrecision=PPrecision; InLin=N[BigLinearization[Pol],PPrecision];
			];	 
	VariableList=DeleteDuplicates[Cases[InLin,_Symbol,Infinity]]; (* List of variables *)
	ConstantMatrix=InLin/.Table[VariableList[[i]]->0,{i,Length[VariableList]}]; 
			(* Constant matrix of the initial (big) linearization \[Rule] K_0 ... *)
	Switch[TType, (* ... and its inverse *)
			"Exact", ConstantMatrixInverse=Simplify[Inverse[ConstantMatrix]];,
			"MachinePrecision", ConstantMatrixInverse
				=LinearSolve[ConstantMatrix,IdentityMatrix[Length[ConstantMatrix]]];,
			"HighPrecision", ConstantMatrixInverse
				=LinearSolve[ConstantMatrix,IdentityMatrix[Length[ConstantMatrix]]];
			];
			
	Switch[TType, (* ... and the (non symmetric) linearization matrices K_i K_0^(-1) *)
			"Exact", LinearizationMatrixList
				=Table[Simplify[D[InLin,VariableList[[i]]] . ConstantMatrixInverse],{i,Length[VariableList]}];,
			"MachinePrecision", LinearizationMatrixList
				=Table[D[InLin,VariableList[[i]]] . ConstantMatrixInverse,{i,Length[VariableList]}];,
			"HighPrecision", LinearizationMatrixList
				=Table[D[InLin,VariableList[[i]]] . ConstantMatrixInverse,{i,Length[VariableList]}];
			];
	Switch[TType, (* ... and matrix E1 *)
			"Exact", E1=SparseArray[{{1}->1},Length[InLin]]//Normal;,
			"MachinePrecision", E1=SparseArray[{{1}->1},Length[InLin]]//Normal;,
			"HighPrecision", E1=N[SparseArray[{{1}->1},Length[InLin]]//Normal,PPrecision];
			];
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* First stage of the minimization procedure: constructing the subspace U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	{AbstractBasis1,VectorBasis1}=GenerateSpace[E1,LinearizationMatrixList,TType];
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* PU is the projection on U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	Switch[TType, 
			"Exact", PUInverse
				=Simplify[Inverse[ConjugateTranspose[VectorBasis1] . VectorBasis1,Method->"DivisionFreeRowReduction"]];,
			"MachinePrecision", PUInverse=
				LinearSolve[
					ConjugateTranspose[VectorBasis1] . VectorBasis1,N[IdentityMatrix[Length[Transpose[VectorBasis1]]]]
					];,
			"HighPrecision", PUInverse
				=LinearSolve[
					ConjugateTranspose[VectorBasis1] . VectorBasis1,IdentityMatrix[Length[Transpose[VectorBasis1]]]
					];
			];
	Switch[TType, (* ... and we compute PU=VectorBasis.(VectorBasis^*.BectorBasis)^{-1}BectorBasis^*,  *)
			"Exact", PU=Simplify[VectorBasis1 . PUInverse . ConjugateTranspose[VectorBasis1]];,
			"MachinePrecision", PU=VectorBasis1 . PUInverse . ConjugateTranspose[VectorBasis1];,
			"HighPrecision", PU=VectorBasis1 . PUInverse . ConjugateTranspose[VectorBasis1];
			];
	Switch[TType, 
			"Exact", PUL
				=Table[Simplify[PU . ConjugateTranspose[LinearizationMatrixList[[i]]]],{i,Length[VariableList]}];,
			"MachinePrecision", PUL=Table[PU . ConjugateTranspose[LinearizationMatrixList[[i]]],{i,Length[VariableList]}];,
			"HighPrecision", PUL=Table[PU . ConjugateTranspose[LinearizationMatrixList[[i]]],{i,Length[VariableList]}]; 
			];
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the indices \mathcal{I}_\tilde{U} that correspond to \tilde{U} *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	Switch[TType, 
			"Exact", DualVector=Simplify[PU . ConstantMatrixInverse . E1];,
			"MachinePrecision", DualVector=PU . ConstantMatrixInverse . E1;,
			"HighPrecision", DualVector=PU . ConstantMatrixInverse . E1; 
			];
	
	{AbstractBasis2,VectorBasis2}=GenerateSpace[DualVector,PUL,TType];
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the basis that allows rewriting the minimal linearization using the K_i and K_0^{-1} *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	VectorBasis3={E1};
	
	For[i=2,i<=Length[AbstractBasis2],i++,
		NewVector=LinearizationMatrixList[[AbstractBasis2[[i]][[1]]]];
		For[j=2,j<=Length[AbstractBasis2[[i]]],j++,
			NewVector=NewVector . LinearizationMatrixList[[AbstractBasis2[[i]][[j]]]];
		];
		NewVector=NewVector . E1;
		VectorBasis3=Append[VectorBasis3,NewVector];
	];
	VectorBasis3=Transpose[VectorBasis3];
	
	A=VectorBasis3; (* this is the matrix A from the notes *)

(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the minimal (not yet canonical) linearization *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	Switch[TType, (* \tilde{K}_0 = A^*.K_0^{-1}.A *)
			"Exact", NewConstantMatrix
				=Simplify[ConjugateTranspose[VectorBasis3] . ConstantMatrixInverse . VectorBasis3];,
			"MachinePrecision", NewConstantMatrix
				=ConjugateTranspose[VectorBasis3] . ConstantMatrixInverse . VectorBasis3;,
			"HighPrecision", NewConstantMatrix
				=ConjugateTranspose[VectorBasis3] . ConstantMatrixInverse . VectorBasis3; 
			];
	Switch[TType,
			"Exact", NewLinearizationMatrixList
				=Table[
				Simplify[
				ConjugateTranspose[VectorBasis3] . ConstantMatrixInverse . LinearizationMatrixList[[i]] . VectorBasis3
				],{i,Length[VariableList]}
				];,
			"MachinePrecision", NewLinearizationMatrixList
				=Table[
				ConjugateTranspose[VectorBasis3] . ConstantMatrixInverse . LinearizationMatrixList[[i]] . VectorBasis3
					,{i,Length[VariableList]}
				];,
			"HighPrecision", NewLinearizationMatrixList
				=Table[
				ConjugateTranspose[VectorBasis3] . ConstantMatrixInverse . LinearizationMatrixList[[i]] . VectorBasis3
					,{i,Length[VariableList]}
				];
			];
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the cannonical linearization *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	Switch[TType, (* W1 gives the norm of the first row/column of the constant matrix of the minimal but 
					not canonical linearization (\tilde{K_0}) *)
			"Exact", W1=Simplify[Norm[NewConstantMatrix[[1]]]];,
			"MachinePrecision", W1=Norm[NewConstantMatrix[[1]]];,
			"HighPrecision", W1=Norm[NewConstantMatrix[[1]]];
			];
	Switch[TType, (* we construct a unitary matrix W such that its first row is parallel 
				to the column of \tilde{K_0} (conjugate of the first row) *)
			"Exact", W=Simplify[Orthogonalize[ExtendToBasis[{Conjugate[NewConstantMatrix[[1]]]/W1}]]];,
			"MachinePrecision", W=Orthogonalize[ExtendToBasis[{Conjugate[NewConstantMatrix[[1]]]/W1}]];,
			"HighPrecision", W=Orthogonalize[ExtendToBasis[{Conjugate[NewConstantMatrix[[1]]]/W1}]];
			];

	W=Transpose[W]; (* now the first column is parallel to the first column *)

	Switch[TType, (* now we get a constant matrix of the canonical minimal linearization .. *)
			"Exact", NewConstantMatrix=Simplify[ConjugateTranspose[W] . NewConstantMatrix . W/(W1*W1)];,
			"MachinePrecision", NewConstantMatrix=ConjugateTranspose[W] . NewConstantMatrix . W/(W1*W1);,
			"HighPrecision", NewConstantMatrix=ConjugateTranspose[W] . NewConstantMatrix . W/(W1*W1);
			];
	Switch[TType, (* and all the rest linearization matrices of the canonical minimal linearization *)
			"Exact", NewLinearizationMatrixList=
				Table[Simplify[ConjugateTranspose[W] . NewLinearizationMatrixList[[i]] . W/(W1*W1)],
						{i,Length[VariableList]}];,
			"MachinePrecision", NewLinearizationMatrixList=
				Table[ConjugateTranspose[W] . NewLinearizationMatrixList[[i]] . W/(W1*W1),
					{i,Length[VariableList]}];,
			"HighPrecision", NewLinearizationMatrixList=
				Table[ConjugateTranspose[W] . NewLinearizationMatrixList[[i]] . W/(W1*W1),
					{i,Length[VariableList]}];
			];
	Switch[TType, (* and compute the canonical minimal linearization*)
			"Exact", NewLinearization=
				Simplify[
			NewConstantMatrix+Total[Table[VariableList[[k]]*NewLinearizationMatrixList[[k]],{k,Length[VariableList]}]]
				];,
			"MachinePrecision", NewLinearization=
			NewConstantMatrix+Total[Table[VariableList[[k]]*NewLinearizationMatrixList[[k]],{k,Length[VariableList]}]];,
			"HighPrecision", NewLinearization=
			NewConstantMatrix+Total[Table[VariableList[[k]]*NewLinearizationMatrixList[[k]],{k,Length[VariableList]}]];
			];
	Return[
		NewLinearization
		];
	];



MinimalIndices[Pol_,TType_,PPrecision_:MachinePrecision]:=
	Block[{InLin,VariableList,ConstantMatrix,ConstantMatrixInverse,LinearizationMatrixList,
			E1,AbstractBasis1,VectorBasis1,PU,PUInverse,PUL,DualVector,AbstractBasis2,
			VectorBasis2,VectorBasis3,A,NewConstantMatrix,NewLinearizationMatrixList,
			NewLinearization,W1,W,NewVector},
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Initial linearization *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
	
	Switch[TType,
	"Exact",
		Block[{},
			InLin=BigLinearization[Pol]; (* InLin : Initial (big) linearization \[Rule] L *) 
			VariableList=DeleteDuplicates[Cases[InLin,_Symbol,Infinity]]; (* List of variables *)
			
			ConstantMatrix=InLin/.Table[VariableList[[i]]->0,{i,Length[VariableList]}]; 
					(* Constant matrix of the initial (big) linearization \[Rule] K_0 ... *)
			ConstantMatrixInverse=Simplify[Inverse[ConstantMatrix]]; (* ... and its inverse *)
			LinearizationMatrixList
				=Table[Simplify[D[InLin,VariableList[[i]]] . ConstantMatrixInverse],{i,Length[VariableList]}]; 
					(* and the (non symmetric) linearization matrices K_i K_0^(-1) *)
			E1=SparseArray[{{1}->1},Length[InLin]]//Normal;
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* First stage of the minimization procedure, constructing of the subspece U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			{AbstractBasis1,VectorBasis1}=GenerateSpace[E1,LinearizationMatrixList,TType];
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* PU is the projection on U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			
			PU=Simplify[
				VectorBasis1 . Inverse[
						ConjugateTranspose[VectorBasis1] . VectorBasis1, Method->"DivisionFreeRowReduction"
								] . ConjugateTranspose[VectorBasis1]
					];(* and PU=VectorBasis.(VectorBasis^*.BectorBasis)^{-1}BectorBasis^*,  *)
			PUL=Table[Simplify[PU . ConjugateTranspose[LinearizationMatrixList[[i]]]],{i,Length[VariableList]}];
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the indices \mathcal{I}_\tilde{U} that correspond to \tilde{U} *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			DualVector=Simplify[PU . ConstantMatrixInverse . SparseArray[{{1}->1},Length[InLin]]//Normal];
			{AbstractBasis2,VectorBasis2}=GenerateSpace[DualVector,PUL,TType];
		
		];
	,
	(*************)
	"MachinePrecision",
	(*************)
		Block[{},
			InLin=N[BigLinearization[Pol]]; (* InLin : Initial (big) linearization \[Rule] L *) 
			VariableList=DeleteDuplicates[Cases[InLin,_Symbol,Infinity]]; 
					(* List of variables *)
			ConstantMatrix=InLin/.Table[VariableList[[i]]->0,{i,Length[VariableList]}]; 
					(* Constant matrix of the initial (big) linearization \[Rule] K_0 ...*) 
			ConstantMatrixInverse
				=LinearSolve[ConstantMatrix,N[IdentityMatrix[Length[ConstantMatrix]]]]; 
						(* ... and its inverse *)
			LinearizationMatrixList
				=Table[D[InLin,VariableList[[i]]] . ConstantMatrixInverse,{i,Length[VariableList]}]; 
					(* Matrices that give the initial non-symmetric linearization \[Rule] K_i .K_0^{-1}*)
			
			E1=SparseArray[{{1}->1},Length[InLin]]//Normal; 
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* First stage of the minimization procedure: constructing the subspece U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			{AbstractBasis1,VectorBasis1}=GenerateSpace[E1,LinearizationMatrixList,TType];
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* PU is the projection on U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			PUInverse
				=LinearSolve[
					ConjugateTranspose[VectorBasis1] . VectorBasis1,
					N[IdentityMatrix[Length[Transpose[VectorBasis1]]]]
					];
			PU=VectorBasis1 . PUInverse . ConjugateTranspose[VectorBasis1];
					(* and PU=VectorBasis.(VectorBasis^*.BectorBasis)^{-1}BectorBasis^*,  *) 
			PUL=Table[PU . ConjugateTranspose[LinearizationMatrixList[[i]]],{i,Length[VariableList]}]; 
	
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the indices \mathcal{I}_\tilde{U} that correspond to \tilde{U} *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			  DualVector=PU . ConstantMatrixInverse . E1;
			
			{AbstractBasis2,VectorBasis2}=GenerateSpace[DualVector,PUL,TType];
			
		];
	,
	(*************)
	"HighPrecision",
	(*************)
		Block[{$MinPrecision=PPrecision},
			InLin=N[BigLinearization[Pol],PPrecision]; (* InLin : Initial (big) linearization \[Rule] L *) 
			VariableList=DeleteDuplicates[Cases[InLin,_Symbol,Infinity]]; (* List of variables *)
			ConstantMatrix=InLin/.Table[VariableList[[i]]->0,{i,Length[VariableList]}]; 
					(* Constant matrix of the initial (big) linearization \[Rule] K_0 ...*) 
			ConstantMatrixInverse=LinearSolve[ConstantMatrix,IdentityMatrix[Length[ConstantMatrix]]]; 
					(* ... and its inverse *)
			LinearizationMatrixList
				=Table[D[InLin,VariableList[[i]]] . ConstantMatrixInverse,{i,Length[VariableList]}]; 
					(* Matrices that give the initial non-symmetric linearization \[Rule] K_i .K_0^{-1}*)
			
			E1=N[SparseArray[{{1}->1},Length[InLin]]//Normal,PPrecision]; 
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* First stage of the minimization procedure: constructing the subspece U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			{AbstractBasis1,VectorBasis1}=GenerateSpace[E1,LinearizationMatrixList,TType];

(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* PU is the projection on U *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			PUInverse
				=LinearSolve[
					ConjugateTranspose[VectorBasis1] . VectorBasis1,
					IdentityMatrix[Length[Transpose[VectorBasis1]]]
					];
			PU=VectorBasis1 . PUInverse . ConjugateTranspose[VectorBasis1];
					(* and PU=VectorBasis.(VectorBasis^*.BectorBasis)^{-1}BectorBasis^*,  *) 
			PUL=Table[PU . ConjugateTranspose[LinearizationMatrixList[[i]]],{i,Length[VariableList]}]; 
			
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
(* Constructing the indices \mathcal{I}_\tilde{U} that correspond to \tilde{U} *)
(*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*)
			  DualVector=PU . ConstantMatrixInverse . E1;
			
			{AbstractBasis2,VectorBasis2}=GenerateSpace[DualVector,PUL,TType];
		
		];
	
	];
	
	Return[AbstractBasis2];
	];

NonCommDegree[Poly_]:=
	Block[{VariableList,PolyDegree,DDegree},
		VariableList=DeleteDuplicates[Cases[Poly,_Symbol,Infinity]]; (* List of variables *)
		PolyDegree=Total[NCMonomialList[Poly,VariableList]]; (* create a new polynomial 
				by replacing all the coefficients by 1 ..*)
		PolyDegree=PolyDegree/.Table[VariableList[[i]]->x,{i,Length[VariableList]}]; 
				(* and replacing all the variables by x; of course degree of PolyDegree 
				is the same as degree of Poly (no cancelations since we changed all the 
				coefficients to 1 *)
		PolyDegree=D[PolyDegree,x]; (* now get the degree by differentiating wrt x*)
		DDegree=0;
		While[(PolyDegree/.x->1)>10^-10,
			DDegree=DDegree+1;
			PolyDegree=D[PolyDegree,x];
			];
		Return[DDegree];
	];

ReconstructPolynomial[CCanLin_,CCheck_:"NoCheck",Poly_:0,PPrecision_:MachinePrecision]:=
	Block[{VariableList,ConstantMatrix,LinearizationMatrixList,ConstantMinor,
	MinorsList,ConstantMinorInverse,MCList,start,end,MinorInverse,ReconstructedPolynomial,
	NewMatrix,DDegree,MonomList,CCheckResult,Nilpotency,MCListExt,MonomListExt,
	LinearizationMatrixListExt,VariableListExt},
	
		Block[{$MinPrecision=PPrecision},
			If[Length[CCanLin]==1,
				ReconstructedPolynomial=CCanLin[[1,1]], (* if the polynomial is linear, 
					then it should be in the (1,1)-component of the 1-dimensional 
					linearization, if not then we have to reconstruct it...*)
				VariableList=DeleteDuplicates[Cases[CCanLin,_Symbol,Infinity]]; (* List 
					of variables *)
				ConstantMatrix
					=SetPrecision[
						CCanLin/.Table[VariableList[[i]]->0,{i,Length[VariableList]}],
						PPrecision
						]; (* now we get a constant matrix of the canonical linearization ..*)
				LinearizationMatrixList
					=Table[
					SetPrecision[D[CCanLin,VariableList[[i]]],PPrecision],{i,Length[VariableList]}
						];(* and all the rest linearization matrices of the canonical linearization  *)
				
				(* Minors *)
				ConstantMinor
					=ConstantMatrix[[2;;Length[ConstantMatrix],2;;Length[ConstantMatrix]]]; (*K0 (minor)*)
				(*ConstantMinorInverse=Inverse[ConstantMinor];(* K0inv (minor) *)*)
				ConstantMinorInverse
					=LinearSolve[
						ConstantMinor,
						SetPrecision[IdentityMatrix[Length[ConstantMinor]],
						PPrecision]
						];(* K0inv (minor) *) 
				
				MinorsList
					=Table[
						SetPrecision[
							LinearizationMatrixList[[i]][[2;;Length[ConstantMatrix],2;;Length[ConstantMatrix]]]
								. ConstantMinorInverse,
								PPrecision
							],
						{i,Length[VariableList]}
						]; (*K_iK0inv minor*)
				
				(* To compute the expansion of the inverse of the minor, we need all the matrix 
					coefficients, all non-zero products of minor(L_i)*minor(L_j)...*)
				MCList={};(*here we will collect all non-zero products ... *)
				MonomList={}; (*and corresponding monomials *)
				For[k=1,k<=Length[VariableList],k++,
					If[VariableList[[k]]=!=z ,
						MCList=Append[MCList,-MinorsList[[k]]]; (* first we have all the 
							"products of length 1", note the minus sign due to 
							(1+x)^-1 = 1-x+x^2-x^3... *)
						MonomList=Append[MonomList,VariableList[[k]]];
					];
				];
				
				(*finding the degree of the polynomial *)
				If[Poly===0,
					start=1;
					end=Length[MCList];
					Nilpotency=Total[Table[Norm[MCList[[l]]],{l,start,end}]];
					While[N[Nilpotency]>=10^(-5), (*in fact DDegree-1 would be enough*)
						For[i=start,i<=end,i++,
							For[j=1,j<=Length[VariableList],j++,
								NewMatrix=-MCList[[i]] . MinorsList[[j]]; (*note the minus 
									sign due to geom expantion (1+x)^-1 = 1 + sum (-1)^k x^k *) 
								If[VariableList[[j]]=!=z 
										&& NewMatrix=!= ConstantArray[0,
												{Length[ConstantMinor],Length[ConstantMinor]}],
									MCList=Append[MCList,NewMatrix]; 
									MonomList=Append[MonomList,MonomList[[i]]**VariableList[[j]]];
								];
							];
						];
						start=end+1;
						end=Length[MCList];
						Nilpotency=Total[Table[Norm[MCList[[l]]],{l,start,end}]];
					];,
					DDegree=NonCommDegree[Poly]
				]; (* if we know the original polynomial, then we can make things more 
					efficient with DDegree=degree(polynomial) *)
				start=1;
				end=Length[MCList];
				For[k=1,k<=DDegree,k++, (*in fact DDegree-1 would be enough*)
					For[i=start,i<=end,i++,
						For[j=1,j<=Length[VariableList],j++,
							NewMatrix=-MCList[[i]] . MinorsList[[j]]; (*note the minus sign 
									due to geom expantion (1+x)^-1 = 1 + sum (-1)^k x^k *) 
							If[VariableList[[j]]=!=z 
								&& NewMatrix=!= ConstantArray[0,
										{Length[ConstantMinor],Length[ConstantMinor]}],
								MCList=Append[MCList,NewMatrix]; 
								MonomList=Append[MonomList,MonomList[[i]]**VariableList[[j]]];
							];
						];
					];
					start=end+1;
					end=Length[MCList];
				];
				
				MCList=Table[ConstantMinorInverse . MCList[[i]],{i,Length[MCList]}]; 
							(* add K0inv on the left*) 
				
				(*MinorInverse=ConstantMinorInverse + Total[Table[MCList[[i]]*MonomList[[i]],{i,Length[MCList]}]];*)
				
				(*Note that the first row/column of the big linearization do not have constant term (except the 1,1-entry)*)
				
				MCListExt = Append[MCList, ConstantMinorInverse]; (*now we define the 
					extended minors list, variable list etc by adding the corresponding 
					constant matrices and 1 (to the monomial list); we use them to reconstruct the polynomial in one
					loop*)
				MonomListExt = Append[MonomList, 1];
				LinearizationMatrixListExt = Append[LinearizationMatrixList, ConstantMatrix];
				VariableListExt = Append[VariableList, 1];
				
				ReconstructedPolynomial
					=NCExpand[
					CCanLin[[1,1]]
					-Total[
							Table[
								LinearizationMatrixListExt[[i]][[1,2;;Length[ConstantMatrix]]]
									. MCListExt[[k]]
									. LinearizationMatrixListExt[[j]][[2;;Length[ConstantMatrix],1]]
									*VariableListExt[[i]]**MonomListExt[[k]]**VariableListExt[[j]]
								,{i,Length[VariableListExt]},{j,Length[VariableListExt]},{k,Length[MCListExt]}
								],Infinity
							]
					];
			];
			
			If[Poly===0,
				VariableList=DeleteDuplicates[Cases[CCanLin,_Symbol,Infinity]];,
				VariableList=DeleteDuplicates[Cases[Poly,_Symbol,Infinity]];
			]; 
				(* to make sure that we have all the variables in case something 
					is wrong with the linearization *)
			CCheckResult
				=Max[Abs/@NCCoefficientList[Poly-ReconstructedPolynomial,VariableList]];
		];
	
	
	Switch[CCheck,
	"Check",
		If[SetPrecision[
				Poly-ReconstructedPolynomial/.Table[VariableList[[i]]->1,{i,Length[VariableList]}],
				PPrecision
			]==0,
			Return[{N[Chop[ReconstructedPolynomial,10^-5],5],
				Style["Reconstruction: OK",FontColor->Green],0}],
			If[CCheckResult<=10^-5,
				Return[{N[Chop[ReconstructedPolynomial,10^-5],5],
					Style["Reconstruction: OK",FontColor->Green],CCheckResult}],
				Return[{N[Chop[ReconstructedPolynomial,10^-5],5],
					Style["Reconstruction: Error",FontColor->Red],CCheckResult}]
			]
		],
	"NoCheck",
		Return[{Chop[N[NCExpand[ReconstructedPolynomial],5],10^-5]
			(*here we can use Chop[] in case we have some noise from computations*)}];
	];
	
	];

    
SuperOperatorToMatrix[SuperOperator_, dim_]:=
	Block[{Matrix,NewLine,E,SE},
		Matrix={};
	For[ss=1, ss<=dim, ss++,
		For[tt=1, tt<=dim, tt++,
			NewLine={};
			E=SparseArray[{{ss,tt}->1},{dim,dim}];
			SE=SuperOperator[E];
			For[k=1,k<=dim,k++,
				NewLine=Join[NewLine,SE[[k]]]
			];
			Matrix=Append[Matrix,NewLine];
		];
	];
	Return[Transpose[Matrix]];
];



End[];
EndPackage[];
