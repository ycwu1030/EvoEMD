
Begin["`EvoProcess`"]


(* Determine the External Particle Types/Numbers *)
GetExternalParticles[process_]:=Flatten[List@@(process)];

GetExternalParticlesByCate[process_,Cate_]:=Select[GetExternalParticles[process],FieldMatchQ[#,Cate]&];

GetExternalFermionNumber[process_]:=Length[GetExternalParticlesByCate[process,F] ];

GetExternalVectorNumber[process_]:=Length[GetExternalParticlesByCate[process,V] ];

GetExternalMasslessVectorNumber[process_]:=Block[{VL,Masses},
VL = GetExternalParticlesByCate[process,V];
Masses=TheMass/@VL;
Count[Masses,0]
]

GetExternalMassiveVectorNumber[process_]:= GetExternalVectorNumber[process]-GetExternalMasslessVectorNumber[process];

GetInitialParticles[process_]:=process[[1]];
GetFinalParticles[process_]:=process[[2]];

(* Add Graph information *)
GraphID /: Conjugate[GraphID[n_] ]:=GraphID[n];
AddGraphID[amp_]:=Block[{HD=Head[amp]},
HD@@(Cases[amp,FeynAmp[GraphID[args__,Number==n_],args1_,res_]:>FeynAmp[GraphID[args,Number==n],args1,GraphID[n] res] ])
]

(* Calculate Amplitudes *)
Mass2Particle:=(Rule[(TheMass[#]), #]) & /@ (Cases[F$Particles, _F | _V | _S]);
SquaredMass2Particle:=(Rule[(TheMass[#])^2, #]) & /@ (Cases[F$Particles, _F | _V | _S]);
SquaredMass2Mass[m2_]:=TheMass[m2//.SquaredMass2Particle];
SquaredMass2Mass[m_^2]:=m;
SquaredMass2Width[0]:=0
SquaredMass2Width[a_]:=Width[SquaredMass2Mass[a] ];

Width/:Conjugate[Width[args__] ]:=Width[args];
GetSquaredME[process_]:=Block[{top,ins,Namp,ampFA,ampWithID,amps,FermionMat,NF,FactorF,ColourMat,MME2,rep,MME2Tmp,ME2PS,ME2F},
ClearProcess[];
top=CreateTopologies[0,Length/@process];
ins=InsertFields[top,process,InsertionLevel->{Particles}];
ampFA=CreateFeynAmp[ins];
Namp=Length[ampFA];
ampWithID=TagDiagrams[ampFA,GraphID];
amps=CalcFeynAmp[ampWithID,FermionChains->Chiral,Invariants->True]/. {Den[S,b_]:>Den[S,b-I SquaredMass2Mass[b]SquaredMass2Width[b] ]};
(* Helicity Amplitude *)
(* If we don't have fermion chain in the amplitude, then it is not necessary to calculate HelicityME *)
FermionMat=If[FreeQ[ampFA,FermionChain[args___] ],{},_Hel=0; HelicityME[amps,amps] ];
NF=GetExternalFermionNumber[process];
FactorF=2^NF;
(* Colour Part *)
(* If we don't have and color index, then no need to calculate ColourME *)
ColourMat=If[FreeQ[ampFA,Index[Colour,___] ],{},ColourME[amps,amps] ];
(* Square the Amplitude *)
{MME2,rep}=SquaredME[amps,amps];
MME2Tmp=FactorF*MME2//.rep//.FermionMat//.ColourMat;
(* Polarization Summation *)
(* Only do the Polarization Summation when there is Polarization Vector in the amplitude *)
ME2PS=If[FreeQ[ampFA,PolarizationVector[args___] ], MME2Tmp, PolarizationSum[MME2Tmp] ];
ME2F=ME2PS//.Abbr[]//.Subexpr[];
COMS=DeleteDuplicates[Sort/@Tuples[Range[Namp],2] ];
ME2List[process]/@((Coefficient[ME2F,GraphID[#[[1]]]GraphID[#[[2]]] ]//Simplify)&/@COMS)
(* Table[Coefficient[ME2F,GraphID[i] GraphID[j] ]//Simplify, {i,1,Namp},{j,1,Namp} ] *)
]

(* Find the S/T/U expressions *)
Ei[sqrts_,m1_,m2_]:=(sqrts^2+m1^2-m2^2)/2/sqrts;
GetSTU[process_]:=Block[{ini,fin,Mini,Mfin,m1,m2,m3,m4,EqualI,EqualF,ResT,ResU,Res},
ini=process[[1]];
fin=process[[2]];
Mini=TheMass/@ini;
Mfin=TheMass/@fin;
EqualI=Mini[[1]]===Mini[[2]];
EqualF=Mfin[[1]]===Mfin[[2]];
m1=Mini[[1]];m2=Mini[[2]];m3=Mfin[[1]];m4=Mfin[[2]];
ResT={T->(m1^2+m3^2-2(Ei[sqrts,m1,m2] Ei[sqrts,m3,m4]-pMag[sqrts,Sequence@@(Sort[{m1,m2}])] pMag[sqrts,Sequence@@(Sort[{m3,m4}])] Cth))};
ResU=If[EqualF,
{U->(m1^2+m4^2-2(Ei[sqrts,m1,m2] Ei[sqrts,m4,m3]+pMag[sqrts,Sequence@@(Sort[{m1,m2}])] pMag[sqrts,Sequence@@(Sort[{m3,m4}])] Cth))},
{U->(m2^2+m3^2-2(Ei[sqrts,m2,m1] Ei[sqrts,m3,m4]+pMag[sqrts,Sequence@@(Sort[{m1,m2}])] pMag[sqrts,Sequence@@(Sort[{m3,m4}])] Cth))}];
Res=Join[{S->sqrts^2},ResT,ResU]//Simplify;
Abbreviate[Res,FreeQ[#,Cth]&&MatchQ[#,_pMag]&,MinLeafCount->1]
]

(* Obtain Numerator and Denominator *)
GetNumerator[amps_]:=amps//._Den->1;
GetDenominator[amps_]:=(Times@@Cases[TestSymbol amps,_Den|_Den^n_])//.TestSymbol->1;

(* Abbreviate the Numerator *)
SimplifyNumerator[expr_]:=Abbreviate[Collect[expr,{Cth}]//FullSimplify,FreeQ[#,Cth|sqrts]& ,MinLeafCount->5];
SimplifyDenominator[expr_]:=Abbreviate[Collect[expr,{Cth}]//FullSimplify,FreeQ[#,Cth|sqrts]& ,MinLeafCount->5];

(* Integral Over Solid Angle *)
IntSingleSquaredDiagram[amp_]:=Block[{process,RepSTU,amprep,num,den,numsp},
process=Head[amp][[1]];
RepSTU=GetSTU[process];
amprep=amp[[1]]//.{Den[S,args_]:>1/(S-args)}//FullSimplify;
num=GetNumerator[amprep]//.RepSTU;
den=GetDenominator[amprep]//.{Den[arg1_,arg2_]:>arg1-arg2}//.RepSTU;
numsp=SimplifyNumerator[num];
densp=SimplifyDenominator[den];
res=2 Pi Int[numsp/densp,{Cth,-1,1}];
Abbreviate[Collect[res,sqrts,Simplify],FreeQ[#,sqrts]&,MinLeafCount->5]
]

GetAllSubs[expr_]:=Block[{AllSubs, SubUsed},
AllSubs=ToExpression[Names["Sub*"] ];
SubinExpr=Intersection[Cases[expr,_Symbol,Infinity],AllSubs];
Subexpressions=Select[Subexpr[], MemberQ[SubinExpr,#[[1]]]&];
SubNew=Intersection[Cases[Subexpressions,_Symbol,Infinity],AllSubs];
NSubs=Length[SubinExpr];
NNext=Length[SubNew];
If[NSubs==NNext,OnePassOrder[OptimizeAbbr[Subexpressions] ],GetAllSubs[Subexpressions] ]
]

End[]
