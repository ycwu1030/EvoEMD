

Begin["`EvoInit`"]

(* Checking whether the particle is in the model and is valid *)
(* These functions come from FeynArts, but they are not exposed, so we copied here *)
ValidIndexQ[{i_, j___}, {fi_, f___}] :=
  If[ IntegerQ[i] && !MemberQ[IndexRange[fi], i], False, ValidIndexQ[{j}, {f}] ]

ValidIndexQ[{}, _] = True

_ValidIndexQ = False


InModelQ[s_. fi_[i_Integer, j_List]] := InModelQ[s fi[i]] &&
  Block[{NoUnfold = Identity}, ValidIndexQ[j, Indices[fi[i]]]]

InModelQ[p_] := MemberQ[F$AllowedFields, p]
(* End of Copy *)

Attributes[ModelHasParticle]={Listable};
ModelHasParticle[part_]:=If[InModelQ[part],True,Print[ToString[part]<>" does not exist in model "<>StringRiffle[$Model, {"(","+",")"}] ];False ];


RelevantProcess[init_->final_, pois_]:=Block[{PartInProcess,GoodProc,ReleProc},
    PartInProcess = Join[Flatten[{init}],Flatten[{final}] ];
    GoodProc = And[ModelHasParticle[PartInProcess] ];
    ReleProc = IntersectingQ[PartInProcess,pois];
    GoodProc && ReleProc
]

(* From a process X to {X, CP[X], T[X], CPT[X]} *)
ExpandProcess[init_->final_]:=Block[{initL, finalL, X, CPX, TX, CPTX},
    initL = Flatten[{init}];
    finalL = Flatten[{final}];
    X = initL -> finalL;
    CPX = (AntiParticle/@initL)->(AntiParticle/@finalL);
    TX = finalL -> initL;
    CPTX = (AntiParticle/@finalL)->(AntiParticle/@initL);
    Union[{X,CPX,TX,CPTX}]
]

SetAttributes[ProcessEquivalentQ,Orderless];
ProcessEquivalentQ[proc1_Rule,proc2_Rule]:=Block[{},
MemberQ[ExpandProcess[proc1],proc2];
]

LoadConfiguration[name_]:=Block[{file},
    file=ToString[name];
    Evo$Model=Evo$Particles=Evo$POIs=Evo$Processes=False;
    Get[file];

    (* Initialize the Model *)
    If[SameQ[Evo$Model,False],
    Print["No Model is specified through Evo$Model, using SM instead"]; Evo$Model=SM;
    ];
    InitializeModel[Evo$Model];

    (* Check whether the POI is valid *)
    If[!And@@(ModelHasParticle/@Flatten[{Evo$POIs}]),Abort[] ];
    Evo$POIandAntiPOIs=Union[Join[Evo$POIs,AntiParticle/@Evo$POIs] ];

    (* Check whether the process is relevant *)
    If[!And@@(RelevantProcess[#,Evo$POIandAntiPOIs]/@Evo$Processes),Abort[] ];
    Evo$UniProcesses = Union[Evo$Processes,SameTest->ProcessEquivalentQ];

    (* Collect All relevant particles *)
    Evo$Particles = Union[Flatten[Evo$Processes//.{Rule[args1__,args2__]:>Flatten[{args1,args2}], s_ f_[args___]:>f[args]}] ];
]


End[]
