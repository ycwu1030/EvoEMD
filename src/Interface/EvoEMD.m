
BeginPackage["EvoEMD`", {"FeynArts`","FormCalc`","Rubi`"}]

LoadConfiguration::usage =
"LoadConfiguration[name], loading configuration file for EvoEMD"

Evo$Model::usage =
"The model to be used with EvoEMD"

Evo$POIs::usage =
"The list of particle of interests supplied by user to construct the Boltzmann equation"

Evo$POIandAntiPOIs::usage =
"The list of POIs and their anti-particles"

Evo$Processes::usage =
"The list of processes supplied by user for construct the collision rate in the Boltzmann equation"

Evo$UniProcesses::usage =
"Similar to Evo$Processes but with CP/T/CPT equivalent one removed"

Evo$Particles::usage =
"Particles involved in Evo$Processes"


$EvoEMD = {0,0,0}

$EvoEMDVersion = "EvoEMD 0.0.0 (17 Nov 2021)"

$EvoEMDDir = DirectoryName[
    $InputFileName /. HoldPattern[$InputFileName] :>
     (File /. FileInformation[System`Private`FindFile[$Input] ] )
]

Print[""];
Print[$EvoEMDVersion];
Print["by Yongcheng Wu"];

(* Needs["FeynArts`"];
Needs["FormCalc`"];
Needs["Rubi`"]; *)
$FAVerbose=0;
$FCVerbose=0;

Block[{$Path=$EvoEMDDir},
<<EvoInit`;
<<EvoProcess`;
]

EndPackage[]
