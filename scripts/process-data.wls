#!/usr/bin/env wolframscript

$K=3;$N=8;


radius[tqpData_]:={tqpData//First,1/2 Total@(#^2&/@tqpData[[2;;(Length[tqpData]-1)/2]])};


angMom[tqpData_]:=Module[{q, p}, (
q = tqpData[[2;;(Length[tqpData]+1)/2]];
p =tqpData[[1+(Length[tqpData]+1)/2;;]];
q =ArrayReshape[q, {$K, $N^2-1}];
p =ArrayReshape[p, {$K, $N^2-1}];

{tqpData//First, Table[1/2 Sum[q[[i, a]]p[[j, a]]- p[[i, a]]q[[j, a]], {a, 1, $N^2-1}], {i, 1, $K}, {j, 1, $K}]}
)];


comm[x1_, x2_, $f_]:=Table[Sum[$f[[a, b, c]]*(x1[[b]]x2[[c]]-x1[[c]]x2[[b]]), {b, 1, $N^2-1}, {c, b+1, $N^2-1}], {a, 1, $N^2-1}];


energy[tqpData_]:=Module[{q, p, cc}, (
q = tqpData[[2;;(Length[tqpData]+1)/2]];
p =tqpData[[1+(Length[tqpData]+1)/2;;]];

q =ArrayReshape[q, {$K, $N^2-1}];

{tqpData//First, 1/4 Total@(#^2&/@p) + 1/4 Sum[Total@(#^2&/@comm[q[[i]], q[[j]], $f]), {i, 1, $K}, {j, i+1, $K}]}
)];


process[filename_]:=Module[{data, d, e, \[ScriptCapitalD], msg}, (
data = Import[filename, "Table"];

e = energy[data//First]//Last;

d = (radius/@data)[[All, 2]];

\[ScriptCapitalD] = SmoothKernelDistribution[d];

Prepend[#[\[ScriptCapitalD]]&/@{Mean, StandardDeviation, Skewness}, e]
)];


$f =Import["./F/N"<>ToString[$N]<>".dat", "Table"]//Flatten//ArrayReshape[#, {$N^2-1, $N^2-1, $N^2-1}]&;

insert[data_] := (StringJoin @@ (ToString[#] <> "    " & /@ data))<>"\n";

(* argv = Rest @ $ScriptCommandLine;

Print[insert@process[argv[[1]]]]; *)

files = FileNames["*.dat", "./dump/"]

stats = OpenAppend["./stats/stats_nolog.dat"];

Do[Print["a="<>ToString[a]]; WriteString[stats, insert@process[files[[a]]]], {a, 1, Length[files]}]

