planes := HjelmslevPlanes(q, vertices, edges);
colors := HjelmslevColors(planes);

level := Length(vertices[1]) - 1;

k12 := AutKernel(planes, colors, 2);
Print("K_1^2 has order ", Order(k12), "\n");
if level > 2 then
  Assert(0, IsAbelian(k12), "Non-abelian kernel not implemented");
  k13lift := LiftSubgroupOfAbelianGroup(planes, colors, 3, k12);
  Print("K_1^3(2) has order ", Order(k13lift[2]), ".\n");
  if everything then
    k23 := AutKernel(planes, colors, 3);
    Print("K_2^3 has order ", Order(k23), ".\n");
    gens := Concatenation(GeneratorsOfGroup(k13lift[1]), GeneratorsOfGroup(k23));
    k13 := Group(gens);
    Print("K_1^3 has order ", Order(k13), ".\n");
  fi;
fi;

g := DecodePermutation(planes, 1, gencoded);

llo := LiftSubgroupOfCyclicGroup(planes, colors, 2, g);
l := llo[1];
lo := llo[2];
Print("N/C^2(1) has order ", Order(lo), ".\n");
if everything then
  ngens := Concatenation(GeneratorsOfGroup(k12), [l]);
  N := Normalizer(Group(ngens), planes[2].group);
  Print("N/C^2 has order ", Order(N), ".\n");
fi;

