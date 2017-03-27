OnBreak := QUIT_GAP;

LoadPackage("grape");

OnHjelmslev := function(pt, g)
  pt := List(pt);
  pt[2] := OnPoints(pt[2]+1, g)-1;
  return pt;
end;

Color := function(prefixes, names)
  local prefix, name, result, names_with_prefix, i;
  result := [];
  for prefix in prefixes do
    names_with_prefix := [];
    for i in [1..Length(names)] do
      if names[i]{[1..Length(prefix)]} = prefix then
        AddSet(names_with_prefix, i);
      fi;
    od;
    Add(result, names_with_prefix);
  od;
  return result;
end;

HjelmslevPlanes := function(q, vertices, edges)
local G, v, vertex, edge, i, level, Incidence, vertices_level, edges_level, planes_level;
level := Length(vertices[1]) - 1;
v := q^2+q+1;
G := Group(PermList(Concatenation([2..v], [1])));
vertices_level := [];
edges_level := [];
planes_level := [];
for i in [1..level] do
  Add(vertices_level, []);
  Add(edges_level, []);
od;
for vertex in vertices do
  for i in [1..level] do
    AddSet(vertices_level[i], vertex{[1..i+1]});
  od;
od;
for edge in edges do
  for i in [1..level] do
    AddSet(edges_level[i], Immutable([edge[1]{[1..i+1]}, edge[2]{[1..i+1]}]));
  od;
od;
for i in [1..level] do
  Incidence := function(pta, ptb)
    if pta[1] = ptb[1] then
      return false;
    fi;
    pta := List(pta);
    ptb := List(ptb);
    if pta[1] = "p" then
      ptb[2] := (ptb[2] - pta[2]) mod v;
      pta[2] := 0;
      return [pta, ptb] in edges_level[i];
    else
      pta[2] := (pta[2] - ptb[2]) mod v;
      ptb[2] := 0;
      return [ptb, pta] in edges_level[i];
    fi;
  end;
  Add(planes_level, GRAPE_Graph(G, vertices_level[i], OnHjelmslev, Incidence));
od;
return planes_level;
end;

HjelmslevColors := function(planes_level)
local i, colors_level;
colors_level := [];
Add(colors_level, Color([["p"],["l"]], planes_level[1].names));
for i in [2..Length(planes_level)] do
  Add(colors_level, Color(planes_level[i-1].names, planes_level[i].names));
od;
return colors_level;
end;

HjelmslevTypeColors := function(planes_level)
local i, colors_level;
colors_level := [];
for i in [1..Length(planes_level)] do
  Add(colors_level, Color([["p"],["l"]], planes_level[i].names));
od;
return colors_level;
end;

AutKernel := function(planes, colors, level)
return AutGroupGraph(planes[level], colors[level]);
end;

DecodePermutation := function(planes, level, enc_perm)
local result, pair;
result := ListWithIdenticalEntries(Length(enc_perm), 0);
for pair in enc_perm do
  result[Position(planes[level].names, pair[1])] := Position(planes[level].names, pair[2]);
od;
return PermList(result);
end;

LiftAut := function(planes, colors, level, perm)
local original_colors, permuted_colors;
original_colors := colors[level];
permuted_colors := Permuted(original_colors, perm);
return GraphIsomorphism(rec(graph := planes[level], colourClasses := original_colors),
                        rec(graph := planes[level], colourClasses := permuted_colors));
end;

LiftSubgroupOfCyclicGroup := function(planes, colors, level, generator)
  local i, p, e, ee, order, fac, ps, es, lift, this_lift, lift_of;
  order := Order(generator);
  fac := PrimePowersInt(order);
  ps := fac{[1,3..Length(fac)-1]};
  es := fac{[2,4..Length(fac)]};
  lift := ();
  lift_of := ();
  for i in [1..Length(ps)] do
    p := ps[i];
    e := es[i];
    for ee in [e,e-1..1] do
      Print("Trying to lift ", generator^(order/(p^ee)), " of order ", p^ee," ...\n");
      this_lift := LiftAut(planes, colors, level, generator^(order/(p^ee)));
      if this_lift <> fail then
        Print("Success.\n");
        lift := lift * this_lift;
        lift_of := lift_of * generator^(order/(p^ee));
        break;
      fi;
    od;
  od;
  return [lift, lift_of];
end;

LiftSubgroupOfAbelianGroup := function(planes, colors, level, gp)
  local lifted_gens, gens, llo, gen;
  lifted_gens := [()];
  gens := [()];
  for gen in GeneratorsOfGroup(gp) do
    llo := LiftSubgroupOfCyclicGroup(planes, colors, level, gen);
    if llo[1] <> () then
      Add(lifted_gens, llo[1]);
      Add(gens, llo[2]);
    fi;
  od;
  return [Group(lifted_gens), Group(gens)];
end;
