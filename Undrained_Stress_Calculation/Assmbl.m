% It assembles the global stiffness matrix (linear)

function Kgel = Assmbl(delU,Tri,dN)
format long
  
  Kgel = zeros(2*Tri(1).nnodes);        % degrees of freedom (2 DOF x node)
  K = zeros(12);
    
  % Global stiffness matrix
  for count=1:Tri(1).nel
      %K = CalcStiff(Tri,dN,count);
      K = Calcstiff_WI(Tri,dN,count);
      %K = CalcStiff_NL(delU,Tri,dN,count);
    for node_i=1:6
      iglobal = Tri(1).ConM(count,node_i);
      for node_j=1:6
      jglobal = Tri(1).ConM(count,node_j);
      % global element stiffness matrix
      Kgel((2*iglobal-1),(2*jglobal-1)) = K((2*node_i-1),(2*node_j-1)) + Kgel((2*iglobal-1),(2*jglobal-1));
      Kgel((2*iglobal-1),(2*jglobal)) = K((2*node_i-1),(2*node_j)) + Kgel((2*iglobal-1),(2*jglobal));
      Kgel(2*iglobal,(2*jglobal-1)) = K(2*node_i,(2*node_j-1)) + Kgel(2*iglobal,(2*jglobal-1));
      Kgel(2*iglobal,2*jglobal) = K(2*node_i,2*node_j) + Kgel(2*iglobal,2*jglobal);
      end
    end
  end
end