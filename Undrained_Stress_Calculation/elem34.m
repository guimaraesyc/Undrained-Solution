% This fuction assembles the triangular local nodes

function [Tri,count] = elem34(Tri,i,j,m,n,count)

% type 3
elem=count;
    
    
    %------------%
    %    type3   %
    %------------%
    
    % local node 1
    temp = (2*i - 1) + (2*j)*(2*m+1);
    Tri(elem).X(1) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(1) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,1) = temp;
    % local node 2
    temp = (2*i - 1) + 2*(j-1)*(2*m+1);
    Tri(elem).X(2) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(2) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,2) = temp;
    % local node 3
    temp = (2*i + 1) + (2*j)*(2*m+1);
    Tri(elem).X(3) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(3) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,3) = temp;
    % local node 4
    temp = (2*i - 1) + (2*j-1)*(2*m+1);
    Tri(elem).X(4) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(4) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,4) = temp;
    % local node 5
    temp = 2*i + (2*j)*(2*m+1);
    Tri(elem).X(5) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(5) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,5) = temp;
    % local node 6
    temp = 2*i + (2*j-1)*(2*m+1);
    Tri(elem).X(4) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(4) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,6) = temp;
    
    
    %------------%
    %    type4   %
    %------------%
    elem = elem + 1;
    count=elem;
    
    % local node 1
    temp = (2*i + 1) + 2*(j-1)*(2*m+1);
    Tri(elem).X(1) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(1) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,1) = temp;
    % local node 2
    temp = (2*i + 1) + (2*j)*(2*m+1);
    Tri(elem).X(2) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(2) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,2) = temp;
    % local node 3
    temp = (2*i - 1) + 2*(j-1)*(2*m+1);
    Tri(elem).X(3) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(3) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,3) = temp;
    % local node 4
    temp = (2*i + 1) + (2*j-1)*(2*m+1);
    Tri(elem).X(4) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(4) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,4) = temp;
    % local node 5
    temp = 2*i + 2*(j-1)*(2*m+1);
    Tri(elem).X(5) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(5) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,5) = temp;
    % local node 6
    temp = 2*i + (2*j-1)*(2*m+1);
    Tri(elem).X(6) = Tri(1).gnodes_0(2*temp-1);
    Tri(elem).Y(6) = Tri(1).gnodes_0(2*temp);
    Tri(1).ConM(elem,6) = temp;
    
end
