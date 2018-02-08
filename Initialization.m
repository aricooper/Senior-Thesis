function [gco,gcn,names,G, cs1] = Initialization(handles)
% generates metric of geodesics for metric met()
computing = true;
%% Generating components

% runmetric

% [x,bco,bcn,gco,gcn,cs1,cs2] = runmetric(met);

[X,names] = met(handles);
% Differentiate the radius vector to get the
% contravariant base vectors.
% bco=simplify([diff(X,names(1)),diff(X,names(2)),diff(X,names(3)),diff(X,names(4))]);

% Obtain the metric tensor components as dot products of
% the base vectors
% gco=simplify(bco.'*bco); gcn=simplify(bcn.'*bcn);
bco=simplify([diff(X,names(1)),diff(X,names(2)),diff(X,names(3)),diff(X,names(4))]);

M = sym(diag(X));  % since metric is diagonal
gco = M;
gcn = inv(gco);

if(computing)
% Compute the Christoffel symbols.
cs1=sym(zeros(4,4,4));
  % Obtain symbols of the first kind by differentiating
  % the covariant metric tensor components.    
  for k=1:4
    for i=1:4, for j=1:i
      cs1(i,j,k)=1/2*gcn(k,k)*(diff(gco(j,k),names(i))+...
             diff(gco(i,k),names(j)) - diff(gco(i,j),names(k))); 
      if j~=i, cs1(j,i,k)=cs1(i,j,k); end   
    end, end
  end
  
%% Generate Geodesics

% initialize four velocity u
u = sym(zeros(1,4));

temp = sym('u',[1 4]);

for i=1:4
    for j=1:4
eqns(i,j) = temp(i)*gco(i,j)*temp(j);
tempsol = solve(sum(eqns(i,:)) == -1, temp(i));
if isempty(tempsol)
    tempsol = 0;
end
sol(i) = tempsol(1);

    end
    
end

for i=1:4
    u(i) = simplify(sol(i));
end
u(2:4) = 0;     % for object at rest (only time component of u)


% Solve Geodesic equations


G = geodesics(cs1, u, handles);


end












