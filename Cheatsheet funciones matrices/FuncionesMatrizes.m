%Metodo de Gauss para la reduccion una matriz a una matriz triangular superior

function [A,b]=Gauss(A,b)
 n=length(b);
 for k=1:n-1
 for i=k+1:n
 m=A(i,k)/A(k,k);
 A(i,k)=0;
 for j=k+1:n
 A(i,j)=A(i,j)-m*A(k,j);
 endfor
 b(i)=b(i)-m*b(k);
 endfor
 endfor
endfunction


%Funcion de resolucion de matriz trialgular superior
function y=backSub(U,b)
n=length(b);
y=zeros(n,1);
for j=n:-1:1
    y(j)=b(j)/U(j,j);
    b(1:j-1)=b(1:j-1)-U(1:j-1,j)*y(j);
endfor
endfunction

function y=backSubMod(U,b)
n=length(b);
y=zeros(n,1);
for j=n:-1:1
    y(j)=b(j)/U(j,j);
    b(1:j-1)=b(1:j-1)-U(1:j-1,j)*y(j);
endfor
k = length(y);
for i=1:k;
  if y(i)<0;
    y(i)=0;
    endif
  endfor
for i=1:k;
  if abs(y(i))<10e-10;
    y(i)=0;
    endif
  endfor
endfunction

function x=backSubInferior(A,B)
  [n n] = size(A);
  A = [A B];
  x=zeros(n,1);
  x(1)=A(1,n+1)/A(1,1);
for i=2:n
  s=0;
  for j=1:n-1
    s = s+A(i,j)*x(j);
  endfor
  x(i)=(A(i,n+1)-s)/A(i,i);
endfor
endfunction

function [U,L] = GaussParcialLU(U,L,k)
  n = length(U(:,1))
  for i = k+1:n
    L(i,k)=U(i,k)/U(k,k);
    U(i,k) = 0;
    for j = k+1:n
      U(i,j)=U(i,j)-L(i,k)*U(k,j);
    endfor
  endfor
endfunction
  
%Funcion para crear una matriz de Vandermonde
function y = vandermonde(v);
  n = length(v);
for k=1:n;
  for i=1:n;
  y(k,i)=(v(k))^(i-1);  
  endfor  
endfor
endfunction

%Funcion para el pivotaje parcial de una matriz
%Pivotaje parcial + resolucion
%pivotaje parcial

function x = GAUSS_ELIM(A, b)
%% Create permutation vector
n = size(A, 1);  % Size of input matrix
r = zeros(n, 1); % Initialize permutation vector
for i = 1 : 1 : n    
   r(i) = i;
endfor
%% Apply Gaussian elimination and rearrange permutation vector
x = zeros(n, 1); % Initialize solution vector
for k = 1 : 1 : n % Go through each element in permutation vector    
    % Compare each element in r(k)th column for the max
    max = abs(A(r(k), r(k)));    
    max_pos = k;    
    for l = k : 1 : n        
        if abs(A(r(l), r(k))) > max            
            max = abs(A(r(l), r(k)));            
            max_pos = l;            
        endif
    endfor
    % Switch the kth r-vector element with max r-vector element
    temp_r = r;
    r(k) = temp_r(max_pos);
    r(max_pos) = temp_r(k);
    % Eliminate A-vector elements in r(k)th column below r(k)th row        
    for i = 1 : 1 : n
        if i ~= k
            zeta = A(r(i), k) / A(r(k), k);
            for j = k : 1 : n
                A(r(i), j) = A(r(i), j) - A(r(k), j) * zeta;                       
            endfor
            b(r(i)) = b(r(i)) - b(r(k)) * zeta;
        endif
    endfor
endfor
% Compute the solution frpm the diagonalized A-matrix
for i = 1 : 1 : n
    x(i) = b(r(i)) / A(r(i), i);
endfor
endfunction

function  [S,b] = separar(A); %Le pasas una matriz ampliada y te retorna las dos matrices separadas
n = length(A);
S = zeros(n-1,n-1);
b = zeros(n-1,1);
for i=1:n-1
  for k=1:n-1
    S(i,k)=A(i,k);   %copiamos todas las columnas menos la ultima
  endfor
endfor
for i=1:n-1
  b(i,1)=A(i,n);
endfor
endfunction

%funcio per realitzar gauss per columna (per implementar maximal)
function [A,b] = GaussC(A,b,k)
  n = length(b);
  for i=k+1:n
    m=A(i,k)/A(k,k);
    A(i,k)=0;
    for j=k+1:n
      A(i,j)=A(i,j)-m*A(k,j);
    endfor
    b(i) = b(i) - m*b(k);
  endfor
endfunction


function [A,b] = Maximal(A,b)
  n = length(b);
  valor = 0;
  
  for j=1:n-1
    for k=j:n
      for i=k:n
        if abs(A(k,k)) < abs(A(i,k)) && abs(A(i,k)) > abs(valor)
          aux = A(j,1:n);
          A(j,1:n) = A(i,1:n);
         A(i,1:n) = aux;
         aux = b(j);
         b(j) = b(i);
          b(i) = aux;
          valor = A(i,k);
        endif
      endfor
   endfor
   for k = j:n
     if abs(A(j,j))<abs(A(j,k))
        aux = A(1:n,j);
        A(1:n,j) = A(1:n,k);
        A(1:n,k) = aux;
      endif
    endfor
    [A,b] = GaussC(A,b,j);
  endfor
endfunction

function b = solucionMarkov(P,i0,n)
T = inv(P);
[A,b] = Gauss(T,i0);
b = backSubMod(A,b);
i=1;
  while i<n
    [A,b] = Gauss(T,i0);
    b = backSubMod(A,b);
    i++
  endwhile
 endfunction

function b = clearError(b,stol)
  n = length(b)
  for i = 1:n;
    if abs(b(i))<stol
      b(i) = 0;
      endif
    endfor
  endfunction

function [L,U]=LUFact(A)
  n=length(A(:,1)); #n files
  m=length(A(1,:)); #n columnes
  if m != n
    printf("bad # of rows and columns");
    return;
  endif
  L=eye(n); #retorna una matriu identitat de n files
  U=A; 
  for k=1:n-1
    for i=k+1:n
      L(i,k)=U(i,k)/U(k,k);
      U(i,k)=0;
      for j=k+1:n
        U(i,j)=U(i,j)-L(i,k)*U(k,j);
      endfor
    endfor
  endfor
endfunction

function [L, U, P] = LU_pivot(A)
    [m, n] = size(A); L=eye(n); P=eye(n); U=A;
    for k=1:m-1
        pivot=max(abs(U(k:m,k)));
        for j=k:m
            if(abs(U(j,k))==pivot);
                ind=j;
                break;
            endif
        endfor
        U([k,ind],k:m)=U([ind,k],k:m);
        L([k,ind],1:k-1)=L([ind,k],1:k-1);
        P([k,ind],:)=P([ind,k],:);
        for j=k+1:m
            L(j,k)=U(j,k)/U(k,k);
            U(j,k:m)=U(j,k:m)-L(j,k)*U(k,k:m);
        endfor
        pause;
    endfor
endfunction

#function [S] = solucionMarkov(P,x0,n)
#  [L,U,P] = LU_pivot(P)
#for i = 0; i < n; i++;
  
#endfor  
#endfunction

function [L U p q] = lucp(A,tol,pm_opt)
%LUCP     LU factorization with complete pivoting.
%
% To compute the LU factorization under default settings:
%
%  [L U p q] = lucp(A)
%
% This produces a factorization such that L*U = A(p,q).  Vectors p and q
% permute the rows and columns, respectively.
%
% The pivot tolerance can be controlled:
%
%  [L U p q] = lucp(A,tol)
%
% The algorithm will terminate if the absolute value of the pivot is less
% than tol.
%
% Permutation matrices can be generated:
%
%  [L U P Q] = lucp(A,tol,'matrix')
%  [L U P Q] = lucp(A,tol,'sparse')
%
% The first will generate full permutation matrices P and Q such that 
% L*U = P*A*Q.  The second generates sparse P and Q.
%
% If A is sparse, L and U will be sparse.
%
% This function works on non-square matrices.
%
% Input:
%  A = matrix
%  tol = pivot tolerance (defualt is 1e-10)
%  pm_opt = permutation output options
%         = 'vector' for permutation vectors, L*U = A(p,q), defualt
%         = 'matrix' for full permutation matrices, L*U = P*A*Q
%         = 'sparse' for sparse permutation matrices, L*U = P*A*Q
%
% Output:
%  L = lower triangular factor
%  U = upper triangular factor
%  p = row permutation
%  q = column permutation
%
% Reference:
%  Algorithm 3.4.2, Matrix Computations, Golub and Van Loan. (3rd ed)
%  
% Other Implementations:
%  Gaussian Elimination using Complete Pivoting.  Alvaro Moraes.
%  http://www.mathworks.com/matlabcentral/fileexchange/25758
%
%  Gauss elimination with complete pivoting.  Nickolas Cheilakos.
%  http://www.mathworks.com/matlabcentral/fileexchange/13451
%  (Does not work with rectangular A)
%
%  Rank Revealing Code.  Leslie Foster.
%  http://www.math.sjsu.edu/~foster/rankrevealingcode.html
%  (Uses mex libraries for computation)
%
%
% 2010-03-28 (nwh) first version.
% 2010-04-14 (nwh) added more references.
% 2010-04-24 (nwh) perform final column swap so U is well conditioned
%
%
% License: see license.txt.
%
% handle optional inputs
if nargin < 2 || isempty(tol)
  tol = 1e-10;
endif
if nargin < 3 || isempty(pm_opt)
  pm_opt = 'vector';
endif
if strcmp(pm_opt,'vector')
  pm_flag = false;
  sp_flag = false;
elseif strcmp(pm_opt,'matrix')
  pm_flag = true;
  sp_flag = false;
elseif strcmp(pm_opt,'sparse')
  pm_flag = true;
  sp_flag = true;
else
  error('lucp:invalidOption','''%s'' is an un recognized option.',pm_opt)
endif
[n m] = size(A);
% pivot vectors
p = (1:n)';
q = (1:m)';
% temp storage
rt = zeros(m,1); % row temp
ct = zeros(n,1); % col temp
t = 0; % scalar temp
for k = 1:min(n-1,m)
  % determine pivot
  [cv ri] = max(abs(A(k:n,k:m)));
  [rv ci] = max(cv);
  rp = ri(ci)+k-1;
  cp = ci+k-1;
  
  % swap row
  t = p(k);
  p(k) = p(rp);
  p(rp) = t;
  rt = A(k,:);
  A(k,:) = A(rp,:);
  A(rp,:) = rt;
  
  % swap col
  t = q(k);
  q(k) = q(cp);
  q(cp) = t;
  ct = A(:,k);
  A(:,k) = A(:,cp);
  A(:,cp) = ct;
  
  if abs(A(k,k)) >= tol
    rows = (k+1):n;
    cols = (k+1):m;
    A(rows,k) = A(rows,k)/A(k,k);
    A(rows,cols) = A(rows,cols)-A(rows,k)*A(k,cols);
  else
    % stop factorization if pivot is too small
    break
  endif
endfor
% final column swap if m > n
if m > n
  % determine col pivot
  [cv ci] = max(abs(A(n,n:m)));
  cp = ci+n-1;
  
  % swap col
  t = q(n);
  q(n) = q(cp);
  q(cp) = t;
  ct = A(:,n);
  A(:,n) = A(:,cp);
  A(:,cp) = ct;
endif
% produce L and U matrices
% these are sparse if L and U are sparse
l = min(n,m);
L = tril(A(1:n,1:l),-1) + speye(n,l);
U = triu(A(1:l,1:m));
% produce sparse permutation matrices if desired
if pm_flag
  p = sparse(1:n,p,1);
  q = sparse(q,1:m,1);
endif
% produce full permutation matrices if desired
if ~sp_flag
  p = full(p);
  q = full(q);
  endif
endfunction