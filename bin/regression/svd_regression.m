function [c, residual, residual_sequence, singular_vals] = svd_regression( A, threshold )
  %{
  PURPOSE:
  The purpose of this function is to find the sparse vector c that
  minimizes a dimensionless residual.

  eta = ||Ac|| / ||c||, where ||c|| is the 2-norm

  ASSUMPTIONS:
  I will assume A(:,1) is the integrated time derivative as a convention.

  INPUT:
  A - a skinny matrix of integrated data. Columns are a certain library
  term being integrated over some huge number of domains.
  threshold - a number between 0 and 1. 
              Determines if a term is worth throwing out in sparsification.
              For example threshold = 0.2 throws out terms if they
              contribute less than 20% to the solution.

  OUTPUT:
  c - a sparse coefficient matrix, which hopefully results in 
      a small residual of norm(A*c)
  residual  - eta as defined above. This tells you how good your sparse model is.
              lower is better
  %}

  %nw is the number of windows that have been integrated
  %nl is the number of library terms (including the time derivative)
  [nw, nl] = size(A);
  
  if( nw < nl )
     fprintf('Error: your matrix you want to perform SVD regression on is possibly too small.\n')
     fprintf('number of windows = %d\n', nw);
     fprintf('number of library terms = %d\n', nl);
     fprintf('Taking a transpose *might* fix your problem, but no promises.\n')
  end
  
  %First compute the svd
  [U, S, V] = svd(A);
  singular_vals = zeros(nl,1);
  for i=1:nl
    singular_vals(i) = S(i,i); 
  end
  
  c = V(:,end); %last column of V is the most singular vector.
  
  %This section computes residual_sequence, which is the residual when
  %terms are turned off one-by-one.
  residual_sequence = norm(A*c)/max( vecnorm(A.*c') );
  
  %B will be A that we remove columns from and recompute svd
  B = A;
  terms = 1:nl;
  %Terms allows me to easilly reconstruct coefficent vector
  
  last_residual = 1e9; %start with an absurdly high starting value.
  c_best = c;
  for i=2:nl %Remove a term this many times
    [~, nl_reduced] = size(B);
    
    residuals = zeros(1, nl_reduced);
    for j=1:nl_reduced
      c_temp = c;
      c_temp(j) = 0;
      residuals(j) = norm(B*c_temp)/max( vecnorm(B.*c_temp') );
      %residuals(j) = norm(B*c_temp)/norm(c_temp);
    end
    %Now that we know the residual when each term is turned off, choose the
    %minimum
    [min_val, idx] = min( residuals );
    
    %Delete this column from B since it isn't needed.
    %Matlab has a useful syntax to do this
    B(:,idx)   = [];
    terms(idx) = [];
    
    %Perform SVD again on this reduced matrix and update our guess c.
    [U, S, V] = svd(B);
    c = V(:,end);
    
    res = norm(B*c)/max( vecnorm(B.*c') );
    %res = norm(B*c)/norm(c);
    
    if res/last_residual < 1 + threshold
      c_best = zeros(nl, 1);
      c_best( terms ) = c;
    end
    
    last_residual = res;
    residual_sequence = [residual_sequence res];
  end
  
  c = c_best;
  
  %output sparse residual
  residual = norm(A*c)/max( vecnorm(A.*c') );
  
  %Now search the orthogonal subspace for single term laws.  
  for i=1:nl
    if( c(i) == 0 )
      %This term is in the orthogonal space
      single_term = zeros(nl,1);
      single_term(i) = 1;
      single_term_residual = norm(A*single_term)/norm(A*c);
      
      if single_term_residual < 1
        c = single_term;
        %Use the single term residual
        residual = norm(A*c)/mean( vecnorm(A) );
      end
      
      %fprintf('i = %d: single term residual = %f\n', i, single_term_residual );
    end
  end
  
  %Normalize c(1)=1 if nonzero
  if(c(1) ~= 0)
    c = c/c(1); 
  end
  
  %fix sum of coefficients > 0
  if( sum(c) < 0 )
    c = -c; 
  end
end