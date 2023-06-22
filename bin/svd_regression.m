function [c, residual, residual_sequence] = svd_regression( A, threshold )
  %{
  PURPOSE:
  This function takes a matrix of integrated library terms and finds the
  best sparse model. Here by "model" I mean a vector c such that c
  is sparse (physically interpretable) and norm(A*c)≈0 in some sense.

  INPUT:
  A         - a skinny matrix of integrated library terms
  threshold - a sparsification parameter in [0,1]. A model will be
              sparsified if 

  OUTPUT:
  c         - the identified sparse coefficient vector
  residual  - a dimensionless measure of the accuracy of the sparse model
              Formally it is the ratio of the residual to the largest term
              in the model.
              
              residual = eta = || Ac || / max(vecnorm(c.*A))
  
  residual_sequence - a typical example is [0.3566, 0.3611, 0.4934, 1]
                      this sequence shows how the residual changes with
                      sparsification. 0.3566 is the residual of the dense
                      c, and it increases as we set terms to 0.
  %}

  assert( threshold > 0 );
  assert( threshold < 1 );

  %nw is the number of windows that have been integrated
  %nl is the number of library terms
  [nw, nl] = size(A);
  assert( nw > nl );
  %^You should be using an overdetermined system 
  
  [U, S, V] = svd(A);
  %last column of V is the most singular vector, which gives the smallest
  %norm(A*c). This vector will not be sparse in general, so sparsification
  %is next.
  c = V(:,end);

  residual_sequence    = zeros(nl, 1);
  residual_sequence(1) = norm(A*c)/max( vecnorm(A.*c') );
  
  B = A;        %B will be A that we remove columns from and recompute svd
  terms = 1:nl; %terms allows me to easilly reconstruct coefficent vector
  
  %Use the actual most recent residual.
  last_residual = norm(B*c)/max( vecnorm(B.*c') );
  
  %We are about to overwrite c, so store it as "c_best" for now.
  c_best = c; 
  for i=2:nl
    [~, nl_reduced] = size(B);
    
    test_residuals = zeros(1, nl_reduced);
    for j=1:nl_reduced
      c_temp = c;
      c_temp(j) = 0;
      test_residuals(j) = norm(B*c_temp)/max( vecnorm(B.*c_temp') );
    end
    
    %Choose to turn off the term that produces the lowest residual
    [min_val, idx] = min( test_residuals );
    
    %Delete this column from B since it isn't needed.
    B(:,idx)   = [];
    terms(idx) = [];
    
    %Perform SVD again on this reduced matrix and update our guess c.
    [U, S, V] = svd(B);
    c = V(:,end);
    
    new_residual = norm(B*c)/max( vecnorm(B.*c') );
    
    %If this sparsification results in a new residual within our threshold,
    %
    if new_residual/last_residual < 1 + threshold
      c_best = zeros(nl, 1);
      c_best( terms ) = c;
    end
    
    last_residual = new_residual;
    residual_sequence(i) = new_residual;
  end
  
  %recover our best c from c_best
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
         %residual = norm(A*c)/mean( vecnorm(A) );
       end
       
       %fprintf('i = %d: single term residual = %f\n', i, single_term_residual );
     end
   end

  %Normalize c(1)=1 if nonzero
  %{
  if(c(1) ~= 0)
    c = c/c(1); 
  end
  %}
  
  %fix the largest coefficient to be positive
  [~, idx] = max( abs(c) );
  if( c(idx) < 0 )
    c = -c; 
  end
end