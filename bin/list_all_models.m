function [models, magnitudes, is_model_identity] = list_all_models( G, G_fake, labels, threshold )
  %{
  PURPOSE:
  List all models derivable from a matrix G.

  INPUT:
  G - 
  labels - 
  threshold - I usually pick like 0.1

  OUTPUT:
  
  %}

  [~,nl] = size(G); %nl = num_library
  is_model_identity = zeros(nl, 1);
  models = cell(nl,1);
  
  for i=1:nl      
    [c, models{i}] = list_models(G, labels, threshold);
    magnitudes(i) = vecnorm( G*c );
    
    if numel(G_fake) ~= 1
      eta_fake = norm(G_fake*c)/max( vecnorm(G_fake.*c') );
      if( eta_fake < 3e-1 )
        is_model_identity(i) = 1;
      end
    end
    
    nonzero = find(c);
    
    %eliminate data from G and labels now
    labels(nonzero(end)) = [];
    G(:,nonzero(end)) = [];
    if numel(G_fake) ~= 1
      G_fake(:,nonzero(end)) = [];
    end
    
    models{i} = is_model_identity(i) + "&" + models{i};
  end
end