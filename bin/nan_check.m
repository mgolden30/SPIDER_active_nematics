function nan_check( data, name )
  %{
  This function checks an array or single value to see if it a NAN.
  %}
  if sum(isnan(data), 'all') ~= 0
    error('Error: nan_check has found that %s has NANs in it.\n', name)
  end
end