function [Q_11out, Q_12out] =  nematic_normalize(Q_11, Q_12)
  norm_q = 2*sqrt( Q_11.*Q_11 + Q_12.*Q_12 );
  Q_11out = Q_11./norm_q;
  Q_12out = Q_12./norm_q;
end
