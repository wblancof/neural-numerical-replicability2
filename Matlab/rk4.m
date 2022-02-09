

function y  = rk4 (dydt, tspan, y0, n)

  y0 = y0(:);

  m = size ( y0, 1 );
  t = zeros ( n + 1, 1 );
  y_temp = zeros ( n + 1, m );

  tfirst = tspan(1);
  tlast = tspan(2);
  dt = ( tlast - tfirst ) / n;

  t(1,1) = tspan(1);
  y_temp(1,:) = y0(:);
  
  for i = 1 : n
    f1 = dydt( y_temp(i,:) );
    f2 = dydt( y_temp(i,:) + dt * f1 / 2.0 );
    f3 = dydt( y_temp(i,:) + dt * f2 / 2.0 );
    f4 = dydt( y_temp(i,:) + dt * f3 );

    t(i+1,1) = t(i,1) + dt;
    y_temp(i+1,:) = y_temp(i,:) + dt * ( f1 + 2.0 * f2 + 2.0 * f3 + f4 ) / 6.0;       
  end
  y = y_temp';
  
  return
end