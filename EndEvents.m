function [value,isterminal,direction]=EndEvents(t,var) % Stop the integration at a certain altitude
rMars=3396190;              % Mars equatorial radius [m]
value = (var(1)-rMars-10000);           % i.e. stop at 10km altitude
isterminal=1;   % Stop the integration
direction=0;    % All end events from either direction
end
