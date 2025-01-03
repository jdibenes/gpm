
function t = twp_compute_t(PA1, PB1, PA2, PB2, R, options)
select_t = options.select_t;

switch (select_t)
case 1,    t = PA2 - R*PA1;
case 2,    t = PB2 - R*PB1;
otherwise, t = ((PA2 + PB2) - R*(PA1 + PB1)) / 2; % least squares
end
end
