function get_y_states(x,y_stat,log_y_spread)

y_low = x
y_high = y_low + log_y_spread
y = [y_low y_high]
mean_y = exp.(y)*y_stat'

crit = (mean_y[1] - 1)^2

return crit, y_high, y_stat

end
