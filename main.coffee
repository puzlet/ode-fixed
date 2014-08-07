# rk methods
# (tn, tl, yl) = (time-now, time-last, output-last)

rk = {}

# Order 1 (Forward Euler)
rk[1] = (tn, tl, yl) ->
    h = tn - tl
    yl + h*(feval tl, yl)

# Order 2 (Improved Euler)
rk[2] = (tn, tl, yl) ->
    h = tn - tl
    F1 = feval tl, yl
    F2 = feval tl+h, yl+h*F1
    yl + (h/2)*(F1 + F2)

# Order 3 (Bogacki-Shampine) 
rk[3] = (tn, tl, yl) ->
    h = tn - tl
    F1 = feval tl, yl
    F2 = feval tl+0.5*h, yl+0.5*h*F1
    F3 = feval tl+0.75*h, yl+0.75*h*F2
    yl + (h/9)*(2*F1 + 3*F2 + 4*F3)

# Order 4 (Classical Runge Kutta)
rk[4] = (tn, tl, yl) ->
    h = tn - tl
    F1 = feval tl, yl
    F2 = feval tl+0.5*h, yl+0.5*h*F1
    F3 = feval tl+0.5*h, yl+0.5*h*F2
    F4 = feval tn, yl+h*F3
    yl + (h/6)*(F1 + 2*F2 + 2*F3 + F4)

# Helpers
zeros = (m,n) -> ((0 for [1..n]) for [1..m])
feval = (t, y) ->  # overridden in ode to allow function params

# Time-step
ode = (rk, f, ts, y0, params...) ->

    feval = (ti, yi) ->
        f(ti, yi, params...)

    L = ts.length # number of time steps
    N = y0.length # number of system states
    
    Y = zeros(L,N)
    Y[0] = y0
    Y[i] = rk(ts[i], ts[i-1], Y[i-1]) for i in [1...L]
    Y

# Export
$blab.ode = {rk, ode}

#!end (coffee)

