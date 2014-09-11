# Test function
f = (t, y, p1, p2) -> [
    p1*t*y[3]*y[0]
    p2*t*y[3]*y[0].pow(5)
    2*t*y[3]
    -2*t*(y[2]-1)
]

# RK
{rk, ode} = $blab.ode # Click to see imported functions
ts = linspace 0, 1, 100 #; Time grid (change step size here)
y0 = [1, 1, 1, 1] # Initial condition
o = 2 # RK order (Choose 1-4)
w = ode(rk[o], f, ts, y0, 2, 10) #; Solve

# Ideal solution (p1=2, p2=10)
z = ([exp(sin(t*t)), exp(5*sin(t*t)), sin(t*t)+1, cos(t*t)] for t in ts) #;

#  Error
plot ts, (w-z).T,
    xlabel: "t"
    ylabel: "Error"
    height: 170
    series:
        shadowSize: 0
        color: "black"
        lines: lineWidth: 1
