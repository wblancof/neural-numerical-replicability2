from SimulationInitialization import zero, int2, int6


def rk4 (dydt, tspan, y0, m, y):
    f0 = [0.0, 0.0, 0.0, 0.0]
    f1 = [0.0, 0.0, 0.0, 0.0]
    f2 = [0.0, 0.0, 0.0, 0.0]
    f3 = [0.0, 0.0, 0.0, 0.0]
    u1 = [0.0, 0.0, 0.0, 0.0]
    u2 = [0.0, 0.0, 0.0, 0.0]
    u3 = [0.0, 0.0, 0.0, 0.0]
    dt = tspan[1] - tspan[0]
    t = tspan[0]

    dydt(t, y0, f0)

    for i in range(m):
        u1[i] = y0[i] + dt * f0[i] / int2
    dydt(t + dt / int2, u1, f1)

    for i in range(m):
        u2[i] = y0[i] + dt * f1[i] / int2
    dydt(t + dt / int2, u2, f2)

    for i in range(m):
        u3[i] = y0[i] + dt * f2[i]
    dydt(t + dt, u3, f3)

    for i in range(m):
        y[i] = y0[i] + dt * (f0[i] + int2 * f1[i] + int2 * f2[i] + f3[i]) / int6
