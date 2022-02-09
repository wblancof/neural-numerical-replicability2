
template <class actualRealType>

void rk4 (void dydt(actualRealType t, const std::array<actualRealType, 4> u, std::array<actualRealType, 4> &f), actualRealType *tspan, std::array<actualRealType, 4> y0, int m, std::array<actualRealType, 4> &y) {
    std::array<actualRealType, 4> f0;
    std::array<actualRealType, 4> f1;
    std::array<actualRealType, 4> f2;
    std::array<actualRealType, 4> f3;
    std::array<actualRealType, 4> u1;
    std::array<actualRealType, 4> u2;
    std::array<actualRealType, 4> u3;
    actualRealType dt = tspan[1] - tspan[0];
    actualRealType t = tspan[0];

    dydt(t, y0, f0);

    for(int i = 0; i < m; i++)
        u1[i] = y0[i] + dt * f0[i] / int2;
    dydt(t + dt / int2, u1, f1);

    for(int i = 0; i < m; i++)
        u2[i] = y0[i] + dt * f1[i] / int2;
    dydt(t + dt / int2, u2, f2);

    for(int i = 0; i < m; i++)
        u3[i] = y0[i] + dt * f2[i];
    dydt(t + dt, u3, f3);

    for(int i = 0; i < m; i++)
        y[i] = y0[i] + dt * (f0[i] + int2 * f1[i] + int2 * f2[i] + f3[i]) / int6;
}
