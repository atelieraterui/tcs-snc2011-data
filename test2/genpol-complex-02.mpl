with (RandomTools);

genpolnonmonic_complex := proc (deg, var, coefrange)
    return Generate(polynom(complex(rational(range=coefrange)),
                            var, degree=deg)) + var ^ deg
end proc;

genpolmonic_complex := proc (deg, var, coefrange)
    return genpolnonmonic_complex(deg - 1, var, coefrange) + var ^ deg
end proc;

gentestpol_complex := proc (degf, degg, degd, var, perturbf, perturbg,
                            svdthreshold, numpol)

    # complex version

    # degf: degree of f
    # degg: degree of g
    # degd: degree of gcd(f,g)
    # var: variable
    # perturbf: norm of perturbation of f in 2-norm
    # perturbg: norm of perturbation of g in 2-norm
    # svdthreshold: threshold of
    #   SingularValues(subresmat(f[i], g[i], degd))[degf + degg - 2*degd]
    # numpol: the number of test polynomials
    local i, d, u, v, ff, gg, fnoize0, gnoize0, fnoize, gnoize, f, g, 
    singularvalue;
    d := Array(1 .. numpol);
    u := Array(1 .. numpol);
    v := Array(1 .. numpol);
    ff := Array(1 .. numpol);
    gg := Array(1 .. numpol);
    fnoize0 := Array(1 .. numpol);
    gnoize0 := Array(1 .. numpol);
    fnoize := Array(1 .. numpol);
    gnoize := Array(1 .. numpol);
    f := Array(1 .. numpol);
    g := Array(1 .. numpol);
    for i to numpol do
        singularvalue := 0;
        while (singularvalue < svdthreshold) do
            d[i] := evalf(genpolmonic_complex(degd, var, -10 .. 10));
            u[i] := evalf(genpolmonic_complex(degf - degd, var, -10 .. 10));
            v[i] := evalf(genpolmonic_complex(degg - degd, var, -10 .. 10));
            ff[i] := expand(d[i] * u[i]);
            gg[i] := expand(d[i] * v[i]);
            fnoize0[i] := evalf(genpolnonmonic_complex(degf - 1, var, -10 .. 10));
            gnoize0[i] := evalf(genpolnonmonic_complex(degg - 1, var, -10 .. 10));
            fnoize[i] := (perturbf * fnoize0[i]) / polynorm(fnoize0[i], var, 2);
            gnoize[i] := (perturbg * gnoize0[i]) / polynorm(gnoize0[i], var, 2);
            f[i] := ff[i] + fnoize[i];
            g[i] := gg[i] + gnoize[i];
            singularvalue :=
            SingularValues(subresmat(f[i], g[i], var, degd))[degf + degg - 2*degd];
        end;
    end:
    return f, g, d, u, v, fnoize, gnoize;
end proc;
