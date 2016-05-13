test_uvgcd := proc(f, g, var, givendeg, initialtheta, case_begin, case_end)

    # test_uvgcd used with uvgcd

    # f: array of polynomials
    # g: array of polynomials
    # var: main variable of polynomials (dummy argument)
    # theta: torrelance given to the algorithm
    # case_begin: index of the beginning of the test cases
    # case_end: index of the end of the test cases
    #
    local i, j, numofcases, st, ddd, u, v, eps,
    perturb, perturb2, ppp, qqq, anew, bnew, vv, deg,
    presenttheta, lasttheta, thetamin, thetamax, lastdeg, bisectional,
    comptime_table, perturb_table, perturb2_table, gcd_table,
    pol_table, deg_table, theta_table;
    comptime_table := Array(case_begin .. case_end);
    perturb_table := Array(case_begin .. case_end);
    perturb2_table := Array(case_begin .. case_end);
    gcd_table := Array(case_begin .. case_end);
    pol_table := Array(case_begin .. case_end);
    deg_table := Array(case_begin .. case_end);
    theta_table := Array(case_begin .. case_end);
    printf("Initial torrelance: %f\n", initialtheta);
    for i from case_begin to case_end do
        printf("Case %d:\n", i);
#        presenttheta := (thetamin + thetamax) / 2.0;
        presenttheta := initialtheta;
        lasttheta := initialtheta;
        thetamin := initialtheta;
        thetamax := initialtheta;
        st := time():
        ddd, u, v, eps := uvGCD(f[i], g[i], var, presenttheta):
        comptime_table[i] := time() - st:
        deg := degree(ddd, var);
        lastdeg := deg;
        bisectional := false;
        printf("theta = %f, Degree of GCD = %d\n", presenttheta, deg);
        # DEBUG();
        while deg <> givendeg do
            if deg < givendeg then
                if bisectional then
                    presenttheta := (thetamax + lasttheta) / 2.0;
                    thetamin := lasttheta;
                else
                    if lastdeg > givendeg then
                        bisectional := true;
                        presenttheta := (thetamax + thetamin) / 2.0;
                        thetamin := lasttheta;
                    else
                        presenttheta := 10 * lasttheta;
                        thetamin := lasttheta;
                        thetamax := presenttheta;
                    end;
                end;
            else
                if bisectional then
                    presenttheta := (thetamin + lasttheta) / 2.0;
                    thetamax := lasttheta;
                else
                    if lastdeg < givendeg then
                        bisectional := true;
                        presenttheta := (thetamax + thetamin) / 2.0;
                        thetamax := lasttheta;
                    else
                        presenttheta := 0.1 * lasttheta;
                        thetamax := lasttheta;
                        thetamin := presenttheta;
                    end;
                end;
            end;
            st := time():
            ddd, u, v, eps := uvGCD(f[i], g[i], var, presenttheta):
            comptime_table[i] := time() - st:
            lasttheta := presenttheta;
            lastdeg := deg;
            deg := degree(ddd, var);
            printf("theta = %f, Degree of GCD = %d\n", presenttheta, deg);
        end do;
        ppp := expand(ddd * u):
        qqq := expand(ddd * v):
        perturb2 := norm(expand(ppp - f[i]), 2) ^ 2 +
        norm(expand(qqq - g[i]), 2) ^ 2;
        perturb2_table[i] := perturb2:
        perturb_table[i] := sqrt(perturb2):
        gcd_table[i] := ddd:
        pol_table[i] := [ppp, qqq]:
        deg_table[i] := deg:
        theta_table[i] := presenttheta:
    end;
    return comptime_table, perturb_table, perturb2_table, gcd_table,
    pol_table, deg_table, theta_table
end proc;
