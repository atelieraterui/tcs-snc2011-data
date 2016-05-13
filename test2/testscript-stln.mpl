test_stln := proc(f, g, var, deg, stopcriterion, numite_bound,
                      case_begin, case_end)
    # test_stln used with con_mulpoly_gcd_2.0-01.mpl

    # f: array of polynomials
    # g: array of polynomials
    # var: main variable of polynomials (dummy argument)
    # deg: degree of common divisor
    # stopcriterion: m = stopcriterion for  with 1e-m
    # numite_bound: upper bound of the number of iterations in gpcd
    # (dummy augument)
    # case_begin: index of the beginning of the test cases
    # case_end: index of the end of the test cases
    #
    local i, j, numofcases, st, ddd, prim_list,
    perturb, perturb2, ppp, qqq, anew, bnew, vv, numite,
    comptime_table, perturb_table, perturb2_table, pol_table, numite_table;
    comptime_table := Array(case_begin .. case_end);
    perturb_table := Array(case_begin .. case_end);
    perturb2_table := Array(case_begin .. case_end);
    pol_table := Array(case_begin .. case_end);
    numite_table := Array(case_begin .. case_end);
    printf("Degree of common divisor: %d\n", deg);
    for i from case_begin to case_end do
        printf("Case %d:\n", i);
        st := time():
        ddd, prim_list, numite :=
        R_con_mulpoly([f[i], g[i]], deg, stopcriterion,
                      'default', 'default', 'cofa'):
        comptime_table[i] := time() - st:
        ppp := expand(ddd * prim_list[1]):
        qqq := expand(ddd * prim_list[2]):
        perturb2 := norm(expand(ppp - f[i]), 2) ^ 2 +
        norm(expand(qqq - g[i]), 2) ^ 2;
        perturb2_table[i] := perturb2:
        perturb_table[i] := sqrt(perturb2):
        pol_table[i] := [ppp, qqq]:
        numite_table[i] := numite:
    end;
    return comptime_table, perturb_table, perturb2_table, pol_table,
    numite_table
end proc;

test_stln_complex := proc(f, g, var, deg, stopcriterion, numite_bound,
                      case_begin, case_end)
    # fixdegtest_C_con_mulpoly used with con_mulpoly_gcd_2.0-02.mpl

    # f: array of polynomials
    # g: array of polynomials
    # var: main variable of polynomials (dummy argument)
    # deg: degree of common divisor
    # stopcriterion: m = stopcriterion for  with 1e-m
    # numite_bound: upper bound of the number of iterations in gpcd
    # (dummy augument)
    # case_begin: index of the beginning of the test cases
    # case_end: index of the end of the test cases
    #
    local i, j, numofcases, st, ddd, prim_list,
    perturb, perturb2, ppp, qqq, anew, bnew, vv, numite,
    comptime_table, perturb_table, perturb2_table, pol_table, numite_table;
    comptime_table := Array(case_begin .. case_end);
    perturb_table := Array(case_begin .. case_end);
    perturb2_table := Array(case_begin .. case_end);
    pol_table := Array(case_begin .. case_end);
    numite_table := Array(case_begin .. case_end);
    printf("Degree of common divisor: %d\n", deg);
    for i from case_begin to case_end do
        printf("Case %d:\n", i);
        st := time():
        ddd, prim_list, numite :=
        C_con_mulpoly([f[i], g[i]], deg, stopcriterion,
                      'default', 'default', 'cofa'):
        comptime_table[i] := time() - st:
        ppp := expand(ddd * prim_list[1]):
        qqq := expand(ddd * prim_list[2]):
        perturb2 := norm(expand(ppp - f[i]), 2) ^ 2 +
        norm(expand(qqq - g[i]), 2) ^ 2;
        perturb2_table[i] := perturb2:
        perturb_table[i] := sqrt(perturb2):
        pol_table[i] := [ppp, qqq]:
        numite_table[i] := numite:
    end;
    return comptime_table, perturb_table, perturb2_table, pol_table,
    numite_table
end proc;
