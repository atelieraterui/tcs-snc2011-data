test_gpgcd := proc(f, g, var, deg, stopcriterion, numite_bound,
                   case_begin, case_end)
    # test_gpcd used with gpgcd-real.mpl
    # in gpnewtonsolve, df includes the differences of polynomials A and B

    # f: array of polynomials
    # g: array of polynomials
    # var: main variable of polynomials
    # deg: degree of common divisor
    # stopcriterion: stopcriterion for gpcd
    # numite_bound: upper bound of the number of iterations in gpcd
    # case_begin: index of the beginning of the test cases
    # case_end: index of the end of the test cases
    #
    local i, j, numofcases, st,
    perturb, perturb2, ddd, ppp, qqq, anew, bnew, vv, numite,
    comptime_table, perturb_table, perturb2_table, cd_table, pol_table,
    cofactor_table, numite_table;

    comptime_table := Array(case_begin .. case_end);
    perturb_table := Array(case_begin .. case_end);
    perturb2_table := Array(case_begin .. case_end);
    cd_table := Array(case_begin .. case_end);
    pol_table := Array(case_begin .. case_end);
    cofactor_table := Array(case_begin .. case_end);
    numite_table := Array(case_begin .. case_end);
    printf("Degree of common divisor: %d\n", deg);
    for i from case_begin to case_end do
        printf("Case %d:\n", i);
        st := time():
        perturb, perturb2, ddd, ppp, qqq, anew, bnew, numite :=
        gpgcd_real(f[i], g[i], var, deg, stopcriterion, numite_bound):
        comptime_table[i] := time() - st:
        perturb_table[i] := perturb:
        perturb2_table[i] := perturb2:
        cd_table[i] := ddd:
        pol_table[i] := [ppp, qqq]:
        cofactor_table[i] := [anew, bnew]:
        numite_table[i] := numite:
    end;
    return comptime_table, perturb_table, perturb2_table, cd_table, pol_table,
    cofactor_table, numite_table
end proc;

test_gpgcd_complex :=
proc(f, g, var, deg, stopcriterion, numite_bound, case_begin, case_end)
    # fixdegtest_gpcd used with gpgcd-complex-01.mpl
    # in gpnewtonsolve, df includes the differences of polynomials A and B

    # f: array of polynomials
    # g: array of polynomials
    # var: main variable of polynomials
    # deg: degree of common divisor
    # stopcriterion: stopcriterion for gpcd
    # numite_bound: upper bound of the number of iterations in gpcd
    # case_begin: index of the beginning of the test cases
    # case_end: index of the end of the test cases
    #
    local i, j, numofcases, st,
    perturb, perturb2, ddd, ppp, qqq, anew, bnew, vv, numite,
    comptime_table, perturb_table, perturb2_table, cd_table, pol_table,
    cofactor_table, numite_table;

    comptime_table := Array(case_begin .. case_end);
    perturb_table := Array(case_begin .. case_end);
    perturb2_table := Array(case_begin .. case_end);
    cd_table := Array(case_begin .. case_end);
    pol_table := Array(case_begin .. case_end);
    cofactor_table := Array(case_begin .. case_end);
    numite_table := Array(case_begin .. case_end);
    printf("Degree of common divisor: %d\n", deg);
    for i from case_begin to case_end do
        printf("Case %d:\n", i);
        st := time():
        perturb, perturb2, ddd, ppp, qqq, anew, bnew, numite :=
        gpgcd_complex(f[i], g[i], var, deg, stopcriterion, numite_bound):
        comptime_table[i] := time() - st:
        perturb_table[i] := perturb:
        perturb2_table[i] := perturb2:
        cd_table[i] := ddd:
        pol_table[i] := [ppp, qqq]:
        cofactor_table[i] := [anew, bnew]:
        numite_table[i] := numite:
    end;
    return comptime_table, perturb_table, perturb2_table, cd_table, pol_table,
    cofactor_table, numite_table
end proc;
