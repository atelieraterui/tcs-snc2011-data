arrayLargeElementIndex := proc (numArray, criterion)

    # Test if there exist number(s) > criterion in numArray. If such
    # element(s) exitsts, outout index(es) of such elements

    # Inputs:
    # numSet (Array):  set of numbers
    # criterion (numeric)

    # Output:
    # indexSet: set of indices {i} satisfying numSet[i] >= criterion

    local i, indexSet;

    indexSet := {};
    if type (numArray, Array) then
        for i in [$ op (2, eval (numArray))] do
            if numArray[i] > criterion then
                indexSet := indexSet union {i};
            end;
        end;
    end;
    if type (numArray, set) or type (numArray, list) then
        for i from 1 to nops (numArray) do
            if numArray[i] > criterion then
                indexSet := indexSet union {i};
            end;
        end;
    end;

    return indexSet;
end proc:
