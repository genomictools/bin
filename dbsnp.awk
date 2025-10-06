#!/usr/bin/awk -f

BEGIN {
    OFS = "\t";
    print "SNP", "Chr", "Position", "PFB";
}

{
    count = 0;
    split($6, arr, ",");
    n = length(arr);
    for (i = 1; i <= n; i++)
        for (j = i + 1; j <= n; j++)
            if (arr[i] < arr[j]) {
                temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
    for (i = 1; i <= n; i++) {
        if (arr[i] != ".") {
            if (++count == 2) {
                print $1, $2, $3, arr[i];
                break;
            }
        }
    }
}