#!/usr/bin/awk -f

BEGIN {
    OFS = "\t";
    print "SNP", "Chr", "Position", "PFB";
    
    # Check if we have the required arguments
    if (ARGC != 3) {
        print "Usage: script snplist_file pfb_file > output" > "/dev/stderr";
        exit 1;
    }
    
    snplist_file = ARGV[1];
    pfb_file = ARGV[2];
    
    # Read snplist file and extract rsids - store ALL matching names, not just one
    while ((getline line < snplist_file) > 0) {
        # Extract rsid from the snp name
        if (match(line, /rs[0-9]+/)) {
            rsid = substr(line, RSTART, RLENGTH);
            # Store multiple SNP names for the same rsid
            if (rsid in snp_names) {
                snp_names[rsid] = snp_names[rsid] "," line;
            } else {
                snp_names[rsid] = line;
            }
        }
    }
    close(snplist_file);
    
    # Set ARGV to only process the pfb_file
    ARGV[1] = pfb_file;
    ARGC = 2;
}

{
    # Extract rsid from PFB file SNP column ($1)
    if (match($1, /rs[0-9]+/)) {
        pfb_rsid = substr($1, RSTART, RLENGTH);
        
        # Check if this rsid exists in our snplist
        if (pfb_rsid in snp_names) {
            # Perform the aggregation (sort and get 2nd highest non-missing value)
            count = 0;
            split($6, arr, ",");
            n = length(arr);
            
            # Sort array in descending order
            for (i = 1; i <= n; i++)
                for (j = i + 1; j <= n; j++)
                    if (arr[i] < arr[j]) {
                        temp = arr[i];
                        arr[i] = arr[j];
                        arr[j] = temp;
                    }
            
            # Find the 2nd highest non-missing value
            pfb_value = "";
            for (i = 1; i <= n; i++) {
                if (arr[i] != ".") {
                    if (++count == 2) {
                        pfb_value = arr[i];
                        break;
                    }
                }
            }
            
            # Output one line for each SNP name that matches this rsid
            split(snp_names[pfb_rsid], names, ",");
            for (j in names) {
                print names[j], $2, $3, pfb_value;
            }
        }
    }
}