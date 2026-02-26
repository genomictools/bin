BEGIN {
    OFS = "\t"

    if (ARGC != 3) {
        print "Usage: awk -f aggregate_cov.awk <cov_list.txt> <stats_list.txt>" > "/dev/stderr"
        exit 1
    }

    cov_list_file = ARGV[1]
    stats_list_file = ARGV[2]

    while ((getline line < cov_list_file) > 0) {
        if (line != "") cov_files[++cov_n] = line
    }
    close(cov_list_file)

    while ((getline line < stats_list_file) > 0) {
        if (line != "") stats_files[++stats_n] = line
    }
    close(stats_list_file)

    if (cov_n != stats_n) {
        print "ERROR: number of cov files and stats files must match" > "/dev/stderr"
        exit 1
    }

    for (i = 1; i <= cov_n; i++) {
        cov_file = cov_files[i]
        stats_file = stats_files[i]

        nf = normalization_factor(stats_file)
        aggregate_sample(cov_file, nf)
    }

    for (k in n) {
        mean = sum[k] / n[k]
        sd   = sqrt((sumsq[k] / n[k]) - (mean * mean))
        print k, mean, sd
    }

    exit
}

function extract_mapped_duplicates(stats_file, out,   line, x, mapped, dups) {
    mapped = -1
    dups   = -1

    while ((getline line < stats_file) > 0) {
        if (mapped < 0 && line ~ /^[0-9][0-9]* \+ [0-9][0-9]* mapped/) {
            x = line
            sub(/ .*/, "", x)
            mapped = x + 0
        }
        if (dups < 0 && line ~ /^[0-9][0-9]* \+ [0-9][0-9]* duplicates/) {
            x = line
            sub(/ .*/, "", x)
            dups = x + 0
        }
    }
    close(stats_file)

    out["mapped"] = mapped
    out["dups"]   = dups
}

function normalization_factor(stats_file, md, nf) {
    delete md
    extract_mapped_duplicates(stats_file, md)

    if (md["mapped"] < 0 || md["dups"] < 0)
        return 0

    nf = (md["mapped"] - md["dups"]) / 100000000
    return (nf > 0 ? nf : 0)
}

function aggregate_sample(cov_file, nf,   line, f, chr, pos, cov, key, v) {
    if (nf <= 0) return

    while ((getline line < cov_file) > 0) {
        if (line == "" || line ~ /^#/) continue
        split(line, f, /[ \t]+/)
        chr = f[1]
        pos = f[2] + 0
        cov = f[3] + 0

        key = chr OFS pos
        v = cov / nf

        n[key]++
        sum[key] += v
        sumsq[key] += (v * v)
    }
    close(cov_file)
}
