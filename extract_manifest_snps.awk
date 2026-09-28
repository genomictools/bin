BEGIN { assay_section=0; name_col=0; chr_col=0 }

# Returns 1 if name matches any entry in a comma-separated candidates string
function matches(name, candidates,    n, arr, i) {
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", name)
  n = split(candidates, arr, ",")
  for (i = 1; i <= n; i++) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", arr[i])
    if (arr[i] == name) return 1
  }
  return 0
}

/^\[Assay\]/ {
  assay_section=1
  getline
  gsub(/\r/, "")
  split($0, cols, ",")
  for (i=1; i<=length(cols); i++) {
    if (matches(cols[i], snp_name)) {
      name_col=i
    }
    if (matches(cols[i], chr)) {
      chr_col=i
    }
  }
  next
}

/^\[Controls\]/ {
  assay_section=0
}

assay_section {
  gsub(/\r/, "")
  if (name_col > 0 && chr_col > 0 && $0 != "") {
    split($0, vals, ",")
    row_name = vals[name_col]
    row_chr  = vals[chr_col]
    chr_norm = tolower(row_chr)
    sub(/^chr/, "", chr_norm)

    if (row_name != "" && chr_norm ~ /^([1-9]|1[0-9]|2[0-2])$/) print row_name
  }
}
