BEGIN { OFS="\t" }
{ gsub(/\r/, "") }

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

$0 == "[Data]" { in_data=1; next }

in_data && !hdr {
  for (i=1; i<=NF; i++) {
    if (matches($i, snp_name))    idx["SNP Name"]=i
    if (matches($i, chr))         idx["Chr"]=i
    if (matches($i, position))    idx["Position"]=i
    if (matches($i, allele_ab_a)) idx["Allele1 - AB"]=i
    if (matches($i, allele_ab_b)) idx["Allele2 - AB"]=i
    if (matches($i, allele_a))    idx["Allele - Top A"]=i
    if (matches($i, allele_b))    idx["Allele - Top B"]=i
  }
  print "SNP Name","Chr","Position", sample_id ".GType", sample_id ".Top Alleles"
  hdr=1
  next
}

in_data && hdr {
  id    = $idx["SNP Name"]
  chr   = $idx["Chr"]
  pos   = $idx["Position"]
  a_ab  = $idx["Allele1 - AB"]
  b_ab  = $idx["Allele2 - AB"]
  a_top = $idx["Allele - Top A"]
  b_top = $idx["Allele - Top B"]

  if (id=="" || chr=="" || pos=="") next
  sub(/^(chr|CHR)/, "", chr)

  gtype = "NC"
  top = "--"

  if (a_ab=="A" && b_ab=="A") gtype = "AA"
  else if ((a_ab=="A" && b_ab=="B") || (a_ab=="B" && b_ab=="A")) gtype = "AB"
  else if (a_ab=="B" && b_ab=="B") gtype = "BB"

  if (a_top ~ /^[ACGT]$/ && b_top ~ /^[ACGT]$/) top = a_top b_top

  print id, chr, pos, gtype, top
}
