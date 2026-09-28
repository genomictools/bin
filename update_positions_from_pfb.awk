BEGIN { FS=OFS="\t" }

NR==FNR {
  if (FNR == 1) next
  pfb_pos[$1] = $3
  next
}

FNR == 1 {
  print
  next
}

{
  if (($1 in pfb_pos) && pfb_pos[$1] != "") {
    $3 = pfb_pos[$1]
  }
  print
}