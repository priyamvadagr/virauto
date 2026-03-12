awk -F',' '
NR==1{
  for(i=1;i<=NF;i++) if($i=="receptor_id") c=i
  if(!c){ print "receptor_id column not found"; exit 1 }
  next
}
{
  total++
  v=$c
  gsub(/^[ \t\r\n"]+|[ \t\r\n"]+$/, "", v)
  if(v!="") non_empty++
}
END{
  printf "total_rows=%d\nnon_empty_receptor_id_rows=%d\nall_receptor_id_empty=%s\n",
         total, non_empty, (non_empty==0?"True":"False")
}' /ix/djishnu/Priyamvada/virauto/data/epitopes/iedb/mhc_i/iedb_epitopes_with_tcr_mhc_i.csv
