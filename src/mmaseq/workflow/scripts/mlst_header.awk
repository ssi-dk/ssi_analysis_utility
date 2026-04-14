#!/usr/bin/awk -f
BEGIN{OFS="\t"}
NR==1{
  n=NF
  printf "file\tscheme\tst"
  for(i=4;i<=n;i++) printf "\tallele%d", i-3
  printf "\n"
}
{ print }