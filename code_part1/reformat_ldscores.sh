for pop in EUR AFR AMR CSA EAS MID; do
  echo $pop
  awk 'BEGIN {
    FS="\t";
    OFS="\t"
  }
  NR==1 {
    print "CHR", "BP", "A0", "A1", $3, $4, $5, $6;
    next
  }
  {
    split($1,a,":");
    CHR=a[1];
    BP=a[2];
    split($2,b,",");
    gsub(/[^[:alpha:]]/, "", b[1]);
    gsub(/[^[:alpha:]]/, "", b[2]);
    print CHR, BP, b[1], b[2], $3, $4, $5, $6
  }' UKBB.${pop}.ldscores_FULL.txt > UKBB.${pop}.ldscores_FULL_adj.txt
done