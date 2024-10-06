### Count total number of variants by effect
```
scratch
cd snpEff/snpEff_out


for f in $(ls *.g.genes.txt); do i=${f%.g.genes.txt}; high=$(cut -f5 $f | tail -n +3 | sum_stdin.py); low=$(cut -f6 $f | tail -n +3 | sum_stdin.py); mod=$(cut -f7 $f | tail -n +3 | sum_stdin.py); modif=$(cut -f8 $f | tail -n +3 | sum_stdin.py); echo -e "$i\t$high\t$low\t$mod\t$modif"; done > /nese/meclab/Katya/snpEff_turtles/all_samples.total_variant_by_effect.tsv

```
