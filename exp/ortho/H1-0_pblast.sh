curl --form 'email=tgw325@alumni.ku.dk' \
     --form 'program=blastp' \
     --form 'matrix=BLOSUM62' \
     --form 'alignments=250' \
     --form 'scores=250' \
     --form 'exp=10' \
     --form 'filter=F' \
     --form 'gapalign=true' \
     --form 'compstats=F' \
     --form 'align=0' \
     --form 'stype=protein' \
     --form 'sequence=>sp|P07305|H10_HUMAN Histone H1.0 OS=Homo sapiens OX=9606 GN=H1-0 PE=1 SV=3
MTENSTSAPAAKPKRAKASKKSTDHPKYSDMIVAAIQAEKNRAGSSRQSIQKYIKSHYKV
GENADSQIKLSIKRLVTTGVLKQTKGVGASGSFRLAKSDEPKKSVAFKKTKKEIKKVATP
KKASKPKKAASKAPTKKPKATPVKKAKKKLAATPKKAKKPKTVKAKPVKASKPKKAKPVK
PKAKSSAKRAGKKK
' \
     --form 'database=uniprotkb_refprotswissprot' \
     https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run