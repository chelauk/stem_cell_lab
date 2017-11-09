#Specify the input file
XLS=$1
#Specify the gene ID column
ID=$2
#Specify the fold change value column
FC=$3
#Specify the raw p-value column
P=$4

sed 1d $XLS | tr -d '"' \
| awk -v I=$ID -v F=$FC -v P=$P '{FS="\t"} $I!="" {print $I, $F, $P}' \
| awk '$2!="NA" && $3!="NA"' \
| awk '{s=1} $2<0{s=-1} {print $1"\t"s*-1*log($3)/log(10)}' \
| awk '{print $1,$2,$2*$2}' | sort -k3gr \
| awk '{OFS="\t"} !arr[$1]++ {print $1,$2}' \
| sort -k2gr > ${XLS}.rnk
