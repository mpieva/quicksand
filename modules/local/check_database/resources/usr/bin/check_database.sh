#!/usr/bin/env bash
#
# Check if the genome- and bedfile-names created by quicksand-build match
# the names that quicksand expects (derived from the taxonomy).
#
# quicksand expects (parse_taxonomy.py + 05_mapbam.nf + 07_bedfilterbam.nf):
#   genome : <genomes>/<Family>/<ScientificName>.fasta
#   bedfile: <bedfiles>/<ScientificName>.masked.bed
#   with ScientificName = names.dmp 'scientific name' of the taxid in seqid2taxid.map,
#   spaces replaced by '_' (nothing else)
#   and Family = scientific name of the nearest ancestor with rank 'family' (as in the kraken report)
#
# quicksand-build names the files (main.nf + parse_taxonomy.py):
#   RefSeq genomes: the highest 'subspecies' in the lineage, else 'species', ' ' and '/' replaced by '_'
#   custom genomes: the accession (= scientific name of the fake taxid)
#
# Usage:
#   check_database.sh -d DB_DIR -g GENOMES_DIR -b BEDFILES_DIR > report.tsv
#
# Output (TSV, one line per taxid in seqid2taxid.map):
#   TaxID Rank Family QuicksandName BuildName BuildCollisions Genome Bedfile Status
# Status:
#   OK             genome and bedfile found
#   MISSING_GENOME no genome file with the expected name -> silently dropped in quicksand
#   MISSING_BED    genome found, but no bedfile -> bedfiltering skipped (with warning)
#   NO_FAMILY      no family in the lineage -> never reaches the mapping step in quicksand
#   NOT_IN_TAXONOMY taxid missing in nodes.dmp/names.dmp
# BuildCollisions: number of taxids that quicksand-build writes to the same file (>1 = files overwrite each other)
#
# Afterwards, genome files that no taxid points to are listed with Status UNUSED_GENOME
# A summary is printed to stderr.

set -euo pipefail

usage() { sed -n '2,/^$/p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

db="" genomes="" bedfiles=""
while getopts "d:g:b:h" opt; do
    case $opt in
        d) db=$OPTARG ;;
        g) genomes=$OPTARG ;;
        b) bedfiles=$OPTARG ;;
        *) usage ;;
    esac
done
[[ -n $db && -n $genomes && -n $bedfiles ]] || usage

for f in "$db/taxonomy/nodes.dmp" "$db/taxonomy/names.dmp" "$db/seqid2taxid.map"; do
    [[ -f $f ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
[[ -d $genomes  ]] || { echo "ERROR: missing directory $genomes"  >&2; exit 1; }
[[ -d $bedfiles ]] || { echo "ERROR: missing directory $bedfiles" >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# genome files as quicksand sees them: ${genomes}/*/*.fasta --> "Family<TAB>Name"
find -L "$genomes" -mindepth 2 -maxdepth 2 -type f -name '*.fasta' \
    | awk -F/ '{ fam=$(NF-1); name=$NF; sub(/\.fasta$/, "", name); print fam "\t" name }' \
    | LC_ALL=C sort > "$tmp/genomes.tsv"

# bedfiles as quicksand sees them: ${bedfiles}/* matching (.+)\.masked\.bed --> "Name"
find -L "$bedfiles" -mindepth 1 -maxdepth 1 -type f -name '*.masked.bed*' \
    | awk -F/ '{ name=$NF; i=0; s=name
                 while ((p=index(s, ".masked.bed")) > 0) { i+=p; s=substr(s, p+1) }   # last occurrence (greedy regex)
                 print substr(name, 1, i-1) }' \
    > "$tmp/beds.txt"

# the taxids with genomes (unique)
cut -f2 "$db/seqid2taxid.map" | tr -d ' \r' | sort -u > "$tmp/taxids.txt"

awk -v GEN="$tmp/genomes.tsv" -v BED="$tmp/beds.txt" -v TAX="$tmp/taxids.txt" '
BEGIN {
    OFS = "\t"
    while ((getline line < GEN) > 0) { split(line, a, "\t"); genome[a[1] "\t" a[2]] = 1; gfiles[++ng] = a[1] "\t" a[2] }
    while ((getline line < BED) > 0) { bed[line] = 1 }
    while ((getline line < TAX) > 0) { if (line != "") taxids[++n] = line }
}
# nodes.dmp: taxid | parent | rank | ... | comment
FILENAME ~ /nodes\.dmp$/ {
    split($0, f, /\t\|\t/)
    t = f[1] + 0; parent[t] = f[2] + 0; rank[t] = f[3]
    if ($0 ~ /created by the quicksand-build pipeline/) custom[t] = 1
    next
}
# names.dmp: same substring check as quicksand parse_taxonomy.py
FILENAME ~ /names\.dmp$/ {
    if (index($0, "scientific name") == 0) next
    split($0, f, /\t\|\t/)
    sci[f[1] + 0] = f[2]
    next
}
function under(s) { gsub(/ /, "_", s); return s }
END {
    print "TaxID", "Rank", "Family", "QuicksandName", "BuildName", "BuildCollisions", "Genome", "Bedfile", "Status"

    for (i = 1; i <= n; i++) {
        t = taxids[i] + 0
        if (!(t in parent) || !(t in sci)) {
            row[i] = t OFS "-" OFS "-" OFS "-" OFS "-" OFS "BUILDNAME_NA" OFS "-" OFS "-" OFS "NOT_IN_TAXONOMY"
            continue
        }
        # walk up the lineage
        fam = ""; sp = ""; subsp = ""; cur = t; steps = 0
        while (steps++ < 200) {
            r = rank[cur]
            if (r == "family"     && fam == "") fam = sci[cur]
            if (r == "species")    sp = sci[cur]      # highest wins, like the dict in quicksand-build
            if (r == "subspecies") subsp = sci[cur]
            if (!(cur in parent) || parent[cur] == cur) break
            cur = parent[cur]
        }
        qname = under(sci[t])
        if (t in custom)      bname = sci[t]
        else if (subsp != "") bname = subsp
        else if (sp != "")    bname = sp
        else                  bname = ""
        gsub(/[ \/]/, "_", bname)
        if (bname != "") collisions[bname]++

        g = (fam != "" && ((fam "\t" qname) in genome)) ? "yes" : "no"
        b = (qname in bed) ? "yes" : "no"
        if (fam != "") used[fam "\t" qname] = 1

        if (fam == "")      status = "NO_FAMILY"
        else if (g == "no") status = "MISSING_GENOME"
        else if (b == "no") status = "MISSING_BED"
        else                status = "OK"

        bn[i] = bname
        row[i] = t OFS rank[t] OFS (fam == "" ? "-" : fam) OFS qname OFS (bname == "" ? "-" : bname)
        rest[i] = g OFS b OFS status
        count[status]++
    }
    for (i = 1; i <= n; i++) {
        if (i in rest) print row[i], (bn[i] == "" ? 0 : collisions[bn[i]]), rest[i]
        else { print row[i]; count["NOT_IN_TAXONOMY"]++ }
    }
    # genome files that no taxid points to
    for (j = 1; j <= ng; j++) if (!(gfiles[j] in used)) {
        split(gfiles[j], a, "\t")
        print "-", "-", a[1], a[2], "-", "-", "yes", ((a[2] in bed) ? "yes" : "no"), "UNUSED_GENOME"
        count["UNUSED_GENOME"]++
    }
    ncol = 0
    for (k in collisions) if (collisions[k] > 1) ncol++

    printf "\nSummary (%d taxids in seqid2taxid.map)\n", n > "/dev/stderr"
    split("OK MISSING_GENOME MISSING_BED NO_FAMILY NOT_IN_TAXONOMY UNUSED_GENOME", order, " ")
    for (j = 1; j <= 6; j++) printf "  %-16s %d\n", order[j], count[order[j]] + 0 > "/dev/stderr"
    printf "  %-16s %d (file names that quicksand-build writes more than once)\n", "BUILD_COLLISION", ncol > "/dev/stderr"
}
' "$db/taxonomy/nodes.dmp" "$db/taxonomy/names.dmp"
