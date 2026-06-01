BEGIN { OFS = "\t" }

function ref_len(cigar,    rest, n, op, len) {
    rest = cigar
    len = 0
    while (match(rest, /^([0-9]+)([MIDNSHP=X])/, m)) {
        n = m[1] + 0
        op = m[2]
        if (op ~ /[MDN=X]/) {
            len += n
        }
        rest = substr(rest, RLENGTH + 1)
    }
    return len
}

$0 !~ /^@/ {
    flag = $2 + 0
    chrom = $3
    pos = $4 + 0
    cigar = $6
    if (chrom == "*" || cigar == "*") {
        next
    }

    if (and(flag, 16)) {
        tss = pos + ref_len(cigar) - 1
        strand = "-"
    } else {
        tss = pos
        strand = "+"
    }

    print chrom, tss - 1, tss, $1, 1, strand
}
