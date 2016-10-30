var fs = require('fs');

var STOP_CODONS = ["TAA", "TAG", "TGA"];
var NUCLEOTIDES = ["A", "C", "G", "T"];

var args = process.argv.slice(2);
var fasta_file = args[0];
var fasta_contents = fs.readFileSync(fasta_file, 'utf8');

// Remove comment line
fasta_contents = fasta_contents.slice(fasta_contents.indexOf("\n") + 1);

// Remove plasmids
fasta_contents = fasta_contents.split(">", 1)[0];

// Remove whitespace
fasta_contents = fasta_contents.replace(/\s/g, "");

// Start indices
for (var start_index = 0; start_index < 3; start_index++) {
    console.log("Aligned at " + start_index);
    var gene_start = start_index;
console.log(gene_start + 1);
    for (var index = start_index; index + 3 <= fasta_contents.length; index = index + 3) {
        var triplet = fasta_contents.substring(index, (index + 3)).toUpperCase();

        // Replace invalid nucleotides with T
        for (var i = 0; i < 3; i++) {
            if (NUCLEOTIDES.indexOf(triplet[i]) == -1)
            {
                //console.log("changed nucleotide " + triplet[i]);
                triplet[i] = "T";
            }
        }

        // Check if the triplet is a stop codon
        if (STOP_CODONS.indexOf(triplet) >= 0)
        {
            if (((index + 3) - gene_start) > 2500) {
                console.log((gene_start + 1) + " " + (index + 3));
            }
            gene_start = index + 3;
        }
    }
}
