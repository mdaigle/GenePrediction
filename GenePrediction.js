var fs = require('fs');
const assert = require('assert');

var STOP_CODONS = ["TAA", "TAG", "TGA"];
var NUCLEOTIDES = ["A", "C", "G", "T"];

var args = process.argv.slice(2);
var fasta_file = args[0];
var gbb_file = args[1];
var fasta_contents = fs.readFileSync(fasta_file, 'utf8');
var gbb_contents = fs.readFileSync(gbb_file, 'utf8');

// Parse fasta file
fasta_contents = fasta_contents.slice(fasta_contents.indexOf("\n") + 1); // Remove comment line
fasta_contents = fasta_contents.split(">", 1)[0]; // Remove plasmids
fasta_contents = fasta_contents.replace(/\s/g, ""); // Remove whitespace

// Parse Genbank file
gbb_contents = gbb_contents.split("ORIGIN")[0];
gbb_contents = gbb_contents.split(/\n/g); // Break into lines
var sequences = [];
gbb_contents.forEach(function(line){
    var regex = /\s+CDS\s+([0-9]+)\.\.([0-9]+)/g;
    var sequence = regex.exec(line);
    if (sequence && sequence.length == 3) {
        sequences.push({
            "start": sequence[1],
            "end": sequence[2]
        });
    }
});
console.log("CDSs Found: " + sequences.length);


var num_short = 0;
var num_long = 0;
var num_orf = 0;

var short_sequences = [];
var long_sequences = [];

// Start indices
for (var start_index = 0; start_index < 3; start_index++) {
    console.log("Aligned at " + start_index);

    var orf_count = 0;
    var gene_start = start_index;

    for (var index = start_index; index + 3 <= fasta_contents.length; index = index + 3) {
        var triplet = fasta_contents.substring(index, (index + 3));

        // Replace invalid nucleotides with T
        for (var i = 0; i < 3; i++) {
            if (NUCLEOTIDES.indexOf(triplet[i].toUpperCase()) == -1)
            {
                //console.log("changed nucleotide " + triplet[i]);
                triplet = triplet.replace(triplet.charAt(i), "T");
            }
        }

        // Continue if not a stop codon
        if (STOP_CODONS.indexOf(triplet) < 0)
            continue;

        //console.log((gene_start + 1) + " " + (index + 3));
        var frame = fasta_contents.substring(gene_start, index);
        orf_count++;

        if (((index + 3) - gene_start) < 50) {
            num_short++;
            short_sequences.push(frame);
        }

        if (((index + 3) - gene_start) > 1400) {
            num_long++;
            long_sequences.push(frame);
        }

        gene_start = index + 3;
    }
    console.log("ORF Count: " + orf_count);
    num_orf += orf_count;
}

console.log("Total ORFs: " + num_orf);
console.log("Total short: " + num_short);
console.log("Total long: " + num_long);
var p_probs = {};
for (var i = 0; i < 4; i++) {
    p_probs[i] = {};
    for (var j = 0; j < 4; j++){
        p_probs[i][j] = {};
        for (var k = 0; k < 4; k++) {
            p_probs[i][j][k] = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0};
        }
    }
}

// Analyze ORFs
/*var p_probs = {
    "A": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    },
    "C": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    },
    "T": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    },
    "G": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    }
}*/

var q_probs = {
    "A": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    },
    "C": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    },
    "T": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    },
    "G": {
        "A": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "C": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "T": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] },
        "G": { "A": [0, 0, 0, 0, 0], "C": [0, 0, 0, 0, 0], "G": [0, 0, 0, 0, 0], "T": [0, 0, 0, 0, 0] }
    }
}

function countSequences(sequence, index, arr, prob_json) {
    for (var i = 0; i + 3 < sequence.length; i++) {
        for (j = 0; j < 3; j++) {
            if (NUCLEOTIDES.indexOf(sequence[i + j].toUpperCase()) == -1) {
                sequence = sequence.replace(sequence.charAt(i + j), "T");
            }
        }

        var first = NUCLEOTIDES.indexOf(sequence.charAt(i));
        var second = NUCLEOTIDES.indexOf(sequence.charAt(i + 1));
        var third = NUCLEOTIDES.indexOf(sequence.charAt(i + 2));

        p_probs[first][second][third][0]++;

        if (i + 4 < sequence.length) {
            var conditional = NUCLEOTIDES.indexOf(sequence.charAt(i + 3)) + 1;
            p_probs[first][second][third][conditional]++;
        }
    }
}

long_sequences.forEach(countSequences, p_probs);
console.log(p_probs);
