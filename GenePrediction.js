var fs = require('fs');
var math = require('mathjs');

math.config({
    number: 'BigNumber',
    precision: 100
});

const assert = require('assert');

var STOP_CODONS = ["TAA", "TAG", "TGA"];
var NUCLEOTIDES = ["A", "C", "G", "T"];

var args = process.argv.slice(2);
var fasta_file = args[0];
var gbb_file = args[1];
var fasta_contents = fs.readFileSync(fasta_file, 'utf8');
var gbb_contents = fs.readFileSync(gbb_file, 'utf8');

// Get interesting stuff from fasta file
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
var all_sequences = [];

// Find ORFs
for (var start_index = 0; start_index < 3; start_index++) {
    var orf_count = 0;
    var gene_start = start_index;

    for (var index = start_index; index + 3 <= fasta_contents.length; index = index + 3) {
        var triplet = fasta_contents.substring(index, (index + 3));

        // Replace invalid nucleotides with T
        for (var i = 0; i < 3; i++) {
            if (NUCLEOTIDES.indexOf(triplet[i].toUpperCase()) == -1)
            {
                var bad_char = triplet.charAt(i);
                triplet = triplet.replace(bad_char, "T");
                fasta_contents = fasta_contents.replace(bad_char, "T");
            }
        }

        // Continue if not a stop codon
        if (STOP_CODONS.indexOf(triplet) < 0)
            continue;

        // Get the whole ORF
        var frame = fasta_contents.substring(gene_start, index + 3);
        orf_count++;

        // Short ORF threshold
        if (frame.length < 50) {
            num_short++;
            short_sequences[gene_start] = {
                "frame": frame,
                "start": gene_start + 1,
                "length": frame.length
            };
        }

        // Long ORF threshold
        if (frame.length > 1400) {
            num_long++;
            long_sequences[gene_start] = {
                "frame": frame,
                "start": gene_start + 1,
                "length": frame.length
            };
        }

        // Update gene start index
        gene_start = index + 3;
    }

    console.log("Alignment " + start_index + " ORF Count: " + orf_count);
    num_orf += orf_count;
}

console.log("Total ORFs: " + num_orf);
console.log("Total short: " + num_short);
console.log("Total long: " + num_long);

// Set up count and probability records
var p_counts = {};
var q_counts = {};

var p_probs = {};
var q_probs = {};

for (var i = 0; i < 4; i++) {
    p_counts[i] = {};
    q_counts[i] = {};

    p_probs[i] = {};
    q_probs[i] = {};
    for (var j = 0; j < 4; j++){
        p_counts[i][j] = {};
        q_counts[i][j] = {};

        p_probs[i][j] = {};
        q_probs[i][j] = {};
        for (var k = 0; k < 4; k++) {
            p_counts[i][j][k] = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0};
            q_counts[i][j][k] = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0};

            p_probs[i][j][k] = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0};
            q_probs[i][j][k] = {"0": 0, "1": 0, "2": 0, "3": 0, "4": 0};
        }
    }
}

long_sequences.forEach(countSequences, p_counts);
short_sequences.forEach(countSequences, q_counts);

calculateProbs(p_counts, p_probs);
calculateProbs(q_counts, q_probs);

console.log("Prob of A given AAA = " + p_probs[0][0][0][1]);

for (var j = 0; j < 4; j++) {
    for (var k = 0; k < 4; k++) {
        console.log(p_counts[0][j][k][1] + ", " + p_counts[0][j][k][0]);
    }
}

console.log("Now Q");

for (var j = 0; j < 4; j++) {
    for (var k = 0; k < 4; k++) {
        console.log(q_counts[0][j][k][1] + ", " + q_counts[0][j][k][0]);
    }
}

long_sequences.forEach(calculateGeneProbability, p_probs);

//console.log("done with long seq probs");

short_sequences.forEach(calculateGeneProbability, q_probs);

//console.log("done calculating gene probs");

printFirstFive(long_sequences);
printFirstFive(short_sequences);

function countSequences(sequence_record, index, arr, prob_array) {
    if (sequence_record == undefined) {
        return;
    }
    var sequence = sequence_record.frame;
    for (var i = 0; i + 3 < sequence.length - 3; i++) {
        for (j = 0; j < 3; j++) {
            if (NUCLEOTIDES.indexOf(sequence[i + j].toUpperCase()) == -1) {
                sequence = sequence.replace(sequence.charAt(i + j), "T");
            }
            assert(NUCLEOTIDES.indexOf(sequence[i + j].toUpperCase()) >= 0);
        }

        var first = NUCLEOTIDES.indexOf(sequence.charAt(i));
        var second = NUCLEOTIDES.indexOf(sequence.charAt(i + 1));
        var third = NUCLEOTIDES.indexOf(sequence.charAt(i + 2));

        this[first][second][third][0]++;

        if (i + 3 < sequence.length - 3) {
            var conditional = NUCLEOTIDES.indexOf(sequence.charAt(i + 3)) + 1;
            this[first][second][third][conditional]++;
        }
    }
    arr[index].frame = sequence;
//}
}

function calculateProbs(count_arr, prob_arr) {
    var num_codons = 0;
    for (var i = 0; i < 4; i++) {
        for (var j = 0; j < 4; j++) {
            for (var k = 0; k < 4; k++) {
                num_codons += count_arr[i][j][k][0];
            }
        }
    }

    for (var i = 0; i < 4; i++) {
        for (var j = 0; j < 4; j++) {
            for (var k = 0; k < 4; k++) {
                prob_arr[i][j][k][0] = math.divide(math.bignumber(count_arr[i][j][k][0]), math.bignumber(num_codons));
                var total_count = math.bignumber(count_arr[i][j][k][0]);
                for (var l = 1; l < 5; l++) {
                    var count = math.bignumber(count_arr[i][j][k][l]);
                    prob_arr[i][j][k][l] = math.divide(count, total_count);
                }
                //assert(sum == total_count, "Sum and total don't match at " + i + " " + j + " " + k + " : " + sum + " " + total_count);
            }
        }
    }
}

function calculateGeneProbability(sequence_record, index, arr) {
    if (sequence_record == undefined) {
        return;
    }
    var sequence = sequence_record.frame;
    //console.log(sequence.length);
    //console.log("probs");

    // Look at first codon
    var p_prob = math.bignumber(1);
    var q_prob = math.bignumber(1);

    if (sequence_record['length'] > 3) {
        var first = NUCLEOTIDES.indexOf(sequence.charAt(0));
        var second = NUCLEOTIDES.indexOf(sequence.charAt(1));
        var third = NUCLEOTIDES.indexOf(sequence.charAt(2));

        assert(p_probs[first][second][third][0] != 0);
        assert(q_probs[first][second][third][0] != 0);

        p_prob = math.multiply(p_prob, p_probs[first][second][third][0]);
        q_prob = math.multiply(q_prob, q_probs[first][second][third][0]);

        //console.log(p_prob);
    }


    for (var i = 3; i < sequence_record['length'] - 3; i++) {
        var first = NUCLEOTIDES.indexOf(sequence.charAt(i - 3));
        var second = NUCLEOTIDES.indexOf(sequence.charAt(i - 2));
        var third = NUCLEOTIDES.indexOf(sequence.charAt(i - 1));
        // Add one to skip past the 3-mer prob
        var curr = NUCLEOTIDES.indexOf(sequence.charAt(i)) + 1;

        assert(curr > 0 && curr <= 4, "char is " + sequence.charAt(i));
        assert(p_probs[first][second][third][curr]);
        assert(q_probs[first][second][third][curr]);

        p_prob = math.multiply(p_prob, p_probs[first][second][third][curr]);
        q_prob = math.multiply(q_prob, q_probs[first][second][third][curr]);
    }

    var div = math.divide(math.bignumber(p_prob), math.bignumber(q_prob));
    var markov_score = math.log(div, math.bignumber(2));

    arr[index].score = markov_score;
}

function printFirstFive(seq_arr) {
    var first_five_sequences = [];
    var index = 0;
    for (var i = 0; i < seq_arr.length && index < 5; i++) {
        if (seq_arr[i] == undefined) {
            continue;
        }
        first_five_sequences[index++] = seq_arr[i];
    }

    for (var i = 0; i < 5; i++) {
        var sequence = first_five_sequences[i];
        console.log(i + " " + sequence['start'] + " " + sequence['length'] + " " + sequence['score']);
    }
}
