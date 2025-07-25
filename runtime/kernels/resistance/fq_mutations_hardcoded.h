#ifndef FQ_MUTATIONS_HARDCODED_H
#define FQ_MUTATIONS_HARDCODED_H

#include <vector>
#include <string>

// Forward declaration if not already included
struct FQResistanceMutation;

// Hardcoded FQ resistance mutations
// Generated from: /home/david/Documents/Code/biogpu/data/quinolone_resistance_mutation_table.csv
// Total mutations: 255
static const std::vector<FQResistanceMutation> HARDCODED_FQ_MUTATIONS = {
    {"Acinetobacter_baumannii", "gyrA", 79, 'S', 'L'},
    {"Acinetobacter_baumannii", "gyrA", 81, 'S', 'L'},
    {"Acinetobacter_baumannii", "parC", 84, 'E', 'K'},
    {"Acinetobacter_baumannii", "parC", 84, 'S', 'W'},
    {"Acinetobacter_baumannii", "parC", 84, 'S', 'F'},
    {"Acinetobacter_baumannii", "parC", 88, 'G', 'C'},
    {"Burkholderia_cenocepacia", "gyrA", 83, 'T', 'I'},
    {"Burkholderia_dolosa", "gyrA", 87, 'D', 'H'},
    {"Campylobacter_coli", "gyrA", 86, 'T', 'I'},
    {"Campylobacter_coli", "gyrA", 90, 'D', 'Y'},
    {"Campylobacter_jejuni", "gyrA", 86, 'P', 'S'},
    {"Campylobacter_jejuni", "gyrA", 90, 'D', 'N'},
    {"Campylobacter_jejuni", "gyrA", 104, 'T', 'A'},
    {"Campylobacter_lari", "gyrA", 86, 'T', 'V'},
    {"Campylobacter_upsaliensis", "gyrA", 86, 'T', 'K'},
    {"Citrobacter_farmeri", "gyrA", 87, 'D', 'G'},
    {"Citrobacter_freundii", "gyrA", 83, 'T', 'I'},
    {"Clostridioides_difficile", "gyrA", 71, 'T', 'I'},
    {"Clostridioides_difficile", "gyrA", 71, 'S', 'A'},
    {"Clostridioides_difficile", "gyrA", 81, 'I', 'R'},
    {"Clostridioides_difficile", "gyrA", 82, 'D', 'N'},
    {"Clostridioides_difficile", "gyrA", 82, 'S', 'V'},
    {"Clostridioides_difficile", "gyrA", 82, 'S', 'A'},
    {"Clostridioides_difficile", "gyrA", 118, 'D', 'V'},
    {"Clostridioides_difficile", "gyrA", 118, 'A', 'S'},
    {"Clostridioides_difficile", "gyrB", 130, 'V', 'I'},
    {"Clostridioides_difficile", "gyrB", 139, 'E', 'K'},
    {"Clostridioides_difficile", "gyrB", 366, 'T', 'A'},
    {"Clostridioides_difficile", "gyrB", 366, 'T', 'V'},
    {"Clostridioides_difficile", "gyrB", 416, 'D', 'G'},
    {"Clostridioides_difficile", "gyrB", 426, 'R', 'K'},
    {"Clostridioides_difficile", "gyrB", 426, 'D', 'V'},
    {"Clostridioides_difficile", "gyrB", 447, 'A', 'T'},
    {"Clostridioides_difficile", "gyrB", 466, 'D', 'N'},
    {"Clostridioides_difficile", "gyrB", 466, 'E', 'V'},
    {"Enterobacter_cloacae", "gyrA", 83, 'S', 'Y'},
    {"Enterobacter_cloacae", "gyrA", 83, 'S', 'I'},
    {"Enterobacter_cloacae", "gyrA", 83, 'D', 'N'},
    {"Enterobacter_cloacae", "gyrA", 87, 'S', 'I'},
    {"Enterobacter_cloacae", "gyrA", 87, 'D', 'G'},
    {"Enterobacter_cloacae", "gyrA", 87, 'S', 'F'},
    {"Enterobacter_cloacae", "parC", 80, 'D', 'H'},
    {"Enterococcus_faecalis", "gyrA", 83, 'S', 'I'},
    {"Enterococcus_faecalis", "gyrA", 83, 'S', 'R'},
    {"Enterococcus_faecalis", "gyrA", 83, 'E', 'K'},
    {"Enterococcus_faecalis", "gyrA", 83, 'S', 'R'},
    {"Enterococcus_faecalis", "gyrA", 87, 'E', 'K'},
    {"Enterococcus_faecalis", "gyrA", 87, 'E', 'G'},
    {"Enterococcus_faecalis", "parC", 80, 'S', 'I'},
    {"Enterococcus_faecalis", "parC", 80, 'E', 'A'},
    {"Enterococcus_faecalis", "parC", 84, 'S', 'Y'},
    {"Enterococcus_faecalis", "parC", 84, 'S', 'N'},
    {"Enterococcus_faecium", "gyrA", 83, 'S', 'I'},
    {"Enterococcus_faecium", "gyrA", 83, 'S', 'R'},
    {"Enterococcus_faecium", "gyrA", 83, 'S', 'R'},
    {"Enterococcus_faecium", "gyrA", 83, 'S', 'I'},
    {"Enterococcus_faecium", "gyrA", 87, 'E', 'G'},
    {"Enterococcus_faecium", "gyrA", 87, 'S', 'Y'},
    {"Enterococcus_faecium", "parC", 80, 'E', 'K'},
    {"Enterococcus_faecium", "parC", 80, 'S', 'N'},
    {"Enterococcus_faecium", "parC", 84, 'E', 'K'},
    {"Enterococcus_faecium", "parC", 84, 'E', 'A'},
    {"Escherichia_coli", "emrR", 64, 'I', 'L'},
    {"Escherichia_coli", "emrR", 113, 'A', 'T'},
    {"Escherichia_coli", "gyrA", 81, 'I', 'T'},
    {"Escherichia_coli", "gyrA", 81, 'S', 'I'},
    {"Escherichia_coli", "gyrA", 83, 'A', 'T'},
    {"Escherichia_coli", "gyrA", 83, 'D', 'E'},
    {"Escherichia_coli", "gyrA", 83, 'S', 'A'},
    {"Escherichia_coli", "gyrA", 84, 'L', 'F'},
    {"Escherichia_coli", "gyrA", 87, 'E', 'G'},
    {"Escherichia_coli", "gyrA", 87, 'S', 'T'},
    {"Escherichia_coli", "gyrA", 87, 'E', 'V'},
    {"Escherichia_coli", "gyrA", 106, 'I', 'F'},
    {"Escherichia_coli", "gyrA", 196, 'S', 'R'},
    {"Escherichia_coli", "gyrB", 426, 'E', 'K'},
    {"Escherichia_coli", "parC", 56, 'S', 'T'},
    {"Escherichia_coli", "parC", 57, 'E', 'D'},
    {"Escherichia_coli", "parC", 78, 'A', 'V'},
    {"Escherichia_coli", "parC", 78, 'G', 'D'},
    {"Escherichia_coli", "parC", 80, 'E', 'A'},
    {"Escherichia_coli", "parC", 80, 'L', 'H'},
    {"Escherichia_coli", "parC", 80, 'D', 'N'},
    {"Escherichia_coli", "parC", 80, 'L', 'P'},
    {"Escherichia_coli", "parC", 84, 'S', 'Y'},
    {"Escherichia_coli", "parC", 84, 'D', 'N'},
    {"Escherichia_coli", "parC", 84, 'G', 'C'},
    {"Escherichia_coli", "parC", 108, 'I', 'F'},
    {"Escherichia_coli", "parC", 108, 'S', 'W'},
    {"Escherichia_coli", "parE", 355, 'S', 'V'},
    {"Escherichia_coli", "parE", 416, 'A', 'T'},
    {"Escherichia_coli", "parE", 444, 'L', 'R'},
    {"Escherichia_coli", "parE", 445, 'D', 'N'},
    {"Escherichia_coli", "parE", 458, 'S', 'A'},
    {"Escherichia_coli", "parE", 458, 'D', 'Y'},
    {"Escherichia_coli", "parE", 460, 'A', 'P'},
    {"Escherichia_coli", "parE", 460, 'S', 'W'},
    {"Escherichia_coli", "parE", 464, 'G', 'D'},
    {"Escherichia_coli", "parE", 475, 'D', 'V'},
    {"Escherichia_coli", "parE", 476, 'Q', 'H'},
    {"Escherichia_coli", "parE", 529, 'G', 'C'},
    {"Escherichia_fergusonii", "gyrA", 83, 'S', 'L'},
    {"Haemophilus_influenzae", "gyrA", 84, 'D', 'N'},
    {"Haemophilus_influenzae", "gyrA", 84, 'D', 'G'},
    {"Haemophilus_influenzae", "gyrA", 88, 'D', 'Y'},
    {"Haemophilus_influenzae", "gyrA", 88, 'S', 'I'},
    {"Haemophilus_influenzae", "gyrA", 88, 'S', 'L'},
    {"Haemophilus_influenzae", "parC", 82, 'D', 'N'},
    {"Haemophilus_influenzae", "parC", 84, 'S', 'Y'},
    {"Haemophilus_influenzae", "parC", 84, 'L', 'F'},
    {"Haemophilus_influenzae", "parC", 88, 'S', 'R'},
    {"Haemophilus_influenzae", "parE", 420, 'P', 'S'},
    {"Haemophilus_influenzae", "parE", 439, 'G', 'C'},
    {"Haemophilus_influenzae", "parE", 502, 'E', 'K'},
    {"Klebsiella_aerogenes", "gyrA", 83, 'S', 'T'},
    {"Klebsiella_michiganensis", "gyrA", 83, 'T', 'I'},
    {"Klebsiella_oxytoca", "parC", 80, 'S', 'R'},
    {"Klebsiella_oxytoca", "parC", 80, 'S', 'I'},
    {"Klebsiella_pneumoniae", "gyrA", 83, 'S', 'I'},
    {"Klebsiella_pneumoniae", "gyrA", 83, 'D', 'N'},
    {"Klebsiella_pneumoniae", "gyrA", 83, 'D', 'A'},
    {"Klebsiella_pneumoniae", "gyrA", 83, 'S', 'F'},
    {"Klebsiella_pneumoniae", "gyrA", 87, 'D', 'Y'},
    {"Klebsiella_pneumoniae", "gyrA", 87, 'D', 'G'},
    {"Klebsiella_pneumoniae", "gyrA", 87, 'S', 'Y'},
    {"Klebsiella_pneumoniae", "gyrA", 87, 'D', 'H'},
    {"Klebsiella_pneumoniae", "gyrA", 87, 'S', 'L'},
    {"Klebsiella_pneumoniae", "parC", 80, 'S', 'I'},
    {"Klebsiella_pneumoniae", "parC", 84, 'E', 'K'},
    {"Neisseria_gonorrhoeae", "gyrA", 91, 'D', 'G'},
    {"Neisseria_gonorrhoeae", "gyrA", 91, 'S', 'R'},
    {"Neisseria_gonorrhoeae", "gyrA", 91, 'S', 'F'},
    {"Neisseria_gonorrhoeae", "gyrA", 95, 'D', 'N'},
    {"Neisseria_gonorrhoeae", "gyrA", 95, 'D', 'A'},
    {"Neisseria_gonorrhoeae", "gyrA", 95, 'E', 'G'},
    {"Neisseria_gonorrhoeae", "parC", 86, 'D', 'N'},
    {"Neisseria_gonorrhoeae", "parC", 86, 'S', 'P'},
    {"Neisseria_gonorrhoeae", "parC", 87, 'G', 'V'},
    {"Neisseria_gonorrhoeae", "parC", 87, 'S', 'I'},
    {"Neisseria_gonorrhoeae", "parC", 87, 'S', 'N'},
    {"Neisseria_gonorrhoeae", "parC", 87, 'S', 'Y'},
    {"Neisseria_gonorrhoeae", "parC", 88, 'E', 'K'},
    {"Neisseria_gonorrhoeae", "parC", 91, 'E', 'Q'},
    {"Neisseria_gonorrhoeae", "parC", 91, 'S', 'C'},
    {"Neisseria_gonorrhoeae", "parC", 91, 'S', 'I'},
    {"Neisseria_gonorrhoeae", "parC", 91, 'E', 'A'},
    {"Neisseria_gonorrhoeae", "parE", 410, 'D', 'G'},
    {"Neisseria_meningitidis", "gyrA", 91, 'T', 'I'},
    {"Neisseria_meningitidis", "gyrA", 95, 'T', 'A'},
    {"Neisseria_meningitidis", "gyrA", 173, 'D', 'N'},
    {"Pseudomonas_aeruginosa", "gyrA", 72, 'S', 'L'},
    {"Pseudomonas_aeruginosa", "gyrA", 81, 'T', 'I'},
    {"Pseudomonas_aeruginosa", "gyrA", 81, 'S', 'F'},
    {"Pseudomonas_aeruginosa", "gyrA", 83, 'D', 'Y'},
    {"Pseudomonas_aeruginosa", "gyrA", 83, 'A', 'V'},
    {"Pseudomonas_aeruginosa", "gyrA", 86, 'S', 'W'},
    {"Pseudomonas_aeruginosa", "gyrA", 87, 'D', 'N'},
    {"Pseudomonas_aeruginosa", "gyrA", 87, 'S', 'Y'},
    {"Pseudomonas_aeruginosa", "gyrA", 87, 'D', 'G'},
    {"Pseudomonas_aeruginosa", "gyrA", 87, 'S', 'G'},
    {"Pseudomonas_aeruginosa", "gyrA", 106, 'E', 'D'},
    {"Pseudomonas_aeruginosa", "gyrB", 466, 'S', 'T'},
    {"Pseudomonas_aeruginosa", "gyrB", 466, 'D', 'H'},
    {"Pseudomonas_aeruginosa", "gyrB", 468, 'T', 'A'},
    {"Pseudomonas_aeruginosa", "nfxB", 39, 'E', 'D'},
    {"Pseudomonas_aeruginosa", "parC", 87, 'G', 'C'},
    {"Pseudomonas_aeruginosa", "parC", 87, 'Q', 'L'},
    {"Pseudomonas_aeruginosa", "parE", 457, 'D', 'G'},
    {"Pseudomonas_aeruginosa", "parE", 457, 'G', 'D'},
    {"Pseudomonas_aeruginosa", "parE", 459, 'T', 'P'},
    {"Pseudomonas_aeruginosa", "parE", 473, 'Y', 'N'},
    {"Raoultella_ornithinolytica", "gyrA", 87, 'D', 'G'},
    {"Salmonella_enterica", "gyrA", 72, 'D', 'N'},
    {"Salmonella_enterica", "gyrA", 81, 'S', 'Y'},
    {"Salmonella_enterica", "gyrA", 81, 'S', 'I'},
    {"Salmonella_enterica", "gyrA", 82, 'D', 'Y'},
    {"Salmonella_enterica", "gyrA", 82, 'S', 'R'},
    {"Salmonella_enterica", "gyrA", 83, 'D', 'N'},
    {"Salmonella_enterica", "gyrA", 83, 'G', 'D'},
    {"Salmonella_enterica", "gyrA", 83, 'A', 'E'},
    {"Salmonella_enterica", "gyrA", 87, 'E', 'G'},
    {"Salmonella_enterica", "gyrA", 87, 'G', 'C'},
    {"Salmonella_enterica", "gyrA", 119, 'S', 'Y'},
    {"Salmonella_enterica", "gyrB", 464, 'E', 'D'},
    {"Salmonella_enterica", "gyrB", 464, 'E', 'K'},
    {"Salmonella_enterica", "gyrB", 464, 'S', 'F'},
    {"Salmonella_enterica", "gyrB", 466, 'V', 'G'},
    {"Salmonella_enterica", "parC", 78, 'H', 'Y'},
    {"Salmonella_enterica", "parC", 80, 'S', 'I'},
    {"Salmonella_enterica", "parC", 80, 'S', 'P'},
    {"Salmonella_enterica", "parC", 84, 'G', 'D'},
    {"Salmonella_enterica", "parC", 84, 'S', 'A'},
    {"Salmonella_enterica", "parC", 106, 'W', 'G'},
    {"Salmonella_enterica", "parE", 458, 'D', 'G'},
    {"Salmonella_enterica", "parE", 461, 'A', 'T'},
    {"Salmonella_enterica", "parE", 462, 'D', 'G'},
    {"Salmonella_enterica", "parE", 499, 'S', 'T'},
    {"Serratia_marcescens", "gyrA", 83, 'S', 'I'},
    {"Serratia_marcescens", "gyrA", 83, 'S', 'R'},
    {"Serratia_marcescens", "gyrA", 87, 'D', 'N'},
    {"Shigella_flexneri", "parC", 84, 'E', 'K'},
    {"Shigella_flexneri", "parE", 439, 'P', 'S'},
    {"Shigella_flexneri", "parE", 458, 'S', 'L'},
    {"Staphylococcus_aureus", "gyrA", 84, 'S', 'L'},
    {"Staphylococcus_aureus", "gyrA", 84, 'S', 'Y'},
    {"Staphylococcus_aureus", "gyrA", 84, 'G', 'D'},
    {"Staphylococcus_aureus", "gyrA", 85, 'E', 'G'},
    {"Staphylococcus_aureus", "gyrA", 88, 'S', 'F'},
    {"Staphylococcus_aureus", "gyrA", 88, 'D', 'H'},
    {"Staphylococcus_aureus", "gyrA", 88, 'P', 'S'},
    {"Staphylococcus_aureus", "gyrA", 106, 'E', 'K'},
    {"Staphylococcus_aureus", "gyrB", 436, 'P', 'S'},
    {"Staphylococcus_aureus", "parC", 80, 'E', 'K'},
    {"Staphylococcus_aureus", "parC", 80, 'R', 'S'},
    {"Staphylococcus_aureus", "parC", 83, 'E', 'V'},
    {"Staphylococcus_aureus", "parC", 84, 'D', 'V'},
    {"Staphylococcus_aureus", "parC", 84, 'E', 'A'},
    {"Staphylococcus_aureus", "parC", 84, 'D', 'N'},
    {"Staphylococcus_aureus", "parC", 84, 'S', 'A'},
    {"Staphylococcus_aureus", "parC", 116, 'D', 'N'},
    {"Staphylococcus_aureus", "parC", 116, 'E', 'G'},
    {"Staphylococcus_aureus", "parE", 432, 'S', 'V'},
    {"Staphylococcus_aureus", "parE", 432, 'D', 'E'},
    {"Staphylococcus_aureus", "parE", 432, 'S', 'P'},
    {"Staphylococcus_aureus", "parE", 443, 'A', 'V'},
    {"Staphylococcus_aureus", "parE", 444, 'A', 'E'},
    {"Staphylococcus_aureus", "parE", 451, 'E', 'L'},
    {"Staphylococcus_aureus", "parE", 470, 'Y', 'N'},
    {"Staphylococcus_aureus", "parE", 585, 'N', 'S'},
    {"Staphylococcus_pseudintermedius", "grlA", 80, 'S', 'L'},
    {"Staphylococcus_pseudintermedius", "grlA", 80, 'S', 'I'},
    {"Staphylococcus_pseudintermedius", "grlA", 84, 'D', 'N'},
    {"Staphylococcus_pseudintermedius", "grlA", 84, 'S', 'R'},
    {"Staphylococcus_pseudintermedius", "grlA", 84, 'D', 'G'},
    {"Staphylococcus_pseudintermedius", "grlA", 84, 'D', 'Y'},
    {"Staphylococcus_pseudintermedius", "gyrA", 84, 'S', 'A'},
    {"Staphylococcus_pseudintermedius", "gyrA", 84, 'E', 'G'},
    {"Staphylococcus_pseudintermedius", "gyrA", 84, 'D', 'H'},
    {"Staphylococcus_pseudintermedius", "gyrA", 85, 'E', 'K'},
    {"Staphylococcus_pseudintermedius", "gyrA", 88, 'S', 'W'},
    {"Staphylococcus_pseudintermedius", "gyrA", 88, 'S', 'A'},
    {"Streptococcus_pneumoniae", "parC", 79, 'S', 'Y'},
    {"Streptococcus_pneumoniae", "parC", 79, 'S', 'F'},
    {"Streptococcus_pneumoniae", "parC", 83, 'D', 'H'},
    {"Vibrio_cholerae", "gyrA", 83, 'S', 'I'},
    {"Vibrio_cholerae", "gyrA", 87, 'S', 'L'},
    {"Vibrio_cholerae", "parC", 85, 'D', 'N'},
    {"Vibrio_cholerae", "parE", 420, 'P', 'S'},
    {"Vibrio_cholerae", "parE", 439, 'D', 'N'},
    {"Vibrio_parahaemolyticus", "gyrA", 83, 'S', 'I'},
    {"Vibrio_parahaemolyticus", "parC", 85, 'S', 'F'},
    {"Vibrio_parahaemolyticus", "parC", 85, 'S', 'L'},
    {"Vibrio_vulnificus", "gyrA", 83, 'S', 'I'},
    {"Vibrio_vulnificus", "gyrA", 83, 'S', 'L'},
    {"Vibrio_vulnificus", "parC", 85, 'S', 'R'}
};

// Total mutations: 255

// Summary by species:
//   Acinetobacter_baumannii: 6
//   Burkholderia_cenocepacia: 1
//   Burkholderia_dolosa: 1
//   Campylobacter_coli: 2
//   Campylobacter_jejuni: 3
//   Campylobacter_lari: 1
//   Campylobacter_upsaliensis: 1
//   Citrobacter_farmeri: 1
//   Citrobacter_freundii: 1
//   Clostridioides_difficile: 18
//   Enterobacter_cloacae: 7
//   Enterococcus_faecalis: 10
//   Enterococcus_faecium: 10
//   Escherichia_coli: 39
//   Escherichia_fergusonii: 1
//   Haemophilus_influenzae: 12
//   Klebsiella_aerogenes: 1
//   Klebsiella_michiganensis: 1
//   Klebsiella_oxytoca: 2
//   Klebsiella_pneumoniae: 11
//   Neisseria_gonorrhoeae: 18
//   Neisseria_meningitidis: 3
//   Pseudomonas_aeruginosa: 21
//   Raoultella_ornithinolytica: 1
//   Salmonella_enterica: 25
//   Serratia_marcescens: 3
//   Shigella_flexneri: 3
//   Staphylococcus_aureus: 26
//   Staphylococcus_pseudintermedius: 12
//   Streptococcus_pneumoniae: 3
//   Vibrio_cholerae: 5
//   Vibrio_parahaemolyticus: 3
//   Vibrio_vulnificus: 3

// Summary by gene:
//   emrR: 2
//   grlA: 6
//   gyrA: 121
//   gyrB: 19
//   nfxB: 1
//   parC: 70
//   parE: 36

#endif // FQ_MUTATIONS_HARDCODED_H
