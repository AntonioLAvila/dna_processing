// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <utility>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include "json.hpp"

// namespace py = pybind11;
using json = nlohmann::json;


void encode(
    const std::set<std::string>& chromosomes,
    const std::string& fasta_file,
    const std::string& output_directory
) {

    std::ifstream file(fasta_file);
    if (!file) {
        std::cerr << "Couldn't open FASTA file: " << fasta_file << std::endl;
        return;
    }

    std::ofstream data_file(output_directory + "/chromosomes.dat", std::ios::binary);
    if (!data_file) {
        std::cerr << "Failed to create chromosomes.dat" << std::endl;
        return;
    }

    std::string line, current_chrom;
    std::ostringstream seq_buffer;
    json index_json;
    size_t current_offset = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Write previous chromosome to .dat if in set
            if (!current_chrom.empty() && chromosomes.contains(current_chrom)) {
                std::string sequence = seq_buffer.str();
                for (char& base : sequence) {
                    base = static_cast<char>(std::toupper(base));
                }

                data_file.write(sequence.data(), sequence.size());

                // Record offset and length in JSON index
                index_json[current_chrom] = { current_offset, sequence.size() };
                current_offset += sequence.size();
            }

            seq_buffer.str("");
            seq_buffer.clear();

            std::istringstream iss(line.substr(1));
            iss >> current_chrom;
        } else {
            if (chromosomes.contains(current_chrom)) {
                seq_buffer << line;
            }
        }
    }

    // Handle last chromosome
    if (!current_chrom.empty() && chromosomes.contains(current_chrom)) {
        std::string sequence = seq_buffer.str();
        for (char& base : sequence) {
            base = static_cast<char>(std::toupper(base));
        }

        data_file.write(sequence.data(), sequence.size());
        index_json[current_chrom] = { current_offset, sequence.size() };
    }

    data_file.close();

    // Write chromosomes.idx
    std::ofstream index_file(output_directory + "/chromosomes.idx");
    if (!index_file) {
        std::cerr << "Failed to create chromosomes.idx" << std::endl;
        return;
    }

    index_file << index_json.dump(2);
}


void mutate(
    const std::string &chromosome,
    const std::vector<std::pair<int, char>>& mutations,
    const std::string& chromosome_data_path,
    const std::string& output_path
) {
    // Load the JSON index
    std::ifstream index_file(chromosome_data_path + "/chromosomes.idx");
    if (!index_file) {
        std::cerr << "Failed to open index file\n";
        return;
    }

    json index;
    index_file >> index;

    if (!index.contains(chromosome)) {
        std::cerr << "Chromosome '" << chromosome << "' not found in index.\n";
        return;
    }

    size_t offset = index[chromosome][0];
    size_t length = index[chromosome][1];

    std::string dat_path = chromosome_data_path + "/chromosomes.dat";
    int fd = open(dat_path.c_str(), O_RDONLY);
    if (fd == -1) {
        perror("open");
        return;
    }

    char* mapped = static_cast<char*>(mmap(nullptr, offset + length, PROT_READ, MAP_PRIVATE, fd, 0));
    if (mapped == MAP_FAILED) {
        perror("mmap");
        close(fd);
        return;
    }
    close(fd);

    const char* chromosome_seq = mapped + offset;

    std::string mutated_seq(chromosome_seq, length);
    for (const auto& [pos, base] : mutations) {
        if (pos < 0 || static_cast<size_t>(pos) >= length) {
            std::cerr << "Warning: Mutation position " << pos << " is out of bounds for " << chromosome << "\n";
            continue;
        }
        mutated_seq[pos] = std::toupper(base);
    }

    std::ofstream out(output_path);
    if (!out) {
        std::cerr << "Failed to open output file: " << output_path << "\n";
        munmap(mapped, offset + length);
        return;
    }

    // wrap lines at 60 bases like FASTA
    const size_t line_width = 60;
    for (size_t i = 0; i < mutated_seq.size(); i += line_width) {
        out << mutated_seq.substr(i, line_width) << '\n';
    }

    munmap(mapped, offset + length);
}


int main() {

    encode({"1", "X"}, "Homo_sapiens.GRCh38.dna.primary_assembly.fa", "test_out");

    mutate(
        "1",
        { {0, 'G'}, {5, 'A'} },
        "test_out",
        "chr1_mutations.txt"
    );
    
    return 0;
}

// PYBIND11_MODULE(dna_processing, m) {
//     m.doc() = "Library for mutating DNA sequences.";
//     m.def(
//         "encode",
//         &encode,
//         "Stores specific chromosomes from a fasta file in a quickly accessible way.",
//         py::arg("chromosomes"),
//         py::arg("fasta_file"),
//         py::arg("output_directory")
//     );
//     m.def(
//         "mutate",
//         &mutate,
//         "Mutates a chromosome from the parsed data and store it where you specify.",
//         py::arg("chromosome"),
//         py::arg("mutations"),
//         py::arg("chromosome_data_directory"),
//         py::arg("output_filepath")
//     );
// }