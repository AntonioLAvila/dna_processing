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


void process_and_write_chromosome(
    const std::string &chrom_name,
    size_t &file_offset,
    size_t chrom_length,
    json &index
) {
    if (chrom_name.empty() || chrom_length == 0) return;
    
    index[chrom_name] = { file_offset, chrom_length };
    
    file_offset += chrom_length;
}


void encode(
    const std::set<std::string> &chromosomes,
    const std::string &fasta_file,
    const std::string &output_directory
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

    json index_json;
    std::string line, current_chrom_name;
    size_t file_offset = 0;
    size_t current_chrom_len = 0;
    bool is_target_chrom = false;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '\r') continue;

        if (line[0] == '>') {
            // Process the prev chromosome
            process_and_write_chromosome(current_chrom_name, file_offset, current_chrom_len, index_json);
            
            // Reset for new chromosome
            current_chrom_len = 0;
            size_t first_space = line.find(' ');
            current_chrom_name = line.substr(1, first_space - 1);

            // Check if this new chromosome is one we want to keep
            is_target_chrom = chromosomes.contains(current_chrom_name);

        } else if (is_target_chrom) {
            // Process the line directly without buffering
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            data_file.write(line.data(), line.size());
            current_chrom_len += line.size();
        }
    }

    // Process last chromosome
    process_and_write_chromosome(current_chrom_name, file_offset, current_chrom_len, index_json);
    data_file.close();

    std::ofstream index_file(output_directory + "/chromosomes.idx");
    if (!index_file) {
        std::cerr << "Failed to create chromosomes.idx" << std::endl;
        return;
    }
    index_file << index_json.dump(2);
}


bool mutate(
    const std::string &chromosome,
    const std::vector<std::tuple<int, std::string, std::string>> &mutations,
    const std::string &chromosome_data_path,
    const std::string &output_path
) {
    std::ifstream index_file(chromosome_data_path + "/chromosomes.idx");
    if (!index_file) { return false; }
    json index;
    index_file >> index;
    if (!index.contains(chromosome)) { return false; }
    size_t offset = index[chromosome][0];
    size_t length = index[chromosome][1];
    std::string dat_path = chromosome_data_path + "/chromosomes.dat";
    int fd = open(dat_path.c_str(), O_RDONLY);
    if (fd == -1) { return false; }
    char* mapped = static_cast<char*>(mmap(nullptr, offset + length, PROT_READ, MAP_PRIVATE, fd, 0));
    if (mapped == MAP_FAILED) { close(fd); return false; }
    close(fd);

    const char* original_seq = mapped + offset;

    std::vector<size_t> mutation_indices(mutations.size());
    std::iota(mutation_indices.begin(), mutation_indices.end(), 0);

    std::sort(mutation_indices.begin(), mutation_indices.end(),
        [&mutations](size_t a_idx, size_t b_idx) {
            return std::get<0>(mutations[a_idx]) < std::get<0>(mutations[b_idx]);
        });

    // Pre-calculating final_size
    size_t final_size = length;
    for (const auto &[pos, type, bases] : mutations) {
        if (type == "INS") final_size += bases.length();
        else if (type == "DEL") final_size--;
    }
    std::string mutated_seq;
    mutated_seq.reserve(final_size);

    size_t last_pos = 0;
    for (size_t index : mutation_indices) {
        const auto& [pos, type, bases] = mutations[index];

        if (pos < 0 || static_cast<size_t>(pos) > length) { /* ... */ continue; }
        
        // Append unchanged chunk
        if (static_cast<size_t>(pos) > last_pos) {
            mutated_seq.append(original_seq + last_pos, pos - last_pos);
        }

        if (type == "SNV") {
            mutated_seq.append(bases);
            last_pos = pos + 1;
        } else if (type == "INS") {
            mutated_seq.append(bases);
            last_pos = pos;
        } else if (type == "DEL") {
            last_pos = pos + bases.length();
        }
    }

    if (last_pos < length) {
        mutated_seq.append(original_seq + last_pos, length - last_pos);
    }
    
    std::ofstream out(output_path);
    if (!out) {
        munmap(mapped, offset + length);
        return false;
    }
    const size_t line_width = 60;
    for (size_t i = 0; i < mutated_seq.size(); i += line_width) {
        out << mutated_seq.substr(i, line_width) << '\n';
    }
    munmap(mapped, offset + length);
    return true;
}


int main() {

    encode({"1", "X"}, "test/Homo_sapiens.GRCh38.dna.primary_assembly.fa", "test/test_out");

    mutate(
        "1",
        { {0, "INS", "A"}, {5, "DEL", ""} },
        "test/test_out",
        "test/chr1_mutations.txt"
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