#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdint>
#include <cstring>
#include <dirent.h>
#include <sched.h>
#include <unistd.h>

int g_debug = 0;
static std::vector<std::vector<int>> g_bank_functions;

// Read bank bit mapping functions from file
void read_bank_map_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open map file %s\n", filename);
        exit(1);
    }

    char line[256];
    g_bank_functions.clear();

    while (fgets(line, sizeof(line), fp)) {
        // Skip empty lines and comments
        if (line[0] == '\n' || line[0] == '#') continue;
        std::vector<int> function_bits;
        char* token = strtok(line, " \t\n");
        while (token != nullptr) {
            int bit = atoi(token);
            function_bits.push_back(bit);
            token = strtok(nullptr, " \t\n");
        }
        if (!function_bits.empty()) {
            g_bank_functions.push_back(function_bits);
        }
    }
    fclose(fp);

    if (g_debug) {
        printf("Loaded %zu bank mapping functions:\n", g_bank_functions.size());
        for (size_t i = 0; i < g_bank_functions.size(); i++) {
            printf("Function %zu: XOR bits ", i);
            for (int bit : g_bank_functions[i]) {
                printf("%d ", bit);
            }
            printf("\n");
        }
    }
}

int paddr_to_color(unsigned long mask, unsigned long paddr)
{
    int color = 0;
    for (size_t func_idx = 0; func_idx < g_bank_functions.size(); func_idx++) {
        int bit_result = 0;
        
        // XOR all the specified bits for this function
        for (int bit_pos : g_bank_functions[func_idx]) {
            bit_result ^= ((paddr >> bit_pos) & 0x1);
        }
        
        // Set the corresponding bit in the color
        if (bit_result) {
            color |= (1 << func_idx);
        }
    }
    return color;
}

// Append L1D/L2/L3 indices using the current CPU's cache geometry.
// Assumes contiguous physical-address set bits; hashing/slices are not modeled.
void print_cache_set_indices(uint64_t address) {
    int cpu = sched_getcpu();
    if (cpu < 0) cpu = 0;
    const std::string cache_path = "/sys/devices/system/cpu/cpu" +
                                   std::to_string(cpu) + "/cache";
    uint64_t line_sizes[4] = {};
    uint64_t set_counts[4] = {};
    DIR* cache_dir = opendir(cache_path.c_str());
    if (cache_dir) {
        struct dirent* entry;
        while ((entry = readdir(cache_dir)) != nullptr) {
            if (strncmp(entry->d_name, "index", 5) != 0) continue;
            const std::string path = cache_path + "/" + entry->d_name + "/";
            unsigned int level;
            std::string type;
            uint64_t line_size, num_sets;
            if (!(std::ifstream(path + "level") >> level) || level < 1 || level > 3 ||
                !(std::ifstream(path + "type") >> type) ||
                (type != "Data" && type != "Unified")) continue;
            if (!(std::ifstream(path + "coherency_line_size") >> line_size) ||
                !(std::ifstream(path + "number_of_sets") >> num_sets)) continue;
            // Contiguous set bits require power-of-two geometry.
            if (line_size == 0 || (line_size & (line_size - 1)) != 0 ||
                num_sets == 0 || (num_sets & (num_sets - 1)) != 0) continue;
            line_sizes[level] = line_size;
            set_counts[level] = num_sets;
        }
        closedir(cache_dir);
    }

    for (int level = 1; level <= 3; ++level) {
        std::cout << ' ';
        if (set_counts[level] != 0) {
            std::cout << (address / line_sizes[level]) % set_counts[level];
        } else {
            std::cout << "N/A";
        }
    }
}

int main(int argc, char *argv[]) {
    bool print_cache = false;
    int option;
    while ((option = getopt(argc, argv, "c")) != -1) {
        if (option == 'c') {
            print_cache = true;
        } else {
            std::cerr << "Usage: " << argv[0]
                      << " [-c] <map file> <physical_address>\n"
                      << "  -c  Include L1D, L2, and L3 set indices\n";
            return 1;
        }
    }
    if (argc - optind != 2) {
        std::cerr << "Usage: " << argv[0]
                  << " [-c] <map file> <physical_address>\n"
                  << "  -c  Include L1D, L2, and L3 set indices\n";
        return 1;
    }
    uint64_t address = std::stoull(argv[optind + 1], nullptr, 0);

    read_bank_map_file(argv[optind]);

    int bank_num = 0;

    bank_num = paddr_to_color(0, address);
    std::cout << "0x" << std::hex << address << ' ' << std::dec << bank_num;
    if (print_cache) print_cache_set_indices(address);
    std::cout << std::endl;
    
    return 0;
}
