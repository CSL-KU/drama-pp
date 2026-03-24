#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdint>
#include <cstring>
#include <algorithm>

int g_debug = 0;

struct RegionMapping {
    struct SelectorTerm {
        int expected = 0;
        std::vector<int> bits;

        bool matches(uint64_t paddr) const {
            int parity = 0;
            for (int bit_pos : bits) {
                parity ^= ((paddr >> bit_pos) & 0x1ULL);
            }
            return parity == expected;
        }
    };

    std::vector<SelectorTerm> selectors;
    std::vector<std::vector<int>> bank_functions;
};

static std::vector<RegionMapping> g_region_mappings;
static int g_max_local_bank_bits = 0;

static int compute_local_color(const RegionMapping& region, uint64_t paddr) {
    int color = 0;
    for (size_t func_idx = 0; func_idx < region.bank_functions.size(); func_idx++) {
        int bit_result = 0;
        for (int bit_pos : region.bank_functions[func_idx]) {
            bit_result ^= ((paddr >> bit_pos) & 0x1ULL);
        }
        if (bit_result) {
            color |= (1 << func_idx);
        }
    }
    return color;
}

// Read bank bit mapping functions from file
void read_bank_map_file(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open map file %s\n", filename);
        exit(1);
    }

    char line[256];
    g_region_mappings.clear();
    g_max_local_bank_bits = 0;
    RegionMapping* current_region = nullptr;

    while (fgets(line, sizeof(line), fp)) {
        // Skip empty lines and comments
        if (line[0] == '\n' || line[0] == '#') continue;

        std::stringstream header_stream(line);
        std::string header_word;
        header_stream >> header_word;
        if (header_word == "region") {
            g_region_mappings.push_back(RegionMapping{});
            current_region = &g_region_mappings.back();

            std::string mask_text;
            std::string value_text;
            if (header_stream >> mask_text >> value_text) {
                uint64_t mask = std::stoull(mask_text, nullptr, 0);
                uint64_t value = std::stoull(value_text, nullptr, 0);
                for (int bit = 0; mask != 0; bit++) {
                    if (mask & 0x1ULL) {
                        RegionMapping::SelectorTerm term;
                        term.expected = (value >> bit) & 0x1ULL;
                        term.bits.push_back(bit);
                        current_region->selectors.push_back(term);
                    }
                    mask >>= 1;
                }
            }
            continue;
        }

        if (header_word == "selector") {
            if (current_region == nullptr) {
                g_region_mappings.push_back(RegionMapping{});
                current_region = &g_region_mappings.back();
            }
            RegionMapping::SelectorTerm term;
            if (!(header_stream >> term.expected)) {
                fprintf(stderr, "Error: malformed selector line in %s\n", filename);
                exit(1);
            }
            int bit = 0;
            bool saw_bit = false;
            while (header_stream >> bit) {
                term.bits.push_back(bit);
                saw_bit = true;
            }
            if (!saw_bit) {
                fprintf(stderr, "Error: selector line without bits in %s\n", filename);
                exit(1);
            }
            current_region->selectors.push_back(term);
            continue;
        }

        if (current_region == nullptr) {
            g_region_mappings.push_back(RegionMapping{});
            current_region = &g_region_mappings.back();
        }

        std::vector<int> function_bits;
        char* token = strtok(line, " \t\n");
        while (token != nullptr) {
            int bit = atoi(token);
            function_bits.push_back(bit);
            token = strtok(nullptr, " \t\n");
        }
        if (!function_bits.empty()) {
            current_region->bank_functions.push_back(function_bits);
        }
    }
    fclose(fp);

    for (const auto& region : g_region_mappings) {
        g_max_local_bank_bits = std::max(g_max_local_bank_bits, static_cast<int>(region.bank_functions.size()));
    }

    if (g_debug) {
        printf("Loaded %zu region mapping(s):\n", g_region_mappings.size());
        for (size_t region_idx = 0; region_idx < g_region_mappings.size(); region_idx++) {
            const auto& region = g_region_mappings[region_idx];
            printf("Region %zu:\n", region_idx);
            for (const auto& selector : region.selectors) {
                printf("  selector %d ", selector.expected);
                for (int bit : selector.bits) {
                    printf("%d ", bit);
                }
                printf("\n");
            }
            for (size_t func_idx = 0; func_idx < region.bank_functions.size(); func_idx++) {
                printf("  Function %zu: XOR bits ", func_idx);
                for (int bit : region.bank_functions[func_idx]) {
                    printf("%d ", bit);
                }
                printf("\n");
            }
        }
    }
}

int paddr_to_color(unsigned long mask, unsigned long paddr)
{
    (void)mask;
    size_t matched_region = static_cast<size_t>(-1);
    for (size_t idx = 0; idx < g_region_mappings.size(); idx++) {
        const auto& region = g_region_mappings[idx];
        bool matched = true;
        for (const auto& selector : region.selectors) {
            if (!selector.matches(paddr)) {
                matched = false;
                break;
            }
        }
        if (matched) {
            if (matched_region != static_cast<size_t>(-1)) {
                fprintf(stderr, "Error: address 0x%lx matches multiple regions\n", paddr);
                exit(1);
            }
            matched_region = idx;
        }
    }

    if (matched_region == static_cast<size_t>(-1)) {
        fprintf(stderr, "Error: address 0x%lx does not match any region\n", paddr);
        exit(1);
    }

    int local_color = compute_local_color(g_region_mappings[matched_region], paddr);
    if (g_region_mappings.size() == 1) {
        return local_color;
    }
    return (static_cast<int>(matched_region) << g_max_local_bank_bits) | local_color;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <map file> <physical_address>" << std::endl;
        return 1;
    }
    uint64_t address = std::stoull(argv[2], nullptr, 0);

    read_bank_map_file(argv[1]);

    int bank_num = 0;

    bank_num = paddr_to_color(0, address);
    std::cout << "Physical address 0x" << std::hex << address << " maps to bank " << std::dec << bank_num << std::endl;
    
    return 0;
}