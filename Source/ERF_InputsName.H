#ifndef ERF_INPUTS_NAME_H_
#define ERF_INPUTS_NAME_H_

#include <string>

extern std::string inputs_name;

namespace
{
    std::string trim_copy (std::string s)
    {
        auto b = s.find_first_not_of(" \t\r");
        if (b == std::string::npos) return "";

        auto e = s.find_last_not_of(" \t\r");
        return s.substr(b, e - b + 1);
    }
}

inline void
CheckForDuplicateInputs (const std::string& inputs_file)
{
    std::ifstream ifs(inputs_file);

    if (!ifs.is_open()) {
        amrex::Abort("Could not open inputs file: " + inputs_file);
    }

    struct EntryInfo {
        int line_number;
        std::string value;
    };

    std::unordered_map<std::string, EntryInfo> seen;

    std::string line;
    int line_number = 0;
    bool found_duplicates = false;

    while (std::getline(ifs, line))
    {
        ++line_number;

        // Remove comments
        auto hash_pos = line.find('#');
        if (hash_pos != std::string::npos) {
            line = line.substr(0, hash_pos);
        }

        line = trim_copy(line);

        // Skip blank lines
        if (line.empty()) continue;

        // Skip lines that are not assignments
        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key   = trim_copy(line.substr(0, eq_pos));
        std::string value = trim_copy(line.substr(eq_pos + 1));

        auto it = seen.find(key);

        if (it != seen.end())
        {
            found_duplicates = true;

            amrex::Print()
                << "\nDuplicate input detected:\n"
                << "    Key         : " << key << "\n"
                << "    First line  : " << it->second.line_number
                << "    Value       : " << it->second.value << "\n"
                << "    Second line : " << line_number
                << "    Value       : " << value << "\n";
        }
        else
        {
            seen.emplace(key, EntryInfo{line_number, value});
        }
    }

    if (found_duplicates)
    {
        amrex::Abort("Duplicate inputs detected in inputs file");
    }
}

#endif
