#ifndef ERF_INPUTS_NAME_H_
#define ERF_INPUTS_NAME_H_

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

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

namespace erf_inputs_detail
{
    struct EntryInfo {
        std::string file;
        int line_number;
        std::string value;
    };

    /**
     * Scan one inputs file, and every file it includes with FILE = <name> (opened by the name
     * as written, as ParmParse does), recording each key with where it was first set. Returns
     * whether a key was set twice anywhere in the deck. A file that includes itself, directly or
     * through others, aborts, as does an included file that cannot be opened.
     */
    inline bool
    scan_inputs (const std::string& inputs_file, std::unordered_map<std::string, EntryInfo>& seen,
                 std::vector<std::string>& open_files)
    {
        std::ifstream ifs(inputs_file);

        if (!ifs.is_open()) {
            amrex::Abort("Could not open inputs file: " + inputs_file);
        }
        for (const auto& f : open_files) {
            if (f == inputs_file) {
                amrex::Abort("The inputs file " + inputs_file + " includes itself (FILE = ...)");
            }
        }
        open_files.push_back(inputs_file);

        bool found_duplicates = false;
        std::string line;
        int line_number = 0;

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

            // FILE = <name> includes another inputs file: its keys belong to the same deck, so a
            // key set there and here, or in two included files, is a duplicate as well
            if (key == amrex::ParmParse::FileKeyword) {
                if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
                    value = value.substr(1, value.size() - 2);
                }
                found_duplicates = scan_inputs(value, seen, open_files) || found_duplicates;
                continue;
            }

            auto it = seen.find(key);

            if (it != seen.end())
            {
                found_duplicates = true;

                amrex::Print()
                    << "\nDuplicate input detected:\n"
                    << "    Key         : " << key << "\n"
                    << "    First line  : " << it->second.line_number << " of " << it->second.file
                    << "    Value       : " << it->second.value << "\n"
                    << "    Second line : " << line_number << " of " << inputs_file
                    << "    Value       : " << value << "\n";
            }
            else
            {
                seen.emplace(key, EntryInfo{inputs_file, line_number, value});
            }
        }

        open_files.pop_back();
        return found_duplicates;
    }
}

/**
 * Abort if any key is set more than once in the inputs file, counting the files it includes with
 * FILE = <name>: a deck may be split into included files, but every key is set in one place.
 * Keys given on the command line are not checked; they override the file's.
 */
inline void
CheckForDuplicateInputs (const std::string& inputs_file)
{
    std::unordered_map<std::string, erf_inputs_detail::EntryInfo> seen;
    std::vector<std::string> open_files;

    if (erf_inputs_detail::scan_inputs(inputs_file, seen, open_files))
    {
        amrex::Abort("Duplicate inputs detected in inputs file");
    }
}

#endif
