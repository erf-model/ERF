#ifndef ERF_INPUTS_NAME_H_
#define ERF_INPUTS_NAME_H_

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <regex>
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
        // the include instances open when the key was set, outermost first
        std::vector<int> includes;
    };

    struct ScanState {
        std::unordered_map<std::string, EntryInfo> seen;
        std::vector<int> includes;                     // include instances now open
        std::vector<std::filesystem::path> open_files; // their files, for the loop check
        int next_include = 0;
        bool found_duplicates = false;
        std::string error;                             // set when the scan cannot go on
    };

    /**
     * The name ParmParse opens for an inputs file: $AMREX_INPUTS_FILE_PREFIX, with a trailing
     * '/' added when missing, in front of the name as written (AMReX_ParmParse.cpp, read_file).
     */
    inline std::string
    prefixed_name (const std::string& name)
    {
        const char* prefix = std::getenv("AMREX_INPUTS_FILE_PREFIX");
        if (prefix == nullptr || prefix[0] == 0) { return name; }
        std::string p(prefix);
        if (p.back() != '/') { p += '/'; }
        return p + name;
    }

    /**
     * The values on the right of '=': whitespace separates them, and a double-quoted value is one
     * value without its quotes.
     */
    inline std::vector<std::string>
    split_values (const std::string& rhs)
    {
        std::vector<std::string> vals;
        std::size_t i = 0;
        while (i < rhs.size()) {
            if (rhs[i] == ' ' || rhs[i] == '\t' || rhs[i] == '\r') { ++i; continue; }
            if (rhs[i] == '"') {
                auto e = rhs.find('"', i + 1);
                if (e == std::string::npos) { e = rhs.size(); }
                vals.push_back(rhs.substr(i + 1, e - i - 1));
                i = e + 1;
            } else {
                auto e = rhs.find_first_of(" \t\r", i);
                if (e == std::string::npos) { e = rhs.size(); }
                vals.push_back(rhs.substr(i, e - i));
                i = e;
            }
        }
        return vals;
    }

    /**
     * Whether a key set again where `includes` are open overrides the setting `first`: it does
     * when the first setting came from a file this one included earlier (the base + override
     * idiom; ParmParse keeps the last value). A key set twice in one file, in two files included
     * side by side, or in a file and then again in a file it includes afterwards is a duplicate.
     */
    inline bool
    is_override (const EntryInfo& first, const std::vector<int>& includes)
    {
        if (first.includes.size() <= includes.size()) { return false; }
        for (std::size_t n = 0; n < includes.size(); ++n) {
            if (first.includes[n] != includes[n]) { return false; }
        }
        return true;
    }

    /**
     * Scan one inputs file, and every file it includes with FILE = <name>, the way ParmParse reads
     * them: names are opened with $AMREX_INPUTS_FILE_PREFIX in front, FILE with more than one value
     * is a parameter, UNSET = <keys> drops the keys, [name] starts a table whose keys are
     * name.<key>, and #if AMREX_SPACEDIM / #ifdef AMREX_USE_GPU regions that are off are skipped.
     * Records each key with where it was set, and sets found_duplicates when a key is set twice
     * (see is_override). An inputs file that cannot be opened, or that includes itself directly or
     * through others, sets error and stops the scan.
     */
    inline void
    scan_inputs (const std::string& inputs_file, ScanState& state)
    {
        const std::string name = prefixed_name(inputs_file);
        std::ifstream ifs(name);

        if (!ifs.is_open()) {
            state.error = "Could not open inputs file: " + name;
            return;
        }
        std::error_code ec;
        auto canonical = std::filesystem::weakly_canonical(name, ec);
        if (ec) { canonical = std::filesystem::absolute(name); }
        for (const auto& f : state.open_files) {
            if (f == canonical) {
                state.error = "The inputs file " + name + " includes itself (FILE = ...)";
                return;
            }
        }
        state.open_files.push_back(canonical);
        state.includes.push_back(state.next_include++);

        // the preprocessor ParmParse applies before it reads the lines (AMReX_ParmParse.cpp, read_file)
        static const std::regex if_regex(R"(^\s*#\s*if\s+\(?\s*AMREX_SPACEDIM\s*(>|<|==|>=|<=|!=)\s*([1-3])\s*\)?\s*$)");
        static const std::regex elif_regex(R"(^\s*#\s*elif\s+\(?\s*AMREX_SPACEDIM\s*(>|<|==|>=|<=|!=)\s*([1-3])\s*\)?\s*$)");
        static const std::regex else_regex(R"(^\s*#\s*else\s*$)");
        static const std::regex endif_regex(R"(^\s*#\s*endif\s*$)");
        static const std::regex if_gpu_regex(R"(^\s*#\s*ifdef\s+AMREX_USE_GPU\s*$)");
        static const std::regex if_not_gpu_regex(R"(^\s*#\s*ifndef\s+AMREX_USE_GPU\s*$)");
        static const std::regex table_regex(R"(^\[\s*([A-Za-z][A-Za-z0-9_.\-]*)\s*\]$)");
#ifdef AMREX_USE_GPU
        constexpr bool use_gpu = true;
#else
        constexpr bool use_gpu = false;
#endif
        auto dim_true = [] (const std::smatch& sm) {
            const std::string op = sm[1].str();
            const int dim = std::stoi(sm[2].str());
            return (op == "<"  && AMREX_SPACEDIM <  dim) || (op == ">"  && AMREX_SPACEDIM >  dim) ||
                   (op == "==" && AMREX_SPACEDIM == dim) || (op == "<=" && AMREX_SPACEDIM <= dim) ||
                   (op == ">=" && AMREX_SPACEDIM >= dim) || (op == "!=" && AMREX_SPACEDIM != dim);
        };
        std::vector<bool> region_on;  // is each open #if region read?
        std::vector<bool> region_had; // has a branch of each open #if region been read?

        std::string table; // the current [name] table, empty at the root
        std::string line;
        int line_number = 0;

        while (std::getline(ifs, line))
        {
            ++line_number;

            std::smatch sm;
            if (std::regex_match(line, if_gpu_regex) || std::regex_match(line, if_not_gpu_regex)) {
                const bool r = std::regex_match(line, if_gpu_regex) ? use_gpu : !use_gpu;
                region_on.push_back(r); region_had.push_back(r);
                continue;
            } else if (std::regex_match(line, sm, if_regex)) {
                const bool r = dim_true(sm);
                region_on.push_back(r); region_had.push_back(r);
                continue;
            } else if (!region_on.empty() && std::regex_match(line, sm, elif_regex)) {
                const bool r = !region_had.back() && dim_true(sm);
                region_on.back() = r;
                if (r) { region_had.back() = true; }
                continue;
            } else if (!region_on.empty() && std::regex_match(line, else_regex)) {
                region_on.back() = !region_had.back();
                region_had.back() = true;
                continue;
            } else if (!region_on.empty() && std::regex_match(line, endif_regex)) {
                region_on.pop_back(); region_had.pop_back();
                continue;
            }
            bool skip = false;
            for (bool r : region_on) { if (!r) { skip = true; } }
            if (skip) { continue; }

            // Remove comments
            auto hash_pos = line.find('#');
            if (hash_pos != std::string::npos) {
                line = line.substr(0, hash_pos);
            }

            line = trim_copy(line);

            // Skip blank lines
            if (line.empty()) continue;

            // A [name] line starts a table: the keys after it are name.<key>
            if (std::regex_match(line, sm, table_regex)) {
                table = sm[1].str();
                continue;
            }

            // Skip lines that are not assignments
            auto eq_pos = line.find('=');
            if (eq_pos == std::string::npos) continue;

            const std::string def   = trim_copy(line.substr(0, eq_pos));
            const std::string value = trim_copy(line.substr(eq_pos + 1));
            const auto vals = split_values(value);

            // FILE = <name> includes another inputs file, whose keys belong to the same deck.
            // ParmParse reads it into a clean table scope.
            if (def == amrex::ParmParse::FileKeyword && vals.size() == 1) {
                scan_inputs(vals.front(), state);
                // stop, closing this file's include on the way out
                if (!state.error.empty()) { break; }
                continue;
            }

            // UNSET = <keys> drops the keys (full names, also inside a table), so they may be set again
            if (def == amrex::ParmParse::UnsetKeyword) {
                for (const auto& k : vals) { state.seen.erase(k); }
                continue;
            }

            const std::string key = table.empty() ? def : table + "." + def;
            auto it = state.seen.find(key);

            if (it == state.seen.end())
            {
                state.seen.emplace(key, EntryInfo{name, line_number, value, state.includes});
            }
            else if (is_override(it->second, state.includes))
            {
                it->second = EntryInfo{name, line_number, value, state.includes};
            }
            else
            {
                state.found_duplicates = true;

                amrex::Print()
                    << "\nDuplicate input detected:\n"
                    << "    Key         : " << key << "\n"
                    << "    First line  : " << it->second.line_number << " of " << it->second.file << "\n"
                    << "    Value       : " << it->second.value << "\n"
                    << "    Second line : " << line_number << " of " << name << "\n"
                    << "    Value       : " << value << "\n";
            }
        }

        state.includes.pop_back();
        state.open_files.pop_back();
    }
}

/**
 * Abort if a key is set twice in the inputs file, counting the files it includes with
 * FILE = <name>. A file may re-set a key that a file it included earlier set (a base deck
 * included first and then overridden); a key set twice in one file, in two files included side by
 * side, or in a file and then again in a file it includes afterwards is refused. UNSET = <key>
 * drops a key so it may be set again. Keys given on the command line are not checked; they
 * override the file's.
 */
inline void
CheckForDuplicateInputs (const std::string& inputs_file)
{
    erf_inputs_detail::ScanState state;
    erf_inputs_detail::scan_inputs(inputs_file, state);

    if (!state.error.empty()) {
        amrex::Abort(state.error);
    }
    if (state.found_duplicates)
    {
        amrex::Abort("Duplicate inputs detected in inputs file");
    }
}

#endif
