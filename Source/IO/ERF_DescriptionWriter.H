#ifndef ERF_DESCRIPTION_WRITER_H_
#define ERF_DESCRIPTION_WRITER_H_

#include <AMReX_BLassert.H>

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace erf_description {

enum class Format { JSON, YAML };

/** Fixed-structure JSON/YAML writer with no parser dependency. YAML output uses block
 * mappings, JSON flow arrays and double-quoted keys and strings, so YAML 1.1 parsers
 * never read them as booleans or numbers.
 */
class Writer
{
public:
    Writer (std::ostream& output, Format format) : os(output), yaml(format == Format::YAML)
    {
        if (!yaml) os << '{';
    }

    void begin_object (const char* name)
    {
        key(name, true);
        os << (yaml ? '\n' : '{');
        first.push_back(true);
    }

    void end_object ()
    {
        AMREX_ASSERT(first.size() > 1);
        first.pop_back();
        if (!yaml) os << '}';
    }

    void finish (bool newline = true)
    {
        AMREX_ASSERT(first.size() == 1);
        if (!yaml) { os << '}'; if (newline) os << '\n'; }
    }

    void value (const char* name, std::nullptr_t)
    {
        key(name); os << "null"; end_value();
    }

    void value (const char* name, bool v)
    {
        key(name); os << (v ? "true" : "false"); end_value();
    }

    void value (const char* name, int v)
    {
        key(name); os << std::to_string(v); end_value();
    }

    void value (const char* name, const char* v)
    {
        key(name);
        if (v == nullptr) os << "null";
        else quoted(v);
        end_value();
    }

    // AMReX build-info getters return "" for unknown values.
    void optional_string (const char* name, const char* v)
    {
        value(name, v != nullptr && *v != '\0' ? v : nullptr);
    }

    void strings (const char* name, const std::vector<std::string>& values)
    {
        key(name);
        os << '[';
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i) os << ", ";
            quoted(values[i].c_str());
        }
        os << ']';
        end_value();
    }

    template <class T, class Emit>
    void objects (const char* name, const std::vector<T>& values, Emit emit)
    {
        key(name);
        os << '[';
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i) os << ", ";
            Writer item(os, Format::JSON);
            emit(item, values[i]);
            item.finish(false);
        }
        os << ']';
        end_value();
    }

private:
    std::ostream& os;
    bool yaml;
    std::vector<bool> first{true};

    void key (const char* name, bool object = false)
    {
        if (yaml) {
            os << std::string((first.size() - 1) * 2, ' ');
            quoted(name);
            os << ':';
            if (!object) os << ' ';
        } else {
            if (!first.back()) os << ',';
            quoted(name);
            os << ':';
        }
        first.back() = false;
    }

    void end_value () { if (yaml) os << '\n'; }

    void quoted (const char* value)
    {
        constexpr char hex[] = "0123456789abcdef";
        os << '"';
        for (const auto* c = reinterpret_cast<const unsigned char*>(value); *c; ++c) {
            if (*c == '"' || *c == '\\') {
                os << '\\' << static_cast<char>(*c);
            } else if (*c < 0x20 || *c == 0x7f) {
                os << "\\u00" << hex[*c >> 4] << hex[*c & 15];
            } else if (*c == 0xc2 && c[1] >= 0x80 && c[1] <= 0x9f) {
                // C1 controls (incl. NEL, a YAML line break)
                ++c;
                os << "\\u00" << hex[*c >> 4] << hex[*c & 15];
            } else if (*c == 0xe2 && c[1] == 0x80 && (c[2] == 0xa8 || c[2] == 0xa9)) {
                // U+2028/9, line breaks to YAML 1.1 parsers
                os << "\\u202" << (c[2] == 0xa8 ? '8' : '9');
                c += 2;
            } else {
                os << static_cast<char>(*c);
            }
        }
        os << '"';
    }
};

} // namespace erf_description
#endif
