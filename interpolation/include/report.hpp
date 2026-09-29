// Minimal table printer, CSV writer and check counter (no test framework).
#pragma once
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace report {

inline std::string results_dir() {
    const char* d = std::getenv("RESULTS_DIR");
    std::string dir = d ? d : "results";
    std::filesystem::create_directories(dir);
    return dir;
}

class Csv {
public:
    Csv(const std::string& name, const std::string& header)
        : out_(results_dir() + "/" + name + ".csv") {
        out_ << header << '\n';
    }
    template <class... A> void row(const A&... a) {
        bool first = true;
        ((out_ << (first ? "" : ",") << a, first = false), ...);
        out_ << '\n';
    }

private:
    std::ofstream out_;
};

// Fixed-width table printed to stdout.
class Table {
public:
    explicit Table(std::vector<std::string> head, int width = 11) : w_(width) { line(head); rule(head.size()); }
    void line(const std::vector<std::string>& cells) {
        for (std::size_t i = 0; i < cells.size(); ++i)
            std::printf(i == 0 ? "%-*s" : " %*s", i == 0 ? 10 : w_, cells[i].c_str());
        std::printf("\n");
    }

private:
    void rule(std::size_t n) { std::printf("%s\n", std::string(10 + (n - 1) * (w_ + 1), '-').c_str()); }
    int w_;
};

inline std::string num(double v, const char* fmt = "%.3g") {
    char buf[64];
    std::snprintf(buf, sizeof buf, fmt, v);
    return buf;
}

// Checks of guaranteed properties: failures make the program exit non-zero.
struct Checks {
    int failed = 0;
    void expect(bool ok, const std::string& what) {
        if (!ok) {
            ++failed;
            std::printf("  CHECK FAILED: %s\n", what.c_str());
        }
    }
    int exit_code() const {
        std::printf(failed ? "\n%d guaranteed-property check(s) FAILED\n" : "\nall guaranteed-property checks passed\n", failed);
        return failed ? 1 : 0;
    }
};

}  // namespace report
