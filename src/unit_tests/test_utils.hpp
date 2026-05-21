#pragma once

#include <iostream>
#include <string>

namespace gene_hmm {

    using namespace std;

    static bool ut_pass(const string& label) { cout << "  [PASS] " << label << "\n"; return true; }
    static bool ut_fail(const string& label) { cout << "  [FAIL] " << label << "\n"; return false; }

    #define CHECK(label, cond) ((cond) ? ut_pass(label) : ut_fail(label))

}
